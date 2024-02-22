import sys, os
import fcntl
import random
import time
import subprocess
from collections import defaultdict

init_size = 10000
est_time = 600 # 10 minutes
max_block_per_batch = 20
max_transfer_per_block = 100
min_interval_in_sec = 2
max_interval_in_sec = 20
log_max_init_balance = 200.0
log_max_transfer_amount = 180.0
max_addr = 2**160 - 1
prob_new_addr = 0.1
unstable_block_num = 10
unstable_rate = 0.2 # 20% of the transfers in the last 10 blocks will change
client_count = 3
init_block = 123456
init_balance_file = "./mock_balance.txt"
log_file = "./tx.log"
db_path = "./mock_db"
try:
    balance = defaultdict(int)
    addr_list = [] # used for sampling, ok to have duplicates (frequent addresses are more likely to be sampled)
    for i in range(init_size):
        addr = random.randint(0, max_addr)
        log_balance = random.uniform(0, log_max_init_balance)
        balance[addr] = int(2**log_balance)
        addr_list.append(addr)

    # create (if not exist) and clear log file
    open(log_file, "w").close()

    with open(init_balance_file, "w") as f:
        first_line = f"{init_block} {max_block_per_batch * max_transfer_per_block} {init_size}\n"
        f.write(first_line)
        for addr, bal in balance.items():
            # write addr as hex and balance as int
            f.write(f"{hex(addr)} {bal}\n")

    # run init_db
    os.system(f"rm -r {db_path}")
    os.system("g++ init_db_balance.cpp -o init_db_balance -L/usr/local/lib -lrocksdb")
    os.system(f"./init_db_balance {db_path} {init_balance_file}")

    # compile client program
    os.system("make client")

    # run mock_server_runner as a subprocess
    mock_server_runner = subprocess.Popen(["./mock_server_runner.sh"])
    time.sleep(20) # wait for the server to start

    transfers = []
    last_scanned_block = init_block + 10
    begin_time = time.time()
    while time.time() - begin_time < est_time:
        new_transfers = []
        num_blocks = random.randint(0, max_block_per_batch)
        for i in range(num_blocks + unstable_block_num):
            block_id = last_scanned_block - unstable_block_num + i + 1
            num_transfers = random.randint(0, max_transfer_per_block)
            for j in range(num_transfers):
                from_addr = random.choice(addr_list)
                if random.random() < prob_new_addr:
                    to_addr = random.randint(0, max_addr)
                else:
                    to_addr = random.choice(addr_list)
                amount = (int)(2**random.uniform(0, log_max_transfer_amount))
                tx_hash = random.getrandbits(256)
                tx_time = "2021-05-01T00:00:00Z" # mock time
                new_transfers.append((block_id, tx_hash, j, from_addr, to_addr, amount, tx_time))
                addr_list.append(to_addr)
                addr_list.append(from_addr)
        
        for tx in transfers:
            if tx[0] > last_scanned_block - unstable_block_num:
                from_addr = tx[3]
                to_addr = tx[4]
                amount = tx[5]
                # undo the previous unstable transfer
                # print(f"Undo the previous transfer: {hex(from_addr)} -> {hex(to_addr)} {amount}")
                balance[from_addr] += amount
                balance[to_addr] -= amount
        for tx in new_transfers:
            assert tx[0] > last_scanned_block - unstable_block_num
            from_addr = tx[3]
            to_addr = tx[4]
            amount = tx[5]
            # apply the transfer
            # print(f"Apply the transfer: {hex(from_addr)} -> {hex(to_addr)} {amount}")
            balance[from_addr] -= amount
            balance[to_addr] += amount
        transfers = new_transfers
        last_scanned_block += num_blocks
        # print("\n\nCurrent balance")
        # for addr, bal in balance.items():
        #     print(f"{hex(addr)} {bal}")
        # print("\n\n")
        with open(log_file, "at") as f:
            try:
                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                # write new transactions to the end of the file
                # print("Write new transactions to the log file")
                for tx in transfers:
                    line = f"{tx[0]}\t{hex(tx[1])}\t{tx[2]}\t{hex(tx[3])}\t{hex(tx[4])}\t{tx[5]}\t{tx[6]}"
                    # print(line)
                    f.write(line + "\n")
                # write a dummy transaction to update the last scanned block
                f.write(f"{last_scanned_block}\t0x0\t9999999\t0x0\t0x0\t0\t0\n")
            finally:
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)
        print(f"Last scanned block: {last_scanned_block}")
        time.sleep(10)
        print("Start testing clients")
        # run multiple clients using subprocess in parallel, make query through standard input, and print the result
        client_processes = []
        for i in range(client_count):
            client_processes.append(subprocess.Popen(["./client"], 
                                stdin=subprocess.PIPE, 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE, 
                                bufsize=1,  # Line buffered
                                universal_newlines=True))
        time.sleep(3)
        for i in range(client_count):
            while True:
                response = client_processes[i].stdout.readline()
                print(response)
                if response.strip() == "Enter an ethereum address:":
                    print("Test ethereum address ")
                    break
        for _ in range(10):
            query = hex(random.choice(addr_list))
            print("query ", query)
            for i in range(client_count):
                client_processes[i].stdin.write(query + "\n")
                client_processes[i].stdin.flush()
            for i in range(client_count):
                p = client_processes[i]
                response = p.stdout.readline()
                print("first response ", response)
                response = p.stdout.readline()
                print("second response ", response)
                balance_str = response.split()[1]
                balance_int = int(balance_str)
                refer_balance = balance[int(query, 16)]
                if refer_balance < 0:
                    refer_balance += 2**256
                assert balance_int == refer_balance
                response = p.stdout.readline()
                print("third response ", response)
            # terminate the client
        for p in client_processes:
            p.terminate()

        # sleep_time = random.uniform(min_interval_in_sec, max_interval_in_sec)
        # time.sleep(sleep_time)
finally:
    mock_server_runner.terminate()
