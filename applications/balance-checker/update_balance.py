# Read a txt file of erc20 balance till a certain block number and transaction id
# Then read multiple log files of erc20 transfer events and update the balance
# Finally, write the updated balance to a new txt file

# The balance file should be in the following format:
# 18784420 230 5574
# 0xb0255cc0ad302545bd1dad0b39af5b1492f7a4b3 4348029885273000000000000
# 0xab697d998be7d4a215b241a31b7fa8615cccea6a 1533670505197973256865817
# ...
# The first line is the block number, transaction id and the number of addresses
# The following lines are the addresses and the balance

# Each line of the log files should be in the following format:
# 19094970        0x7822aa2b164b4b7dbb0a49dc8bde35b7275c30a23f59ad335a027dbb02a42e13      79      0x974caa59e49682cda0ad2bbe82983419a2ecc400      0x08ccadb7724eef6fa927796a120b87cd0c82addb       1159000000      2024-01-27T02:34:47
# The first column is the block number, the second column is the transaction id, the third column is the log index, the fourth column is the from address, the fifth column is the to address, the sixth column is the amount, the seventh column is the timestamp
# The log files are sorted by block number and transaction id
# There could be duplicate lines in the log files, so if the block number is smaller than a previous line, the line should be ignored (if the block number is the same, compare the transaction id)

# The updated balance file should be in the same format as the original balance file
# The log files could be very large, so the program should be efficient (don't store the logs in memory, update the balance in a streaming fashion)

def read_initial_balances(balance_file):
    with open(balance_file, 'r') as f:
        lines = f.readlines()
        
    metadata = lines[0].strip().split()
    balances = {}
    
    for line in lines[1:]:
        address, balance = line.strip().split()
        balances[address] = int(balance)
        
    return int(metadata[0]), int(metadata[1]), int(metadata[2]), balances

def update_balances_from_logs(balance_file, log_files, updated_balance_file):
    # Read initial balances
    last_block_number, last_log_index, num_addresses, balances = read_initial_balances(balance_file)
    
    for log_file in log_files:
        with open(log_file, 'r') as f:
            for line in f:
                log_block_number, transaction_id, log_index, from_address, to_address, amount, timestamp = line.strip().split()
                log_block_number = int(log_block_number)
                log_index = int(log_index)
                # transaction_id = int(transaction_id)
                amount = int(amount)
                
                # Check if the log entry is relevant (i.e., after the last processed block/transaction)
                if log_block_number < last_block_number or (log_block_number == last_block_number and log_index <= last_log_index):
                    continue
                
                # Update balances
                if from_address in balances:
                    balances[from_address] -= amount
                else:
                    balances[from_address] = -amount
                
                if to_address in balances:
                    balances[to_address] += amount
                else:
                    balances[to_address] = amount
                
                # Update last processed block and transaction ID
                last_block_number = log_block_number
                last_log_index = log_index

    # remove addresses with zero balance
    balances = {address: balance for address, balance in balances.items() if balance > 0}
    
    # Write updated balances to a new file, sorted by balance (descending)
    with open(updated_balance_file, 'w') as f:
        f.write(f"{last_block_number} {last_log_index} {len(balances)}\n")
        for address, balance in sorted(balances.items(), key=lambda x: x[1], reverse=True):
            f.write(f"{address} {balance}\n")

# Example usage:
balance_file = 'usdt_balance.txt'
log_files = ['usdt_tx.log', 'usdt_tx2.log']  # List of log file paths
updated_balance_file = 'updated_usdt_balances.txt'

update_balances_from_logs(balance_file, log_files, updated_balance_file)
