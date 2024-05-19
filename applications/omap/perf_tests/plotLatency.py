import re
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogLocator
# Function to parse the text file and extract relevant information
def parse_text_file(file_path):
    data = {'map_size': [], 'thread_count': [], 'batch_size': [], 'find_time': [], 'init_time': [], 'insert_time': [], 'erase_time': []}

    with open(file_path, 'r') as file:
        lines = file.readlines()

    for line in lines:
        map_size_match = re.search(r'mapSize = (\d+)', line)
        thread_count_match = re.search(r'threadCount = (\d+)', line)
        batch_size_match = re.search(r'batchSize = (\d+)', line)
        find_time_match = re.search(r'oram find time (\d+\.\d+) us', line)
        insert_time_match = re.search(r'oram insert time (\d+\.\d+) us', line)
        erase_time_match = re.search(r'oram erase time (\d+\.\d+) us', line)
        init_time_match = re.search(r'oram init time (\d+\.\d+) s', line)

        if map_size_match:
            data['map_size'].append(int(map_size_match.group(1)))
        if thread_count_match:
            data['thread_count'].append(int(thread_count_match.group(1)))
        if batch_size_match:
            data['batch_size'].append(int(batch_size_match.group(1)))
        if find_time_match:
            data['find_time'].append(float(find_time_match.group(1)))
        if insert_time_match:
            data['insert_time'].append(float(insert_time_match.group(1)))
        if erase_time_match:
            data['erase_time'].append(float(erase_time_match.group(1)))
        if init_time_match:
            data['init_time'].append(float(init_time_match.group(1)))
    return data


def create_and_save_plot_init_single_thread(data, file_path, num_sub_omap=1):
    plt.figure()
    if num_sub_omap == 1:
        plt.title(f'Initialization Time\n(32-byte key, 32-byte value, 64 GB EPC, swap to SSD)')
    else:
        plt.title(f'Initialization Time ({num_sub_omap} sub-omaps)\n(32-byte key, 32-byte value, 64 GB EPC, swap to SSD)')
    plt.xlabel('Map Size')
    plt.ylabel('Initialization Time (s)')
    
    y_values = []
    sorted_map_sizes = sorted(set(data['map_size']), key=lambda x: x)
    x_values = sorted_map_sizes
    plotted_map_sizes = set()
    for map_size in sorted_map_sizes:
        for i in range(len(data['map_size'])):
            if data['map_size'][i] == map_size and map_size not in plotted_map_sizes:
                y_values.append(data['init_time'][i])  # Convert us to seconds
                plotted_map_sizes.add(map_size)
        # print(x_values, y_values)
    plt.plot(x_values, y_values)

    # plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.yscale('log')
    plt.xscale('log')
    # plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, subs=[1.0, 2.0, 5.0]))
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=[1.0, 2.0, 5.0]))
    # plt.gca().xaxis.set_major_formatter(ScalarFormatter())
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    # show all the x ticks
    plt.savefig(file_path + 'Init_Time.jpg', bbox_inches='tight')  # Save the plot as a JPG image

def create_and_save_plot_find_single_thread(data, file_path, batch_size=1):
    plt.figure()
    if batch_size == 1:
        plt.title(f'Single Thread Throughput\n(32-byte key, 32-byte value, 64 GB EPC, Swapped to SSD)')
    else:
        plt.title(f'Latency of different operations in batch of {batch_size}\n(32-byte key, 32-byte value, 64 GB EPC, swap to SSD)')
    plt.xlabel('Number of key-value pairs')
    if batch_size == 1:
        plt.ylabel('Latency (us)')
        rate = 1
    elif batch_size <= 1e5:
        plt.ylabel('Latency (ms)')
        rate = 1e3
    else:
        plt.ylabel('Latency (s)')
        rate = 1e6
        
    
    
    find_latency = []
    insert_latency = []
    erase_latency = []
    sorted_map_sizes = sorted(set(data['map_size']), key=lambda x: x)
    x_values = sorted_map_sizes
    for map_size in sorted_map_sizes:
        for i in range(len(data['map_size'])):
            if data['map_size'][i] == map_size and data['batch_size'][i] == batch_size:
                find_latency.append(data['find_time'][i] * batch_size / rate)
                insert_latency.append(data['insert_time'][i] * batch_size / rate)
                erase_latency.append(data['erase_time'][i] * batch_size / rate)
        
    # plt.plot(x_values, y_values)
    plt.plot(x_values, find_latency, label='Find')
    plt.plot(x_values, insert_latency, label='Insert')
    plt.plot(x_values, erase_latency, label='Erase')

    # plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    # if batch_size > 1:
    plt.yscale('log')
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=[1.0, 2.0, 5.0]))
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.xscale('log')
    # show legend
    plt.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0), borderaxespad=0.5)
   
    plt.savefig(file_path + f'Latency{batch_size}.jpg', bbox_inches='tight')  # Save the plot as a JPG image

# Example usage
# file_path = 'erc20SingleThread/'
file_path = 'pipeline/'
file_name = 'omap_perf_disk_64g.log'

parsed_data = parse_text_file(file_path + file_name)

# for map_size in set(parsed_data['map_size']):
#     create_and_save_plot_find(parsed_data, map_size, file_path)

# create_and_save_plot_init(parsed_data, file_path)

create_and_save_plot_find_single_thread(parsed_data, file_path)
create_and_save_plot_init_single_thread(parsed_data, file_path)

file_path = 'pardisk/'
file_name = 'par_30t_64g.log'
parsed_data = parse_text_file(file_path + file_name)
create_and_save_plot_find_single_thread(parsed_data, file_path, 1000)
create_and_save_plot_find_single_thread(parsed_data, file_path, 100000)
create_and_save_plot_init_single_thread(parsed_data, file_path, 15)
