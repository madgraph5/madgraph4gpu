import platform
import psutil
import subprocess
import re
import multiprocessing
import os

cpu_info = {}


def get_cpu_name():
    try:
        output = subprocess.check_output(['lscpu'])
        lines = output.decode('utf-8').strip().split('\n')
        for line in lines:
            if 'Model name:' in line:
                _, cpu_name = line.split(':', 1)
                return cpu_name.strip()
    except subprocess.CalledProcessError:
        pass

    return None


cpu_name = get_cpu_name()
if cpu_name:
    cpu_info["CPU Name"] = cpu_name
else:
    cpu_info["CPU Name"] = ("Unable to determine CPU name.")


def system_info():
    cpu_info["Platform processor"] = platform.processor()
    cpu_info["Platform architecture"] = platform.architecture()
    cpu_info["Machine type"] = platform.machine()
    cpu_info["Machine platform information"] = platform.platform()
    cpu_info["Operating System"] = platform.system()


def get_cpu_info():
    cpu_info["CPUcorecount"] = multiprocessing.cpu_count()
    cpu_info["CPUusage"] = psutil.cpu_percent(interval=0.01)


def print_memory_info():
    mem = psutil.virtual_memory()
    total_memory = mem.total / (1024 ** 3)
    available_memory = mem.available / (1024 ** 3)
    used_memory = total_memory - available_memory
    cpu_info["memorytot"] = total_memory
    cpu_info["memoryused"] = used_memory
    cpu_info["memoryavailable"] = available_memory
    cpu_info["memorypercent"] = mem.percent


def print_disk_info():
    disk = psutil.disk_usage('/')
    total_space = disk.total / (1024 ** 3)
    used_space = disk.used / (1024 ** 3)
    free_space = disk.free / (1024 ** 3)
    cpu_info["totspace"] = total_space
    cpu_info["usedspace"] = used_space
    cpu_info["freespace"] = free_space
    cpu_info["diskutil"] = disk.percent


def print_network_info():
    net_io = psutil.net_io_counters()
    cpu_info["bytessent"] = net_io.bytes_sent
    cpu_info["bytesreceived"] = net_io.bytes_recv


def extract_vectoring_info():
    command = 'lscpu'
    output = subprocess.check_output(command, shell=True, universal_newlines=True)

    cpuinfo = ''
    lines = output.splitlines()
    for line in lines:
        if line.startswith('CPU op-mode(s):'):
            cpuinfo += line + '\n'
        elif line.startswith('Flags:'):
            cpuinfo += line + '\n'

    with open('cpuinfo.txt', 'w') as file:
        file.write(cpuinfo)

    with open('cpuinfo.txt', 'r') as file:
        cpuinfo = file.read()

    os.remove('cpuinfo.txt')

    avx_matches = re.findall(r'avx([\d_]+)', cpuinfo)
    avx_numbers = [float(num.replace('_', '.')) for num in avx_matches]

    if avx_numbers:
        max_avx_number = max(avx_numbers)
        max_avx_string = 'avx' + str(int(max_avx_number))
        cpu_info["maxavxstring"] = max_avx_string
    else:
        sse_matches = re.findall(r'sse([\d_]+)', cpuinfo)
        sse_numbers = [float(num.replace('_', '.')) for num in sse_matches]
        if sse_numbers:
            max_sse_number = max(sse_numbers)
            max_sse_string = 'sse' + str(max_sse_number).replace('.', '_')
            cpu_info["maxssestring"] = max_sse_string
        else:
            cpu_info["sseavxnotfound"] = "AVX and SSE have not been found"


def print_cpu_info():
    global cpu_info  # Declare cpu_info as a global variable
    get_cpu_name()
    system_info()
    get_cpu_info()
    print_memory_info()
    print_disk_info()
    print_network_info()
    extract_vectoring_info()
    print(f' The plaform processor is {cpu_info["CPU Name"]} \n Platform architecture is {cpu_info["Platform architecture"]} \n')
    print(f' The Machine type is {cpu_info["Machine type"]} \n Machine platform information: {cpu_info["Machine platform information"]} \n ')
    print(f' The operating system is {cpu_info["Operating System"]} \n CPU core count: {cpu_info["CPUcorecount"]} \n CPU usage:  {cpu_info["CPUusage"]} %\n')
    print(f' total memory: {cpu_info["memorytot"]} GB \n  Memory used: {cpu_info["memoryused"]} GB \n Memory available: {cpu_info["memoryavailable"]} GB \n')
    print(f' Memory usage {cpu_info["CPUusage"]} % \n Total disk space: {cpu_info["totspace"]} GB \n Used disk space: {cpu_info["usedspace"]} GB \n ')
    print(f' Free disk space:  {cpu_info["freespace"]} GB \n Disk utilisation: {cpu_info["diskutil"]} GB  \n Network Information: ')
    print(f' Bytes Sent: {cpu_info["bytessent"]} \n Bytes received: {cpu_info["bytesreceived"]}')
    print(f' Vectoring level: {cpu_info.get("maxavxstring", "")} {cpu_info.get("maxssestring", "")} 
{cpu_info.get("sseavxnotfound", "")}')


print_cpu_info()


