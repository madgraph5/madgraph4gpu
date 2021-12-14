import os,sys,subprocess,time,math


threads_per_block_range = [2**x for x in range(3,9)]
total_threads_max = 22
iterations = 10

for threads_per_block in threads_per_block_range:
    blocks_max = int(math.log(2**total_threads_max / threads_per_block,2))
    blocks_range = [2**x for x in range(0,blocks_max)]
    print(f'{threads_per_block} == {blocks_range}')

    for blocks in blocks_range:
        cmd = f'{os.getcwd()}/ccheck.exe -j {blocks} {threads_per_block} {iterations}'
        print('running ',cmd)
        start = time.time()
        procout = subprocess.run(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        duration = time.time() - start
        print(f'completed: {procout.returncode} in {duration:10.2} seconds')


