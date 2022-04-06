import os,sys,subprocess,time,math,argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t','--total-threads-max',help='maximum total number of threads to test',default=22,type=int)
parser.add_argument('-i','--iterations',help='repeat measurement this many times in the program',default=10,type=int)
parser.add_argument('-e','--exe',help='executable name',default='ccheck.exe')
parser.add_argument('-a','--tpb-min',help='sets A in "2^A" which is the minimum number of threads per block (or team size) to test.',default=5,type=int)
parser.add_argument('-b','--tpb-max',help='sets B in "2^B" which is the maximum number of blocks (or league size) to test.',default=9,type=int)
args = parser.parse_args()

threads_per_block_range = [2**x for x in range(args.tpb_min,args.tpb_max)]
total_threads_max = args.total_threads_max
iterations = args.iterations

exe = args.exe

for threads_per_block in threads_per_block_range:
    blocks_max = int(math.log(2**total_threads_max / threads_per_block,2))
    blocks_range = [2**x for x in range(0,blocks_max)]
    print(f'{threads_per_block} == {blocks_range}')

    for blocks in blocks_range:
        cmd = f'{os.getcwd()}/{exe} -j {blocks} {threads_per_block} {iterations}'
        print('running ',cmd)
        start = time.time()
        procout = subprocess.run(cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        duration = time.time() - start
        print(f'completed: {procout.returncode} in {duration:10.2} seconds')


