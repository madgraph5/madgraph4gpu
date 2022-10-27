import sys
import os
import subprocess
import datetime
import argparse

# Parser arguments defaults
absLayer = "SYCL"
branch = "br_golden_epochX4"

mgProcesses = ["ee_mumu", "gg_tt", "gg_ttg", "gg_ttgg", "gg_ttggg"]

doublePrecisionConstant = 2560
iterations = 10
threadsPerBlock = [32, 64, 128, 256]
blocksPerGrid = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]

# Parser
parser = argparse.ArgumentParser(description='A program for profiling GPUs using MadGraph.')

parser.add_argument("-l", "--layer", help="Choose which abstraction layer you want to use (CUDA/SYCL).", default=absLayer)
parser.add_argument("-b", "--branch", help="Choose which branch the madgraph4gpu repo is in.", default=branch)

args = parser.parse_args()

# How many runs in total the program made
count = 0

for process in mgProcesses:
    for TPB in threadsPerBlock:
        for BPG in blocksPerGrid:
            if (TPB * BPG > doublePrecisionConstant):

                if args.l.upper() == 'SYCL':
                    args = ["./buildSYCLProcess.sh", "-n",  process, "-i",  str(iterations), "-t",  str(TPB), "-b", str(BPG)]

                elif args.l.upper() == 'CUDA':

                    # Used in br_golden_epochX4 branch
                    if args.b == 'br_golden_epochX4':
                        if ".sa" not in process:
                            process = process + ".sa"
                    
                    args = ["./buildCUDAProcess.sh", "-n",  process, "-i",  str(iterations), "-t",  str(TPB), "-b", str(BPG)]

                else: sys.exit("No abstraction layer matching the supplied string!")

                print(str(datetime.datetime.now().strftime("%H:%M:%S")) + " Started " + process + " with TPB("+ str(TPB) +") * BPG("+ str(BPG) +"): " + str(TPB * BPG) + "!")
                
                build = subprocess.run(args, stdout=subprocess.DEVNULL)
                if build.returncode != 0:
                    print(str(datetime.datetime.now().strftime("%H:%M:%S")) + " " + process + " FAILED!, threadsPerBlock: " + str(TPB) + ", blocksPerGrid: " + str(BPG) + ", Product: " + str(TPB * BPG))
                else:
                    print(str(datetime.datetime.now().strftime("%H:%M:%S")) + " " + process + " COMPLETED!, threadsPerBlock: " + str(TPB) + ", blocksPerGrid: " + str(BPG) + ", Product: " + str(TPB * BPG))
                count += 1

print("Builded " + str(count) + " processes!")