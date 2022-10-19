import sys
import subprocess
import datetime

# Required info
absLayer = "CUDA"
doublePrecisionConstant = 2560
mgProcesses = ["ee_mumu", "gg_tt", "gg_ttg", "gg_ttgg", "gg_ttggg"]
iterations = 10
threadsPerBlock = [32, 64, 128, 256]
blocksPerGrid = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]

# How many runs in total the program made
count = 0

for process in mgProcesses:
    for TPB in threadsPerBlock:
        for BPG in blocksPerGrid:
            if (TPB * BPG > doublePrecisionConstant):

                if absLayer.upper() == 'SYCL':
                    args = ["./buildMadGraphFile.sh", "-n",  process, "-i",  str(iterations), "-t",  str(TPB), "-b", str(BPG)]

                elif absLayer.upper() == 'CUDA':

                    ### Used in br_golden_epochX4 branch
                    #if ".sa" not in process:
                    #    process = process + ".sa"
                    ###
                    
                    args = ["./buildCUDAFile.sh", "-n",  process, "-i",  str(iterations), "-t",  str(TPB), "-b", str(BPG)]

                else: sys.exit("No abstraction layer selected")

                print(str(datetime.datetime.now().strftime("%H:%M:%S")) + " Started " + process + " with TPB("+ str(TPB) +") * BPG("+ str(BPG) +"): " + str(TPB * BPG) + "!")
                
                build = subprocess.run(args, stdout=subprocess.DEVNULL)
                if build.returncode != 0:
                    print(str(datetime.datetime.now().strftime("%H:%M:%S")) + " " + process + " FAILED!, threadsPerBlock: " + str(TPB) + ", blocksPerGrid: " + str(BPG) + ", Product: " + str(TPB * BPG))
                else:
                    print(str(datetime.datetime.now().strftime("%H:%M:%S")) + " " + process + " COMPLETED!, threadsPerBlock: " + str(TPB) + ", blocksPerGrid: " + str(BPG) + ", Product: " + str(TPB * BPG))
                count += 1

print("Builded " + str(count) + " processes!")