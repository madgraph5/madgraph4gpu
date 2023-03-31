#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
#   __  __               _    ____                          _       _  _      ____   ____    _   _ 
#  |  \/  |   __ _    __| |  / ___|  _ __    __ _   _ __   | |__   | || |    / ___| |  _ \  | | | |
#  | |\/| |  / _` |  / _` | | |  _  | '__|  / _` | | '_ \  | '_ \  | || |_  | |  _  | |_) | | | | |
#  | |  | | | (_| | | (_| | | |_| | | |    | (_| | | |_) | | | | | |__   _| | |_| | |  __/  | |_| |
#  |_|  |_|  \__,_|  \__,_|  \____| |_|     \__,_| | .__/  |_| |_|    |_|    \____| |_|      \___/ 
#                                                  |_|                                             
#
#
#   Python script for performance profiling using the MadGraph5_aMC@NLO GPU development framework
#
#   Author: Jorgen Teig, CERN 2023
#

import sys
import subprocess
import datetime
import argparse

# Parser arguments defaults
ABS_LAYER = "SYCL"
BRANCH = "master"

# Physics processes
MG_PROCESSES_SA = ["ee_mumu.sa", "gg_tt.sa", "gg_ttg.sa", "gg_ttgg.sa", "gg_ttggg.sa"]

DOUBLE_PRECISION_CONSTANT = 2560
ITERATIONS = 10
THREADS_PER_BLOCK = [256]
#THREADS_PER_BLOCK = [32, 64, 128, 256]
BLOCKS_PER_GRID = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]

# Parser
parser = argparse.ArgumentParser(description='A program for profiling GPUs using MadGraph.')

parser.add_argument("-l", help="Choose which abstraction layer you want to use (CUDA/SYCL).", default=ABS_LAYER)
parser.add_argument("-b", help="Choose which branch the madgraph4gpu repo is in.", default=BRANCH)

pyArgs = parser.parse_args()

# How many runs in total the program made
count = 0

for process in MG_PROCESSES_SA:
    for TPB in THREADS_PER_BLOCK:
        for BPG in BLOCKS_PER_GRID:
            if TPB * BPG > DOUBLE_PRECISION_CONSTANT:

                if pyArgs.l.upper() == 'SYCL':

                    # There is no .sa in br_golden_epochX4
                    # so it makes sure that .sa is included in everything other than that branch
                    # if pyArgs.b != 'br_golden_epochX4':
                    #if ".sa" not in process:
                    #    process = process + ".sa"

                    bashArgs = ["./buildSYCLProcess.sh",
                                "-n", process,
                                "-i", str(ITERATIONS),
                                "-t", str(TPB),
                                "-b", str(BPG),
                                "-r", str(pyArgs.b).lower()]

                elif pyArgs.l.upper() == 'CUDA':

                    bashArgs = ["./buildCUDAProcess.sh",
                                "-n", process,
                                "-i", str(ITERATIONS),
                                "-t", str(TPB),
                                "-b", str(BPG),
                                "-r", str(pyArgs.b).lower()]

                else: sys.exit("No abstraction layer matching the supplied string!")

                time = str(datetime.datetime.now().strftime("%H:%M:%S"))

                print(time + " Started " + process + " with TPB("+ str(TPB) +") * BPG("+ str(BPG) +"): " + str(TPB * BPG) + "!")

                build = subprocess.run(bashArgs, check=True)#, stdout=subprocess.DEVNULL)
                if build.returncode != 0:
                    print(time + " " + process +
                          " FAILED!, threadsPerBlock: " + str(TPB) +
                                    ", blocksPerGrid: " + str(BPG) +
                                    ", Product: " + str(TPB * BPG))
                else:
                    print(time + " " + process +
                          " COMPLETED!, threadsPerBlock: " + str(TPB) +
                                    ", blocksPerGrid: " + str(BPG) +
                                    ", Product: " + str(TPB * BPG))

                count += 1

print("Builded " + str(count) + " processes!")