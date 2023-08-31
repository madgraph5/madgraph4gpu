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
#   Python script for sending generated reports from performance profiling to InfluxDB instance 
#   using the MadGraph5_aMC@NLO GPU development framework
#
#   Author: Jorgen Teig, CERN 2023
#

import os
import glob
import json
import re
import logging
import subprocess
import datetime
import argparse
import sys

# Parameter defaults
URL = 'https://dbod-madgraph4gpu-db.cern.ch:8082/api/v2/write?bucket=ProfilerData'
secret = os.environ.get('MADGRAPH4GPU_DB_SECRET')
AUTH = ['db_user', secret]
PHYS_PROCESSES = ['ee_mumu', 'gg_ttggg', 'gg_ttgg', 'gg_ttg', 'gg_tt']
ABS_LAYERS = ['SYCL', 'CUDA', 'HIP']
BRANCH = 'master'
FIELDS = ['EvtsPerSec[MatrixElems] (3)', 'EvtsPerSec[MECalcOnly] (3)']

# Default reportPath (Useful for testing)
REPORT_PATH = 'C:\\Users\\jteig\\cernbox\\Documents\\test\\22-12-07_cudacpp_Xeon-Silver-4216_v100s_gcc-11.3_cuda-11.6.2_master'

# Argument parser
parser = argparse.ArgumentParser(description='A script for sending data from profiler to InfluxDB.')

parser.add_argument('-r', '--reportPath', help="Path for the reports that is being put into the database.", default=REPORT_PATH)
parser.add_argument('-f', '--fields', help="Fields in the JSON to be put into the database.", default=FIELDS)
parser.add_argument('-a', '--absLayer', help="Abstraction layer used when profiling.", default=ABS_LAYERS[0])
parser.add_argument('-b', '--branch', help="Branch the profiler data is in.", default=BRANCH)
parser.add_argument('-p', '--profiler', help="Enable CI profiling defaults.", default='0')

args = parser.parse_args()

#
#   Main
#
if __name__=='__main__':

    # Sets report path for extracting the reports generated from performanceProfiler.py 
    if args.profiler == '1':

        if args.absLayer.upper() == "SYCL":

            syclNamePrefix = os.getenv('SYCL_NAME_PREFIX')
            
            if syclNamePrefix is None:
                logging.error('Sycl name prefix has not been set!')
                sys.exit(1)

            reportfolder= "workspace_mg4gpu/" + datetime.datetime.now().strftime('%y-%m-%d') + '_' + syclNamePrefix + '_' + args.branch

            if not os.path.exists(reportfolder):
                logging.error('SYCL report path does not exist!')
                sys.exit(1)

        elif args.absLayer.upper() == "CUDA":

            cudaNamePrefix = os.getenv('CUDA_NAME_PREFIX')

            if cudaNamePrefix is None:
                logging.error('Cuda name prefix has not been set!')
                sys.exit(1)

            reportfolder= "workspace_mg4gpu/" + datetime.datetime.now().strftime('%y-%m-%d') + '_' + cudaNamePrefix + '_' + args.branch

            if not os.path.exists(reportfolder):
                logging.error('CUDA report path does not exist!')
                sys.exit(1)

        elif args.absLayer.upper() == "HIP":

            hipNamePrefix = os.getenv('HIP_NAME_PREFIX')

            if cudaNamePrefix is None:
                logging.error('HIP name prefix has not been set!')
                sys.exit(1)

            reportfolder= "workspace_mg4gpu/" + datetime.datetime.now().strftime('%y-%m-%d') + '_' + hipNamePrefix + '_' + args.branch

            if not os.path.exists(reportfolder):
                logging.error('CUDA report path does not exist!')
                sys.exit(1)

        else:
            logging.error('No abstraction layer that is supported has been selected!')
            sys.exit(1)

    else:
        reportfolder = args.reportPath

    filePath = []
    filePath.append(glob.glob(reportfolder + '/test_*.json'))
    filePath.append(glob.glob(reportfolder + '/*/test_*.json'))

    # Flatten the list
    files = [p for sublist in filePath for p in sublist]

    for file in files:

        with open(file, "r", encoding='utf-8') as f:

            fileContents = f.read()

            if fileContents != '':
                data = json.loads(fileContents)

                fileName = (os.path.basename(file))

                for process in PHYS_PROCESSES:
                    if process in fileName.lower():
                        physicsProcess = process
                        break

                fileNameParts = fileName.split('_')

                CPU = fileNameParts[4]

                GPU = fileNameParts[5]

                GCCVersion = fileNameParts[6].split('-')[1]

                GPUVersion = fileNameParts[7].split('-')[1]

                gridsize = data[0]["NumThreadsPerBlock"] * data[0]["NumBlocksPerGrid"]

                DBdata = f'{physicsProcess},CPU={CPU},GPU={GPU},AbstractionLayer={args.absLayer},GCCVersion={GCCVersion},GPUVersion={GPUVersion},NumThreadsPerBlock={data[0]["NumThreadsPerBlock"]},NumBlocksPerGrid={data[0]["NumBlocksPerGrid"]},NumIterations={data[0]["NumIterations"]} Gridsize={gridsize}'

                for field in FIELDS:
                    value = float(re.findall(r'[\d.]+',data[0][field])[0])
                   
                    DBdata = DBdata + ',' + args.absLayer + "_" + field.replace(" ", "_") + '=' + str(value)

                requestInfo = ["curl", "-i", "-k",  '-XPOST', "-i",  URL, "--header",  "Authorization: Token "+AUTH[0]+":"+AUTH[1], "--data-raw", DBdata]
               
                request = subprocess.run(requestInfo, stdout=subprocess.DEVNULL, check=True)

                f.close()
               
                if request.returncode != 0:
                    print(str(datetime.datetime.now().strftime("%H:%M:%S")) + " Request FAILED! Data: " + DBdata)
                else:
                    print(str(datetime.datetime.now().strftime("%H:%M:%S")) + " Request COMPLETED! Data: " + DBdata)


            else: logging.error('No information/fields in the JSON report!')
