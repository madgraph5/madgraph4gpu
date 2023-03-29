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
secret = 'fV8dKViWTVdnA3Rw*qCeA@MYtZki@q'
Auth = ['db_user', secret]
physicsProcesses = ['ee_mumu', 'gg_ttggg', 'gg_ttgg', 'gg_ttg', 'gg_tt']
absLayers = ['SYCL', 'CUDA']
branch = 'master'
fields = ['EvtsPerSec[MatrixElems] (3)', 'EvtsPerSec[MECalcOnly] (3)']
reportPath = 'C:\\Users\\jteig\\cernbox\\Documents\\test\\22-12-07_cudacpp_Xeon-Silver-4216_v100s_gcc-11.3_cuda-11.6.2_master'

# Argument parser
parser = argparse.ArgumentParser(description='A script for sending data from profiler to InfluxDB.')

parser.add_argument('-r', '--reportPath', help="Path for the reports that is being put into the database.", default=reportPath)
parser.add_argument('-f', '--fields', help="Fields in the JSON to be put into the database.", default=fields)
#parser.add_argument('-g', '--gpu', help="GPU used when profiling.", default=GPU)
#parser.add_argument('--GCCVersion', help="GCC version used when profiling.", default=GCCVersion)
#parser.add_argument('--CUDAVersion', help="CUDA version used when profiling.", default=CUDAVersion)
parser.add_argument('-a', '--absLayer', help="Abstraction layer used when profiling.", default=absLayers[0])
parser.add_argument('-b', '--branch', help="Branch the profiler data is in.", default=branch)

# Fix this
parser.add_argument('-p', '--profiler', help="Enable CI profiling defaults.", default='0')

args = parser.parse_args()

#
#   Main
#
if __name__=='__main__':

    # Fix this
    if args.profiler == '1':

        if args.absLayer.upper() == "SYCL":

            syclNamePrefix = os.getenv('SYCL_NAME_PREFIX')
            
            if syclNamePrefix == None:
                logging.error('Sycl name prefix has not been set!')
                sys.exit(1)

            # Fix the branch detection from the file name here
            reportfolder= "workspace_mg4gpu/" + datetime.datetime.now().strftime('%y-%m-%d') + '_' + syclNamePrefix + '_' + args.branch
            print(reportfolder)

            if not os.path.exists(reportfolder):
                logging.error('SYCL report path does not exist!')
                sys.exit(1)

        elif args.absLayer.upper() == "CUDA":

            cudaNamePrefix = os.getenv('CUDA_NAME_PREFIX')
            if cudaNamePrefix == None:
                logging.error('Cuda name prefix has not been set!')
                sys.exit(1)

            reportfolder= "workspace_mg4gpu/" + datetime.datetime.now().strftime('%y-%m-%d') + '_' + cudaNamePrefix + '_' + args.branch

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

        with open(file, "r") as f:

            fileContents = f.read()

            if fileContents != '':
                data = json.loads(fileContents)

                fileName = (os.path.basename(file))

                for process in physicsProcesses:
                    if process in fileName.lower():
                        physicsProcess = process
                        break

                fileNameParts = fileName.split('_')

                CPU = fileNameParts[4]

                GPU = fileNameParts[5]

                GCCVersion = fileNameParts[6].split('-')[1]

                CUDAVersion = fileNameParts[7].split('-')[1]

                gridsize = data[0]["NumThreadsPerBlock"] * data[0]["NumBlocksPerGrid"]

                DBdata = f'{physicsProcess},CPU={CPU},GPU={GPU},AbstractionLayer={args.absLayer},GCCVersion={GCCVersion},CUDAVersion={CUDAVersion},NumThreadsPerBlock={data[0]["NumThreadsPerBlock"]},NumBlocksPerGrid={data[0]["NumBlocksPerGrid"]},NumIterations={data[0]["NumIterations"]} Gridsize={gridsize}'

                for field in fields:
                    value = float(re.findall(r'[\d.]+',data[0][field])[0])
                    
                    DBdata = DBdata + ',' + args.absLayer + "_" + field.replace(" ", "_") + '=' + str(value)

                requestInfo = ["curl", "-i", "-k",  '-XPOST', "-i",  URL, "--header",  "Authorization: Token "+Auth[0]+":"+Auth[1], "--data-raw", DBdata]
                
                request = subprocess.run(requestInfo, stdout=subprocess.DEVNULL)

                f.close()
                
                if request.returncode != 0:
                    print(str(datetime.datetime.now().strftime("%H:%M:%S")) + " Request FAILED! Data: " + DBdata)
                else:
                    print(str(datetime.datetime.now().strftime("%H:%M:%S")) + " Request COMPLETED! Data: " + DBdata)


            else: logging.error('No information/fields in the JSON report!')