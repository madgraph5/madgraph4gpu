import subprocess
import numpy as np
import os

base = os.getcwd()


paths = [
'ee_mumu/SubProcesses/P1_epem_mupmum',
'gg_tt/SubProcesses/P1_gg_ttx',
'gg_ttg/SubProcesses/P1_gg_ttxg',
'gg_ttgg/SubProcesses/P1_gg_ttxgg',
'gg_ttggg/SubProcesses/P1_gg_ttxggg',
]

for path in paths:
   p = subprocess.Popen(f'cd {path}; mpirun -n 56 ./check',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   stdout,stderr = p.communicate()
   stdout = stdout.decode("utf-8") 

   value = []
   error = []
   for line in stdout.split('\n'):
      if 'time' in line:
         #print(line)
         try:
            values = line.split('=')[1]
            v,e = values.split('+-')
            value.append(float(v))
            error.append(float(e))
         except: continue

   value = np.array(value)
   error = np.array(error)

   mean = np.mean(value)
   std = np.sqrt(np.sum(np.square(error))) / error.shape[0]

   print("%50s: %10.6e +/- %10.6e  rate:  %10.6e" % (path,mean,std,1/mean*error.shape[0]))

