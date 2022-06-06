import subprocess
import numpy as np
import os

base = os.getcwd()


paths = [

]

p = subprocess.Popen('mpirun -n 56 ./check',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
stdout,stderr = p.communicate()
stdout = stdout.decode("utf-8") 

value = []
error = []
for line in stdout.split('\n'):
   if 'time' in line:
      print(line)
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

print("%10.6e %10.6e" % (mean,std))

