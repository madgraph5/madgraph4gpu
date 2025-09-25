#Block for General Code
#This code will only contain the parts of Intel that actually work
import re
import subprocess
gpu_info = {}

#_____________________________________________________
#GPU VENDOR
#Block for finding GPU Vendor
def get_gpu_vendor():
    lspci_output = subprocess.check_output(['lspci', '-nnk'], universal_newlines=True)
    if 'NVIDIA' in lspci_output:
        gpu_info["vendor"] = "NVIDIA"
        return "nvidia"
    elif 'AMD' in lspci_output:
        gpu_info["vendor"] = "AMD"
        return "amd"
    elif 'Intel' in lspci_output:
        gpu_info["vendor"] = "Intel"
        return "intel" 
    else:
       gpu_info["None"]
#________________________________________________
# CARD TYPE
def getgpucardtype(vendor):
    if vendor == "intel":
       gpu_info["card_type"] = "dysfunctional (Intel)"
       return None
    try:

        # Execute the 'lspci' command to list PCI devices
        lspci_output = subprocess.check_output(['lspci'], universal_newlines=True)

        # Search for NVIDIA GPUs and VGA compatible controllers
        nvidia_lines = [line for line in lspci_output.split('\n') if 'NVIDIA' in line]
        vga_lines = [line for line in lspci_output.split('\n') if 'VGA compatible controller' in line]

        # Extract the GPU card types from the NVIDIA lines
        for line in nvidia_lines:
            card_type = line.split('NVIDIA')[1]
            gpu_info["card_type"] = card_type

        # Extract the GPU card types from the VGA lines
        for line in vga_lines:
            card_type = line.split('VGA compatible controller: ')[1]
            gpu_info["card_type"] = card_type

        # Execute the 'lspci' command with the appropriate arguments to list AMD GPUs
        amd_output = subprocess.check_output(['lspci', '-nnk', '-d', '1002:'], universal_newlines=True)

        # Search for AMD GPUs in the output
        amd_lines = [line for line in amd_output.split('\n') if 'VGA compatible controller' in line]

        # Extract the GPU card types from the AMD lines
        for line in amd_lines:
            card_type = line.split('VGA compatible controller: ')[1]
            gpu_info["card_type"] = card_type 
    except subprocess.CalledProcessError:
        pass

    return None

#________________________________________________


def get_gpu_count(vendor):
    gpu_count = 0
    if vendor == "intel":
       gpu_info['gpucount'] = "dysfunctional (intel code)"
       return None
    try:
        # Check for NVIDIA GPUs using the 'nvidia-smi' command
        nvidia_smi_output = subprocess.check_output(['nvidia-smi', '--list-gpus'], universal_newlines=True)
        gpu_count += len(nvidia_smi_output.strip().split('\n'))
				
    except subprocess.CalledProcessError:
        pass

    try:
        # Check for AMD GPUs using the 'lspci' command
        lspci_output = subprocess.check_output(['lspci', '-nnk'], universal_newlines=True)
        amd_lines = [line for line in lspci_output.split('\n') if 'AMD' in line and 'VGA' in line]
        gpu_count += len(amd_lines)
				
    except subprocess.CalledProcessError:
        pass
    gpu_info['gpucount'] = gpu_count  # Add the GPU count to the dictionary

 

   


 #_____________________________________________________

#Block for GPU Memory
def get_gpu_memory(vendor):
    if vendor == "amd":
        command = "lspci -nnk | grep -A 2 'VGA' | grep 'Memory'"
        output = subprocess.check_output(command, shell=True, universal_newlines=True)
        memory_info = output.split(":")[1].strip()
        gpu_info["memory"] = memory_info
    elif vendor == "intel":
        gpu_info["memory"] = "Intel code dysfunctional"
        #No answer
        return None
    elif vendor == "nvidia":
        command = ["nvidia-smi", "--query-gpu=memory.total", "--format=csv,noheader"]
        process = subprocess.run(command, capture_output=True, text=True)
        output = process.stdout.strip()
        gpu_info["memory"] = output
     
    else:
        gpu_info["memory"] = "Unknown"

#____________________________________________________

#Block for Number of GPU streaming multiprocessors
def get_num_sm(vendor):
    if vendor == "amd":
        command = "lspci -nnk | grep -A 2 'VGA' | grep 'device'"
        output = subprocess.check_output(command, shell=True, universal_newlines=True)
        num_sm_info = output.split(":")[1].strip()
        gpu_info["num_sm"] = num_sm_info
    elif vendor == "intel":
        gpu_info["num_sm"] = "dysfunctional (Intel) "
        gpu_info["cuda"] = "N/A"
        return None
    elif vendor == "nvidia":
        command = ["/usr/local/cuda-12.0/extras/demo_suite/deviceQuery"]
        process = subprocess.run(command, capture_output=True, text=True)
        output = process.stdout.strip()
    # Search for the line containing the number of processors and CUDA cores
        pattern = r"\((\d+)\) Multiprocessors, \(\s*(\d+)\) CUDA Cores/MP"
        match = re.search(pattern, output)
        if match:
            num_processors = int(match.group(1))
            cuda_cores_per_processor = int(match.group(2))
            gpu_info["num_sm"] = num_processors
            gpu_info["cuda"] = cuda_cores_per_processor
    else:
        gpu_info["num_sm"] = "Unknown"
    #Will have to see about this


#_________________________________________________

def get_gpu_info(vendor):
    if vendor == "intel":
      gpu_info["tensorcores"] = "Intel Dysfunctional"
      return None
    tensor_count = 0
    try:
        # Check for GPU information using the 'lspci' command
        lspci_output = subprocess.check_output(['lspci', '-nnk'], universal_newlines=True)
        lines = lspci_output.split('\n')
        if vendor.lower() == 'amd' or vendor.lower() == 'intel':
           tensor_count = "None (AMD/Intel) device"
           return None
        elif vendor.lower() == 'nvidia':
            # Find NVIDIA GPU lines in the lspci output
            nvidia_lines = [line.strip() for line in lines if 'NVIDIA' in line]
            for line in nvidia_lines:
                if 'a100' in line.lower():
                    tensor_count = 108
                elif 'v100' in line.lower():
                    tensor_count = 640
                elif 'h100' in line.lower():
                     pass
                     #In the future
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    gpu_info["tensorcores"] = tensor_count
#________________________________________________


#This Block will print out all the necessary info
def printcode():
   vendor = get_gpu_vendor()
   get_gpu_vendor()
   getgpucardtype(vendor)
   get_gpu_count(vendor)
   get_gpu_memory(vendor) 
   get_num_sm(vendor)
   get_gpu_info(vendor)
   print(f' The GPU Vendor is {gpu_info["vendor"]} \n The GPU card type is {gpu_info["card_type"]} \n The number of GPUS are {gpu_info["gpucount"]} \n The number of cuda cores per Multi-processor is {gpu_info["cuda"]}\n The memory information is {gpu_info["memory"]} \n The number of streaming multiprocessors {gpu_info["num_sm"]} \n The number of tensor cores are {gpu_info["tensorcores"]}  ')
printcode()
