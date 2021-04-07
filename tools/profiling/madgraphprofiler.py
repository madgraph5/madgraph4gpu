#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 10:25:46 2021

@author: andy
"""
import platform
import subprocess
import os
import configparser
import sys
import shutil
import time

class PerformanceTesting:
    root_path=None
    configvalues = {}
    cpu = None
    gpus = None
    
    
    def rootfolfdercreation(self):
        print('Creating root folder\n')
        PerT.root_path = '/tmp/profile/root'
        if os.path.exists(PerT.root_path):
            print('root directory %s already existing.\n' % PerT.root_path)
            self.delete=input('Do you want to refresh the folder structure? y/n\n')
            if self.delete=='y':
                shutil.rmtree(PerT.root_path)
                os.makedirs(PerT.root_path) #Create root directory
            elif self.delete=='n':
                pass
            else:
                print('invalid input. Keep current folder structure')
                
        else:
        #Create root directory
            try:
                os.makedirs(PerT.root_path) #Create root directory
            except OSError:
                print('Creation of the directory %s failed' % PerT.root_path)
                print('Root directory can neither be found nor created.\nExecution failed. ')
                sys.exit()
            else:
                print('Successfully created the directory %s\n' % PerT.root_path)        
   
    def readConfig(self):
        print('Including values from configuration file\n')
        self.config = configparser.ConfigParser()
        self.file='profileconfig.ini'
        self.config.read(self.file)
        self.sections=self.config.sections()
        for i in range(len(self.sections)): 
            PerT.configvalues[self.sections[i]] = dict(self.config.items(self.sections[i]))
        print('Configuration completed\n')
        
    def CPU(self):
        print('Investigating CPU\n')
        cpu=platform.processor()
        print(f'%s is used as CPU device\n' %cpu)
    
    def GPU(self):
        print('Investigating GPU\n')
        self.gpu_hw_list=['NVIDIA','AMD','Asus','Intel','EVGA','Inno','Gigabyte','Zotac']
        
        
        self.getHW = subprocess.run(['lshw', '-C', 'display'], capture_output=True)
        self.display_list=str(self.getHW.stdout).split('\\n')
        #check which GPUs are installed and print it out
        self.found = False            
        for self.gpu in self.gpu_hw_list:
            
            vendor = [desc_item for desc_item in self.display_list if self.gpu in desc_item]
            if(len(vendor) == 0):
                print(f'Not found: {self.gpu}')
            else:
                print(f'Found:     {self.gpu}')
                PerT.gpus = self.gpu
        
    def cuda(self):
        print('Installing Cuda version 11.1')
        self.curdir=os.getcwd()
        print(self.curdir)
        
        self.curdir=str(os.getcwd())+'/GPU.sh'
        print(self.curdir)
        try:
            self.cuda = subprocess.run(['sh',self.curdir],capture_output=True)
            print('Call %s' %self.curdir+ '/GPU.sh')   
            #Check installed Cuda Version
        except OSError:
            print('\nMission failed\nGPU.sh not executed')
 
        k = None
        while k != 1:
            
            if os.getcwd() == PerT.root_path:
                self.clone=subprocess.run(['git','clone','https://github.com/madgraph5/madgraph4gpu.git'],capture_output=True)
                print(self.clone)
                k=1
            else:
                os.chdir(PerT.root_path)
                print('Changed directory to '+PerT.root_path)
                
    def check_gcheck(self):
        #Create folder structure according to configuration file and execute check or gcheck
        self.cg=PerT.configvalues['runfiles']['sys'] #Get check or gcheck
        print('Assemble path according to configuration')
        self.ex_path=PerT.assemblingpath() # Assemble path 
        print('Hop in '+self.ex_path+'\n') 
        os.chdir(self.ex_path)
        print('Execute make file.\nIt may fail if it has already been executed in this session')
        self.make=subprocess.run('make')
        
        #Offer path for json files    
        self.path_json = 'perf/data'            
        if os.path.exists(self.path_json): #If path already exists do nothing
            pass
        else:
            os.makedirs('perf/data') #Create path if not yet existing
        #Call function for check and gcheck    
        if self.cg=='check.exe':
            PerT.ex_check(self.ex_path)
        elif self.cg=='gcheck.exe':
            PerT.ex_gcheck(self.ex_path)
        else:
            print('Something went wrong. gcheck or check could not be executed')
            sys.exit()
                
    def ex_check(self,path): #In Progress
        print('Execute check.exe')
        
        self.maxevents = int(PerT.configvalues['variations']['numevents'])
        self.blocks =1
        self.thread_start = int(PerT.configvalues['variations']['threads_start'])
        self.thread_max = int(PerT.configvalues['variations']['threads_max'])       
        self.iters = None
        
        i =1
        for self.threads in [x*self.thread_start for x in[2**y for y in range(self.thread_max)]]:
            while self.blocks * self.threads < self.maxevents:
                self.iters = int(self.maxevents/(self.blocks * self.threads))
                
                print('%i %i %i' %(self.blocks, self.threads, self.iters))
                
                self.name = PerT.makename(self.blocks, self.threads, self.iters)
                
                subprocess.run(['./'+PerT.cg,'-j',str(self.blocks),str(self.threads),str(self.iters),self.name,str(i)])
                self.blocks = self.blocks *2
                i+=1
            self.blocks =1 
        
    def ex_gcheck(self,path):
        #execute gcheck.exe
        #shutil.rmtree('/tmp/andy/root/madgraph4gpu/epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/perf/data/')
        print('Execute gcheck.exe')
        
        self.maxevents = int(PerT.configvalues['variations']['numevents'])
        self.blocks =1
        self.thread_start = int(PerT.configvalues['variations']['threads_start'])
        self.thread_max = int(PerT.configvalues['variations']['threads_max'])       
        self.iters = None
        
        i =1
        for self.threads in [x*self.thread_start for x in[2**y for y in range(self.thread_max)]]:
            while self.blocks * self.threads < self.maxevents:
                self.iters = int(self.maxevents/(self.blocks * self.threads))
                
                print('%i %i %i' %(self.blocks, self.threads, self.iters))
                
                self.name = PerT.makename(self.blocks, self.threads, self.iters)
                
                subprocess.run(['./'+PerT.cg,'-j',str(self.blocks),str(self.threads),str(self.iters),self.name,str(i)])
                
                self.blocks = self.blocks *2
                i+=1
            self.blocks =1 
        
        
        #self.exe=subprocess.run(['./',PerT.cg,'-j','1024', '32','32','20210325', '1'],capture_output=True)
        
    def assemblingpath(self):
        print('Assembling path for execution\n')
        #Get values from configuration file to create execution path
        self.epoch  =   PerT.configvalues['runfiles']['epoch']+'/'
        self.lang   =   PerT.configvalues['runfiles']['lang']+'/'
        self.mumu   =   PerT.configvalues['runfiles']['mumu']+'/'
        self.sigma  =   PerT.configvalues['runfiles']['sigma']+'/'
        self.path   =   PerT.root_path+'/madgraph4gpu/'+self.epoch+self.lang+self.mumu+'SubProcesses/'+self.sigma
        return self.path
        
    def makename(self, blocks, threads, iters):
        
        self.name= str(blocks)+str(threads)+str(iters)
        return self.name
        
        
if __name__ == '__main__':
    
    PerT = PerformanceTesting()
    PerT.rootfolfdercreation()
    PerT.readConfig()
    PerT.CPU()
    PerT.GPU()
    PerT.cuda()
    PerT.check_gcheck()
    
    
