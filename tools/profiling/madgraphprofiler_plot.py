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
import json
import pandas as pd
import matplotlib.pyplot as plt
import re


class PerformanceTesting:
    root_path=None
    configvalues = {}
    cpu = None
    gpus = None
    
    
    def rootfolfdercreation(self):
        print('Creating root folder\n')
        PerT.profile_path = '/tmp/profile'
        PerT.root_path = '/tmp/profile/root'
        PerT.data_store= '/tmp/profile/data'
        if os.path.exists(PerT.root_path):
            print('root directory %s already existing.\n' % PerT.root_path)
            self.delete=input('Do you want to refresh the folder structure? y/n\n')
            if self.delete=='y':
                shutil.rmtree(PerT.profile_path)
                os.makedirs(PerT.root_path,0o777) #Create root directory
                os.makedirs(PerT.data_store,0o777) #Create data directory
            elif self.delete=='n':
                pass
            else:
                print('invalid input. Keep current folder structure')
                
        else:
        #Create root directory
            try:
                os.makedirs(PerT.root_path,0o777) #Create root directory
                os.makedirs(PerT.data_store,0o777) #Create data directory
            except OSError:
                print('Creation of the directory %s failed' % PerT.root_path)
                print('Root directory can neither be found nor created.\nExecution failed. ')
                sys.exit()
            else:
                print('Successfully created the directory %s\n' % PerT.root_path)    
                print('Successfully created the directory %s\n' % PerT.data_store)
                
    def cuda(self):
        print('NVIDIA has been detected.\nSet up cuda environment version 11.1')
        self.curdir=os.getcwd()
        
        
        self.curdir=str(os.getcwd())+'/GPU.sh'
        
        try:
            self.cuda = subprocess.run(['sh','.',self.curdir],capture_output=True)
            print('Call %s' %self.curdir)   
            print(self.cuda)
            #Check installed Cuda Version
        except OSError:
            print('\nMission failed\nGPU.sh not executed')
 
        # k = None
        # while k != 1:
            
        #     if os.getcwd() == PerT.root_path:
        #         self.clone=subprocess.run(['git','clone','https://github.com/madgraph5/madgraph4gpu.git'],capture_output=True)
        #         print(self.clone)
        #         k=1
        #     else:
        #         os.chdir(PerT.root_path)
        #         print('Changed directory to '+PerT.root_path)
    
    def gitclone(self):
        print('Execute git clone')
        k = None
        while k != 1:
            
            if os.getcwd() == PerT.root_path:
                print('Changed directory to '+PerT.root_path)
                self.clone=subprocess.run(['git','clone','https://github.com/madgraph5/madgraph4gpu.git'],capture_output=True)
                print(self.clone)
                k=1
            else:
                os.chdir(PerT.root_path)
                print('Changed directory to '+PerT.root_path)
        pass
    
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
                return self.gpu
        
    
                
    def check_gcheck(self,epoch):
        #Create folder structure according to configuration file and execute check or gcheck
        self.cg=PerT.configvalues['runfiles']['sys'] #Get check or gcheck
        #print('Assemble path according to configuration')
        print('Assembling path for execution\n')
        self.ex_path=PerT.assemblingpath(epoch) # Assemble path
        #returns something like: /tmp/profile/root/madgraph4gpu/epoch1/cuda/gg_tt/SubProcesses/P1_Sigma_sm_gg_ttx

        print('Hop in '+self.ex_path+'\n') 
        os.chdir(self.ex_path)
        print('Execute make file.\nIt may fail if it has already been executed in this session')
        input('make changes in mgOnGPUConfig.h and press return to continue\n')
        print('execute make\n')
        self.make=subprocess.run('make')
        
        #Offer path for json files    
        self.path_json = 'perf/data'
                  
        if os.path.exists(self.path_json): #If path already exists do nothing
            print('directory perf/data/ already exists. Datas will be overwritten' )
            
        else:
            print('create perf/data')
            os.makedirs('perf/data') #Create path if not yet existing
        #Call function for check and gcheck    
            print('Create path for data\n'+self.path)  
        if self.cg=='check.exe':
            PerT.ex_check(self.ex_path,epoch)
        elif self.cg=='gcheck.exe':
            PerT.ex_gcheck(self.ex_path,epoch)
        else:
            print('Something went wrong. gcheck or check could not be executed')
            sys.exit()
                
    def ex_check(self,path): #In Progress
        print('Execute check.exe')
        f=open('TestInProgress.txt','w+')
        self.maxevents = int(PerT.configvalues['variations']['numevents'])
        self.blocks =1
        self.thread_start = int(PerT.configvalues['variations']['threads_start'])
        self.thread_max = int(PerT.configvalues['variations']['threads_max'])       
        self.iters = None
        if os.path.exists(path+PerT.cg):
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
            os.remove('TestInProgress.txt')
        else:
            print('There is no gcheck')        
        
    def ex_gcheck(self,path,epoch):
        #execute gcheck.exe
        #shutil.rmtree('/tmp/andy/root/madgraph4gpu/epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/perf/data/')
        print('Execute gcheck.exe')
        f=open('TestInProgress.txt','w+')
        
        self.maxevents = int(PerT.configvalues['variations']['numevents'])
        self.blocks =1
        self.thread_start = int(PerT.configvalues['variations']['threads_start'])
        self.thread_max = int(PerT.configvalues['variations']['threads_max'])       
        self.iters = None
        if os.path.exists(path+PerT.cg):
            print(path+PerT.cg)
           
            i =1
            for self.threads in [x*self.thread_start for x in[2**y for y in range(self.thread_max)]]:
                while self.blocks * self.threads < self.maxevents:
                    #self.iters = int(self.maxevents/(self.blocks * self.threads))
                    self.iters=128
                    
                    print('%i %i %i' %(self.blocks, self.threads, self.iters))
                    
                    self.name = PerT.makename(self.blocks, self.threads, self.iters)
                    
                    subprocess.run(['./'+PerT.cg,'-j',str(self.blocks),str(self.threads),str(self.iters),self.name,str(i)])
                    # grid size 32 
                    # grid size 64
                        #32 2
                        #64 1
                    # 128     
                    self.blocks = self.blocks *2
                    i+=1
                self.blocks =1
            os.remove('TestInProgress.txt')
            print('Execution of ' +PerT.cg +' complete.\nCopy files to /tmp/profile/data')
            #Copy directory with json files to /tmp/profile/data
            shutil.copytree(path+'perf/data',PerT.data_store+'/'+
                            PerT.cg+'_'+
                            epoch+'_'+
                            PerT.configvalues['runfiles']['abstr_layer']+'_'+
                            PerT.configvalues['runfiles']['process'])
            print('try to copy from '+ path+ 'perf/data' +' '+'to '+ PerT.data_store)
            
        else:
            print('There is no gcheck')
            sys.exit()
        
        
        #self.exe=subprocess.run(['./',PerT.cg,'-j','1024', '32','32','20210325', '1'],capture_output=True)
        
    def assemblingpath(self,epoch):
        
        #Get values from configuration file to create execution path
        #self.epoch  =   PerT.configvalues['runfiles']['epoch']+'/'
        self.epoch = epoch+'/'
        self.lang   =   PerT.configvalues['runfiles']['abstr_layer']+'/'
        self.mumu   =   PerT.configvalues['runfiles']['process']+'/'
        self.sigma  =   PerT.configvalues['runfiles']['sigma']+'/'
        self.path   =   PerT.root_path+'/madgraph4gpu/'+self.epoch+self.lang+self.mumu+'SubProcesses/'+self.sigma
        return self.path
        
    def makename(self, blocks, threads, iters):
        
        self.name= str(blocks)+str(threads)+str(iters)
        return self.name
        


class Evaluation:
    
    list_results=[]     #List results
    Data=pd.DataFrame() #To store all results in one DataFrame
    plot_confi={}
    
    def readConfig(self):
        print('Including values from configuration file\n')
        self.config = configparser.ConfigParser()
        self.config.optionxform=str
        self.file='profileconfig.ini'
        self.config.read(self.file)
        self.sections=self.config.sections()
        for i in range(len(self.sections)): 
            Ev.plot_confi[self.sections[i]] = dict(self.config.items(self.sections[i]))
        print('Configuration completed\n')
        
    def loadfiles(self,directory):
        
        #Loads json files and returns a Dataframe
        #You can take a look into the shape of the returned df with the commnd df.info
        print('Load files directory:'+directory)
        os.chdir(directory)

        self.listfiles=os.listdir() #Get a list of all files stored in target directory
        
        '''Load files in DataFrame'''
        for i in range(len(self.listfiles)):
            
              self.file=open(self.listfiles[i],'r')        #Open file
              self.file_load=json.load(self.file)          #Load from file  
              Ev.list_results.append(self.file_load[0])    #Take the Data of the file that is stored in first list entry
              if i==0:                           #
                  Ev.Data=pd.DataFrame.from_dict(Ev.list_results[0],orient='index')    #initialize DataFrame
              Ev.Data[i]=pd.DataFrame.from_dict(Ev.list_results[i],orient='index')     #Add columns from list_results to df
        Data2=Ev.Data.T
        return Data2            
    
    def convertunits(self, dataframetoconvert):
        #Change the values and units from the values that are set to "on" in the config file
        #Return the a new dataframe and a list with elements to be plotted
        self.list_to_be_print=[] #list with elements to be plotted
   
        print('Change datatype of:\n')
        for item in Ev.plot_confi['plots']:
            
            
            if Ev.plot_confi['plots'][item] == 'on': #item = column name as a str
               
                self.list_to_be_print.append(item)
                print(item)
             
                for val in range(len(dataframetoconvert[item])):
                    temp_string=dataframetoconvert.loc[val,item]
                    #123456789.123456 sec^-1
                    #
                    #dataframetoconvert.loc[val,item]=int(temp_string[:temp_string.rfind('.')])
                    dataframetoconvert.loc[val,item]=float(re.findall(r'[\d.]+',temp_string)[0])
        dataframetoconvert['NumThreadsPerBlock*NumBlocksPerGrid']=dataframetoconvert['NumThreadsPerBlock']*dataframetoconvert['NumBlocksPerGrid']        
        dataframetoconvert['NumThreadsPerBlock*NumBlocksPerGrid']=dataframetoconvert['NumThreadsPerBlock*NumBlocksPerGrid'].astype('int')
        return dataframetoconvert,self.list_to_be_print

    def EventsPerSecMatrixElement(self):pass
        
        
    def plots(self,df,plotlist):
        
        for yaxis in plotlist:   
            print(yaxis)
            plt.style.use('ggplot')
            fig,ax = plt.subplots()
            
            fig.set_size_inches(18.5,10.5)
            fig.suptitle(yaxis,fontsize=50)
        
            plt.ylabel(yaxis,fontsize=30)
            plt.xlabel('NumThreadsPerBlock*NumBlocksPerGrid',fontsize=30)
            
            ax.set_yscale('log')
            ax.set_facecolor('darkgray')
        
            ax.set_xscale('log')
            ax.set_xticklabels(df['NumThreadsPerBlock*NumBlocksPerGrid'],rotation=75)
            
            #print(type(yaxis))
        
            for k in sorted(df['NumThreadsPerBlock'].unique(),reverse=True):
                    
                    d=df.loc[df['NumThreadsPerBlock']==k]
                    ax.scatter(
                        d['NumThreadsPerBlock*NumBlocksPerGrid'],
                        d[yaxis],
                        label=k,
                        s=k*6,
                        c=Ev.color(k),
                        alpha=Ev.al(k),                        
                        edgecolors='black'
                        )
                    ax.legend(loc='upper left',title='Threads Per Block',prop={'size':30})
                    plt.xticks(df['NumThreadsPerBlock*NumBlocksPerGrid'],fontsize=25)
                    plt.rcParams['ytick.labelsize']=25          
            plt.show()
            fig.savefig(PerT.assemblingpath()+yaxis)


    def al(self,value):
        if value == 32:
            self.alpha = .5
        elif value == 64:
            self.alpha = 0.75
        elif value == 128:
            self.alpha = 0.75
        elif value == 256:
            self.alpha=1
        return self.alpha
    
    def color(self,value):
        if value == 32:
            self.col = 'midnightblue'
        elif value == 64:
            self.col = 'sienna' 
        elif value == 128:
            self.col = 'aqua'
        elif value == 256:
            self.col='lightcyan'
        return self.col

    
if __name__ == '__main__':
    
    PerT = PerformanceTesting()
    Ev = Evaluation()
    
    PerT.readConfig()
    Ev.readConfig()
    
    PerT.rootfolfdercreation()
    
    
    PerT.CPU()
    PerT.GPU()
    PerT.gitclone()
    epoch = PerT.configvalues['runfiles']['epoch'].split(';')
    for epoch in epoch:
        
        PerT.check_gcheck(epoch)
    
    
    df_raw=Ev.loadfiles(PerT.assemblingpath()+'perf/data')
    #need path of datas
    df_adj_units, plotlist = Ev.convertunits(df_raw)
    Ev.plots(df_adj_units,plotlist)
    
    
