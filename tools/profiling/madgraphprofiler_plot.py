#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 10:25:46 2021

@author: andy
"""
from __future__ import print_function
import platform
import subprocess
import os
import configparser
import sys
import shutil
import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import re
import argparse
from datetime import datetime, date
from statistics import mean
import tempfile as tmp
from distutils.dir_util import copy_tree



class PerformanceTesting:
    root_path=None
    configvalues = {}
    cpu = None
    gpus = None
    product=str()
    today=date.today()
    date="{:%Y_%m_%d_%H_%M_%S}".format(datetime.now())
    throughput_dict={}
    storage_path=None
    
    def rootfolfdercreation(self,results_path):
        # results_path='profiling_results/'
        print('Creating profiling results directory')
        if os.path.exists(results_path):
            print('Result directory already existing.\nRefreshing result directory\n')
            shutil.rmtree(results_path)
            os.makedirs(results_path,0o777)
            os.chmod(results_path, 0o777)
        else:
            print('Create result directory\n'+os.getcwd()+results_path)
            os.makedirs(results_path,0o777)
            os.chmod(results_path, 0o777)
            
        
        print('Creating temporary root folder\n')
        #PerT.profile_path = '/tmp/profile'
        PerT.root_path = tmp.mkdtemp(prefix='profile_')
        os.chmod(PerT.root_path,0o777)
        print('created temporary directory'+PerT.root_path)
        PerT.data_store= PerT.root_path+'/data'
            
        # else:
        #     print('Create none tempoary profiler root directory')
        #     PerT.root_path = '/tmp/profile/'
        #     PerT.data_store= '/tmp/profile/data'
        #     if os.path.exists(PerT.root_path):
        #         print('root directory %s already existing.\n' % PerT.root_path)
        #         self.delete=input('Do you want to refresh the folder structure? y/n\n')
        #         if self.delete=='y':
        #             shutil.rmtree(PerT.root_path)
        #             os.makedirs(PerT.root_path,0o777)
        #             os.chmod(PerT.root_path, 0o777)
                    
        #             print('Refresh directory')
        #             os.makedirs(PerT.root_path,0o777) #Create root directory
        #             print('Change rights')
        #             os.chmod(PerT.root_path,0o777)
        #             per=os.stat(PerT.root_path)
        #             print(per)
        #             os.makedirs(PerT.data_store,0o777) #Create data directory
        #             os.chmod(PerT.data_store,0o777)
        #         elif self.delete=='n':
        #             pass
        #         else:
        #             print('invalid input. Keep current folder structure')
                
            # else:
            # #Create root directory
            #     try:
            #         os.makedirs(PerT.profile_path,0o777)
            #         os.chmod(PerT.profile_path, 0o777)
            #         os.makedirs(PerT.root_path,0o777) #Create root directory
            #         os.chmod(PerT.root_path,0o777)
            #         os.makedirs(PerT.data_store,0o777) #Create data directory
            #         os.chmod(PerT.data_store,0o777)
            #     except OSError:
            #         print('Creation of the directory %s failed' % PerT.root_path)
            #         print('Root directory can neither be found nor created.\nExecution failed. ')
            #         sys.exit()
            #     else:
            #         print('Successfully created the directory %s\n' % PerT.root_path)    
            #         print('Successfully created the directory %s\n' % PerT.data_store)
                
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
        with open('/proc/cpuinfo') as f:
            for line in f:
                if line.strip():
                    if line.rstrip('\n').startswith('model name'):
                        model_name=line.strip('\n').split(':')[1]
                        print(model_name)
    
    def GPU(self):
        print('Investigating GPU\n')
        self.gpu_hw_list=['NVIDIA','AMD','Asus','Intel','EVGA','Inno','Gigabyte','Zotac']
        
        
        self.getHW = subprocess.run(['lshw', '-C', 'display'], capture_output=True)
        self.display_list=str(self.getHW.stdout).split('\\n')
        print(self.display_list)
        #check which GPUs are installed and print it out
        self.found = False            
        for self.gpu in self.gpu_hw_list:
            
            vendor = [desc_item for desc_item in self.display_list if self.gpu in desc_item]
            if(len(vendor) == 0):
                print(f'Not found: {self.gpu}')
            else:
                print(f'Found:     {self.gpu}')
                PerT.gpus = self.gpu
                #return self.gpu
        # Get the product i.e. V100 for NVIDIA
        for self.gpu in self.gpu_hw_list:
            # Define index counter
            self.index =0
            for display_entry in self.display_list:
                if self.gpu in display_entry:
                    # Product is usually one list entry above the vendor entry
                    PerT.product=self.display_list[self.index-1]
                    PerT.product=PerT.product.replace('       product: ','')
                    print('Mounted card: '+PerT.product)
                self.index+=1
    
                
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
        input('Adjust settings in Makefile if needed.\nPress retun\n')
        print('execute make\n')
        
        self.make=subprocess.run('make')
        
        # Create json directory path for json files    
        self.path_json = 'perf/data'
        #shutil.rmtree(self.path_json)
           
        if os.path.exists(self.path_json): #If path already exists do nothing
            shutil.rmtree('perf/data')
            print('refresh perf/data' )
            os.makedirs('perf/data',0o777)
            os.chmod('perf/data',0o777)
        else:
            print('create perf/data')
            os.makedirs('perf/data',0o777) #Create path if not yet existing
            os.chmod('perf/data',0o777)
            # Call function for check and gcheck    
            print('Create path for data\n'+self.path)  
        
        
        if args.check != None: PerT.check_for_plot(epoch)
        os.chdir(self.ex_path)
        shutil.rmtree('perf/data')
        os.makedirs('perf/data',0o777)
        os.chmod('perf/data',0o777)
        if self.cg=='check.exe':
            PerT.ex_check(self.ex_path,epoch)
        elif self.cg=='gcheck.exe' or self.cg=='ccheck.exe':
            PerT.ex_gcheck(self.ex_path,epoch)
        else:
            print('Something went wrong. Neither gcheck nor check could not be executed')
            sys.exit()
    
    def check_for_plot(self,epoch):
        # If CPU values are required this function will bhe executet.
        # It executes check.exe and calculates the throughput to be added in the plot
        print('CPU was triggered by command line. Execute check.exe ' )
        print('WE ARE NOW IN '+os.getcwd())
        
        self.i=1
        while self.i <= args.check:
            print('Execute check.exe for plot run %s of ' %(self.i))
            subprocess.run(['./'+'check.exe','-j','32','32','128','0815',str(self.i)])
            self.i +=1
        # Open json files and take the data
        # Change directory 
        os.chdir(os.getcwd()+'/perf/data')
        # List items that will get plotted according to config file
        list_to_be_print=[item for item in Ev.plot_confi['plots'] if Ev.plot_confi['plots'][item] == 'on'] #item = column name as a str: #item = column name as a str]
        
        listoffiles = os.listdir()
        
        for print_item in list_to_be_print:
            print(list_to_be_print)
            temp_list=[]
            for i in range(len(listoffiles)):
                file = open(listoffiles[i])
                file_load=json.load(file)
                temp_string=file_load[0][print_item]               
                temp_list.append(int(float(re.findall(r'[\d.]+',temp_string)[0])))
            PerT.throughput_dict[print_item] = temp_list
            PerT.throughput_dict[print_item].append(int(mean(temp_list)))
        os.chdir(PerT.assemblingpath(epoch))
            #print('Avergae is %i ' %(throughput_dict[print_item][-1]))
        #print(throughput_dict)
    
        
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
                    
                    # subprocess.run(['./'+PerT.cg,'-j',str(self.blocks),str(self.threads),str(self.iters),self.name,str(i)])
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
                            PerT.configvalues['runfiles']['process']+'_'+
                            self.date+'/'+'jsonfiles')
            print('Copy files from '+ path+ 'perf/data' +' '+'to '+ PerT.data_store)
            
            #clean perf/data folder for next run
            self.folder=os.getcwd()+'/perf/data'
            for filename in os.listdir(self.folder):
                file_path = os.path.join(self.folder, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print('Failed to delete %s. Reason: %s' % (filename, e))
                            
                            
        else:
            print('check.exe can not be found')        
        
    def ex_gcheck(self,path,epoch):
        #execute gcheck.exe
        #shutil.rmtree('/tmp/andy/root/madgraph4gpu/epoch2/cuda/ee_mumu/SubProcesses/P1_Sigma_sm_epem_mupmum/perf/data/')
        print('Execute '+PerT.cg)
        f=open('TestInProgress.txt','w+')
        #today=date.today()
        self.date = PerT.date
        
        self.maxevents = int(PerT.configvalues['variations']['numevents'])
        self.blocks =1
        self.thread_start = int(PerT.configvalues['variations']['threads_start'])
        self.thread_max = int(PerT.configvalues['variations']['threads_max'])       
        self.iters = None
        if os.path.exists(path+PerT.cg):
            print('Now we are in '+path+PerT.cg)
            
            i =1
            for self.threads in [x*self.thread_start for x in[2**y for y in range(self.thread_max)]]:
                while self.blocks * self.threads < self.maxevents:
                    #self.iters = int(self.maxevents/(self.blocks * self.threads))
                    self.iters=128
                    
                    print('%i %i %i' %(self.blocks, self.threads, self.iters))
                    print('gridsize= '+str(self.blocks * self.threads))
                    self.name = PerT.makename(self.blocks, self.threads, self.iters)
                    #subprocess.run(['./'+PerT.cg,'-j',str(512),str(128),str(128),self.name,str(i)])
                    subprocess.run(['./'+PerT.cg,'-j',str(self.blocks),str(self.threads),str(self.iters),self.name,str(i)])
                       
                    self.blocks = self.blocks *2
                    i+=1
                self.blocks =1
            os.remove('TestInProgress.txt')
            print('Execution of ' +PerT.cg +' complete.\nCopy files to /tmp/profile/data')
            #Copy directory with json files to /tmp/profile/data
            PerT.storage_path = PerT.root_path+'/'+'data/'+PerT.cg+'_'+epoch+'_'+PerT.configvalues['runfiles']['abstr_layer']+'_'+PerT.configvalues['runfiles']['process']+'_'+self.date
            try:
                         
                shutil.copytree(path+'perf/data',PerT.storage_path+'/'+'jsonfiles')
                print('Copy files from '+ path+ 'perf/data' +' '+'to '+ PerT.storage_path+'/'+'jsonfiles')
                
                # shutil.copytree(path+'perf/data',profiler_path+'/'+results_dir+
                #                 PerT.cg+'_'+
                #                 epoch+'_'+
                #                 PerT.configvalues['runfiles']['abstr_layer']+'_'+
                #                 PerT.configvalues['runfiles']['process']+'_'+
                #                 self.date+'/'+'jsonfiles')
            except os.error:
                print('We are in error case. Let us see, if program keeps working')
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
        dataframetoconvert['gridsize']=dataframetoconvert['NumThreadsPerBlock']*dataframetoconvert['NumBlocksPerGrid']        
        dataframetoconvert['gridsize']=dataframetoconvert['gridsize'].astype('int')
        return dataframetoconvert,self.list_to_be_print

    def EventsPerSecMatrixElement(self):pass
        
        
    def plots(self,df,plotlist):
        plot_path=PerT.data_store+'/'+PerT.cg+'_'+epoch+'_'+\
            PerT.configvalues['runfiles']['abstr_layer']+'_'+\
            PerT.configvalues['runfiles']['process']+'_'+PerT.date+'/'+\
            'plots/'
        os.makedirs(plot_path,0o777) #Create root directory
        os.chmod(plot_path, 0o777)
                            
        for yaxis in plotlist:   
            print('Creation of '+ yaxis)
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
            fig.savefig(PerT.root_path+'/'+'data/'+
                                PerT.cg+'_'+
                                epoch+'_'+
                                PerT.configvalues['runfiles']['abstr_layer']+'_'+
                                PerT.configvalues['runfiles']['process']+'_'+
                                self.date+'/'+'plots'+yaxis+'.png')
            print('Plotting complete\n\n')
            
    def plots_(self,df,plotlist):
        plot_path=PerT.data_store+'/'+PerT.cg+'_'+epoch+'_'+\
            PerT.configvalues['runfiles']['abstr_layer']+'_'+\
            PerT.configvalues['runfiles']['process']+'_'+PerT.date+'/'+\
            'plots/'
        os.makedirs(plot_path,0o777)
        os.chmod(plot_path, 0o777)
        # plt.cla()
        # plt.clf()
                    
        for yaxis in plotlist:   
            print('Creation of '+ yaxis)
            
            
            fig = plt.figure()
            ax = plt.subplot()
            
            # Set figure size
            fig.set_size_inches(18,12)
           
            #remove frame
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            
            #enable grid
            plt.rcParams['grid.linestyle']=':'
            ax.yaxis.grid()
            
            #setup x-axis
            ax.set_xscale('log')
            plt.xticks(df['gridsize'])
            ax.set_xticklabels(df['gridsize'],rotation=75,fontsize=13.5)
            
            #setup y-axis
            if args.log:
                #logarithmic
                plt.ylim((min(df[yaxis])-min(df[yaxis])*0.30),max(df[yaxis])*3)
                ax.set_yscale('log')
                plt.yticks(size=13.5)    
                
            else:    
         
                plt.ylim(0,max(df[yaxis])*1.3) #linear scale
            
            
            #Labels and title
            
            plt.xlabel('Gridsize',fontsize=15)
            plt.ylabel('Troughput\n'+yaxis,fontsize=13.5)
            plt.title(yaxis,fontsize=15)
            if args.check != None:
                df['cpu']=PerT.throughput_dict[yaxis][-1]
                #label cpu value
                max_cpu=df['cpu'][0]
                length=len(str(max_cpu))-1
                label_maxima=str(round(max_cpu*10**-(length),3))+'e'+str(length)
                # Set new ylim
                #plt.ylim((min(df[yaxis])-min(df[yaxis])*0.30),max(df[yaxis])*3)
                #min((temp_df['NumIterations'].min(),temp_df['gridsize'].min()))
                plt.ylim((min(df[yaxis].min(),df['cpu'].min())*0.3),max(df[yaxis])*3)
                ax.plot(df['gridsize'],df['cpu'],'-ok',label='cpu: '+label_maxima,c='b')
            for k in sorted(df['NumThreadsPerBlock'].unique(),reverse=True):
                    
                    d=df.loc[df['NumThreadsPerBlock']==k]
                    ax.scatter(
                        d['gridsize'],
                        d[yaxis],
                        label=k,
                        s=k*2,
                        c=Ev.color(k),
                        alpha=Ev.al(k),                        
                        edgecolors='black'
                        )
                    ax.legend(loc='upper left',title='Threads Per Block',prop={'size':15})
                    plt.xticks(df['gridsize'])
             
            # Include TextBox with product information
            text_box=AnchoredText(PerT.gpus+' '+PerT.product, loc=1,prop=dict(size=15),frameon=False)
            plt.gca().add_artist(text_box)
            
            plt.show()
            # Andere Formate fuer LaTeX?????
            fig.savefig(PerT.storage_path+'/'+'plots/'+yaxis+'.png')
            #fig.savefig()
            print('Plotting complete\n\n')
        print('plots are ready to be investigated.\nSee path  %s' %plot_path)

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
    
    # Get command line options
    parser = argparse.ArgumentParser()
    # Optional parser if plots are wished
    parser.add_argument('-p','--plots',help='Generates performance plots and stores them in seperate directory',
                        action='store_true')
    # Optional parser if check.exe shall be executed 5 times
    parser.add_argument('-c','--check', type=int,help='Defines how many runs check.exe will have and calculates the average of the throughput to include it in the plot\n For example -c 5')
    
    # Optional parser for log scare
    parser.add_argument('-l','--log', help='Prints the results in log scale if enabled',action='store_true')
    
    # # Optional parser keep profiler directory or generate tmpdir
    # parser.add_argument('-t', '--tmpdir', help='If tmdir is enable the program deletes all directories and stores output in directory where program is stored',
    #                     action='store_true')
    
    args=parser.parse_args()
    
    PerT = PerformanceTesting()
    Ev = Evaluation()
    profiler_path = os.getcwd()
    
    
    # Check if config file is available
    if os.path.isfile('profileconfig.ini'):
        PerT.readConfig()
        Ev.readConfig()
    else:
        print('profileconfig.ini can not be found.\nPlease make sure that the configuration file is stored in the same directory as the profiler')
        sys.exit()
  
    results_dir=profiler_path+'/profiling_results/'
    PerT.rootfolfdercreation(results_dir)
    
    
    # Get CPU information
    PerT.CPU()
    
    # Get GPU information
    PerT.GPU()
    
    
    # Execute git clone from https://github.com/madgraph5/madgraph4gpu.git 
    PerT.gitclone()
    
    # Get a list of epoch that has to be investigated
    epoch = PerT.configvalues['runfiles']['epoch'].split(';')
    
    # Performe permutations and store results in /tmp/profile/data/
    for epoch in epoch:
        
        PerT.check_gcheck(epoch)
        print('directory before plot: '+os.getcwd())
        if args.plots:
            # Load json files in a DataFrame
            
            df_raw=Ev.loadfiles(PerT.assemblingpath(epoch)+'perf/data')
            
            
            # Get a list of plots from config file
           
            df_adj_units, plotlist = Ev.convertunits(df_raw)
           
            
            # Plot results
            
            Ev.plots_(df_adj_units,plotlist)
            
    print('Copy results to '+results_dir)
    copy_tree(PerT.data_store,results_dir)
    
    
   
    
