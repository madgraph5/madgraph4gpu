#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 09:59:03 2021

@author: andy
"""
import json
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style
#import matplotlib.image as mpimg
#from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
import seaborn as sns
import configparser
import re
import math




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
        
    def load_df(self,path):
        df_dict={}
        
        os.chdir(path)
        listfolders = os.listdir()
        
        for datafolder in listfolders:
            os.chdir(path+'/'+datafolder)   #Jump in datafolder
            df_dict[datafolder]=pd.DataFrame()
            Data=pd.DataFrame()
            list_results =[]
            i=0
            while i in range(len(os.listdir())):
                for files in os.listdir():
                    file=open(files,'r')            #open json file
                    file_load=json.load(file)       #load json file
                    list_results.append(file_load[0])  #fill list with results of the loaded json file
                    #Data.append(list_results[i],ignore_index=True)
                    if i==0:
                        Data=pd.DataFrame.from_dict(list_results[0],orient='index')
                    Data[i]=pd.DataFrame.from_dict(list_results[i],orient='index')
                    
                 
                   
                    i+=1
                df_dict[datafolder] = Data.T
    
        return df_dict
        
    def loadfiles(self,directory):
        
        #Loads json files and returns a Dataframe
        #You can take a look into the shape of the returned df with the commnd df.info
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
                    if type(dataframetoconvert.loc[val,item]) !='float':
                        #dataframetoconvert.loc[val,item]=int(temp_string[:temp_string.rfind('.')])
                        dataframetoconvert.loc[val,item]=float(re.findall(r'[\d.]+',temp_string)[0])
                    else:
                        dataframetoconvert.loc[val,item]=(re.findall(r'[\d.]+',temp_string)[0])
                    
        dataframetoconvert['NumThreadsPerBlock*NumBlocksPerGrid']=dataframetoconvert['NumThreadsPerBlock']*dataframetoconvert['NumBlocksPerGrid']        
        dataframetoconvert['NumThreadsPerBlock*NumBlocksPerGrid']=dataframetoconvert['NumThreadsPerBlock*NumBlocksPerGrid'].astype('int')
        return dataframetoconvert,self.list_to_be_print

    
    def convertunits_2(self):
        
        dataframes_threated={}
    
        for df in dataframes:
            # reduce dataframes
            # iter_1024 = dataframes['iteration_fix_1024'][['NumIterations','NumThreadsPerBlock', 'NumBlocksPerGrid','EvtsPerSec[MatrixElems] (3)']]
            temp_df =pd.DataFrame()
            temp_df = dataframes[df][['NumIterations','NumThreadsPerBlock', 'NumBlocksPerGrid',
                                      'EvtsPerSec[MatrixElems] (3)','EvtsPerSec[Rnd+Rmb+ME](123)',
                                      'EvtsPerSec[Rmb+ME] (23)']]
            
            columns_to_convert = ['EvtsPerSec[MatrixElems] (3)','EvtsPerSec[Rnd+Rmb+ME](123)',
                                      'EvtsPerSec[Rmb+ME] (23)']
            for column in columns_to_convert:
                
                for val in range(len(temp_df[column])):
                    temp_string=temp_df.loc[val,column]
                    
                    temp_df.loc[val,column]=int(float(re.findall(r'[\d.]+',temp_string)[0]))
            temp_df['gridsize']=temp_df['NumThreadsPerBlock']*temp_df['NumBlocksPerGrid']
            temp_df['gridsize']=temp_df['gridsize'].astype(int)
            dataframes_threated[df]=temp_df
            
        return dataframes_threated
        
    def plots(self,df,plotlist):
        plt.cla()
        plt.clf()
        
        for yaxis in plotlist:   
           
            fig = plt.figure()
            ax = plt.subplot()
            
            #set figure size
            fig.set_size_inches(18,10)
            
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
            #logarithmic
            plt.ylim((min(df[yaxis])-min(df[yaxis])*0.30),max(df[yaxis])*3)
            ax.set_yscale('log')
            plt.yticks(size=13.5)
            #linear scale
            #plt.ylim(0,max(df[yaxis])*1.3)
            
            #Labels and titel
            plt.xlabel('Gridsize',fontsize=15)
            plt.ylabel('Troughput\n'+yaxis,fontsize=13.5)
            plt.title(yaxis,fontsize=15)
        
            # plt.ylabel(yaxis,fontsize=30)
            # plt.xlabel('NumThreadsPerBlock*NumBlocksPerGrid',fontsize=30)
            # plt.rcParams['axes.facecolor']='whitesmoke'
            # plt.rcParams['ytick.labelsize']=25
            # #ax.set_yscale('log')
            # ax.set_facecolor('darkgray')        
            # ax.set_xscale('log')#(1)
            # ax.set_xticklabels(df['gridsize'],rotation=75)
            # ax.set_facecolor('silver')
            # #print(type(yaxis))

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
                    #(1)
            #plt.rcParams['legend.title_fontsize']='large'
            #plt.text(16400, 250000, 'Here we have space for some\nfurther information like:\n\nCuda\nepoch2\ngg_ttgg',fontsize=25)
            plt.show()
            fig.savefig('/home/andy/cernbox/data/raw/data'+'epoch2_ee_mumu_gcheck_float'+yaxis)

    def data_compare(self,df_dict,compare_list,stat):
        #This function takes the dictinary of data frames and plots the selected df from the list
        #stat means the statistic value like "max", "min","avg"
        
        
        
        df_to_be_plotted=pd.DataFrame()
        gridsizelist=sorted(dataframes_conv[list(dataframes_conv.keys())[0]]['gridsize'].astype(int).unique())
        
        #DataFrame treatment: 
        
        for gs in gridsizelist: #get the gridsize values
            for df in compare_list:
                #self.temp_df=test_i=test[(test['gridsize']==1024)]
                temp_df=df_dict[df][(df_dict[df]['gridsize']==gs)]
                temp_df['process']=df
                if temp_df.empty:
                    pass
                else:
                    df_to_be_plotted=df_to_be_plotted.append(temp_df[(temp_df['EvtsPerSec[MatrixElems] (3)']
                                                               ==eval(stat)(temp_df['EvtsPerSec[MatrixElems] (3)']))])
                    df_to_be_plotted=df_to_be_plotted.astype({'gridsize':int}) 
        
        
        fig,ax = plt.subplots()
        
        
        ax.set_xscale('log')
        plt.ylim((0,1e9))
        plt.xticks(df_to_be_plotted['gridsize'])
        plt.suptitle('space for a smart title')
        plt.rcParams['ytick.labelsize']=10
        max_tp=max(df_to_be_plotted['EvtsPerSec[MatrixElems] (3)'])
        max_tp=round(max_tp/(1e9),3)   
        plt.text(64, 200000000,'max Throughput\n'+str(max_tp)+'E9 s-1')
        #plt.text(32, 900000000, 'max throughput from profiler: %i s^-1'%epoch2_ee_mumu_double)
        #ax.legend(loc='upper left',title='yolo')
        
        
        sns.scatterplot(data=df_to_be_plotted,x='gridsize',y='EvtsPerSec[MatrixElems] (3)',
                        hue='process',legend='full',style='process',palette='deep')
        ax.set_xticklabels(df_to_be_plotted['gridsize'],rotation=75)
        
        
        return df_to_be_plotted
        
    
    def data_compare2(self,df_dict,compare_list):
        #Takes a dictionary with dataframes and plots it in the same scatter plot
        
        fig = plt.figure()
        ax1 = plt.subplot()
        
        fig.set_size_inches(12.5,5.5)
        
        #remove frame
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        
        #enable grid
        plt.rcParams['grid.linestyle']=':'
        plt.grid()
        
        #setup x axis
        ax1.set_xscale('log')
        plt.xticks(df_dict[list(df_dict.keys())[0]]['gridsize'])
        ax1.set_xticklabels(df_dict[list(df_dict.keys())[0]]['gridsize'],rotation=75)
        
        #setup y axis
        #get maximum value of all df for ylim
        max_y = [max(df_dict[df]['EvtsPerSec[MatrixElems] (3)']) for df in df_dict]
        #plt.ylim(-0.1*10**9,max(max_y)*1.3)
        plt.ylim(10**5,max(max_y)*10)
        ax1.set_yscale('log')
        
        #Add labels and title
        plt.ylabel('Throughput\nMatrix Elements [s-1]')
        plt.xlabel('Gridsize')
        plt.title('Cuda throughput for ee_mumu on NVIDIA T4\n')
        
        #Change colormap. More info here https://matplotlib.org/stable/tutorials/colors/colormaps.html 
        cmap=plt.get_cmap('Set1')
        
        i=1
        for data in compare_list:
            #Get maximum values for each dataset
            maxima_y=max(df_dict[data]['EvtsPerSec[MatrixElems] (3)'])
            maxima_x=df_dict[data].loc[df_dict[data]['EvtsPerSec[MatrixElems] (3)']==maxima_y,'gridsize'].item()
            
            #label maximum values
            length=len(str(maxima_y))-1
            label_maximas=str(round(maxima_y*10**-(length),3))+'e'+str(length)
            
            #plot datasets
            ax1.scatter(df_dict[data]['gridsize'].to_list(),df_dict[data]['EvtsPerSec[MatrixElems] (3)'].to_list(),
                        label=data+ ' (max = %s)'%label_maximas,
                        color=cmap(i),
                        s=150,alpha=0.9)
            #Get next cmap color
            i+=1 
            
            #plot max values
            ax1.scatter(maxima_x,maxima_y,c='r',marker='o',s=50)
            
            
            
        ax1.legend(loc='upper left')  
        
        plt.show()
        
    def dataframes_statistical_transfomation(self,df_dict,stat):
        #This functions takes a dictionary of dataframes and returns a dictionary with dataframes
        #where the dataframes contain only the gridsize values according to the statistical expression
        #which was given to the function like max, min, etc..
        df_dict_to_return={}
        
        #Get a list of gridsize
        gridsize = sorted(df_dict[list(df_dict.keys())[0]]['gridsize'].astype(int).unique())
        for df in df_dict:
            df_dict_to_return[df]=pd.DataFrame()
            for gs in gridsize:
                temp_df=df_dict[df][(df_dict[df]['gridsize']==gs)]
        
                if temp_df.empty:
                    pass
                else:
                    df_dict_to_return[df]=df_dict_to_return[df].append(temp_df[(temp_df['EvtsPerSec[MatrixElems] (3)']
                                                                      ==eval(stat)(temp_df['EvtsPerSec[MatrixElems] (3)']))])
                    df_dict_to_return[df]=df_dict_to_return[df].astype({'gridsize':int}) 
        return df_dict_to_return

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

    
if __name__=='__main__':
    
    Ev = Evaluation()
    Ev.readConfig()
    #logo=mpimg.imread('/home/andy/cernbox/Madgraph/profiler/Logo/Logo_CERN.png')
    #imagebox=OffsetImage(logo)
    path='/home/andy/cernbox/data/Andrea'
    dataframes=Ev.load_df(path) #returns a directory that contains df for all data given in the path
    plotlist= [item for item in Ev.plot_confi['plots']if Ev.plot_confi['plots'][item] == 'on']
    
    

    dataframes_conv=Ev.convertunits_2() #returns a df directory with converted units
    dataframes_statisical=Ev.dataframes_statistical_transfomation(dataframes_conv,'max')
  
    '''
    Ev.plots(dataframes_conv['gcheck.exe_epoch1_cuda_ee_mumu_double'],plotlist)
    '''
    #max(df_adj_units['EvtsPerSec[MatrixElems] (3)'])
    # To be done
    list_to_compare=['gcheck.exe_epoch1_cuda_ee_mumu_float','gcheck.exe_epoch1_cuda_ee_mumu_double']
    #test_df=Ev.data_compare(dataframes_conv,list_to_compare,'max')
    
    Ev.data_compare2(dataframes_statisical,list_to_compare)
    
    dataframes_statisical[list(dataframes_statisical.keys())[0]]
    dataframes_statisical[list(dataframes_statisical.keys())[0]]['gridsize']
    dataframes_statisical['gcheck.exe_epoch1_cuda_ee_mumu_float'].dtypes
    
    
    