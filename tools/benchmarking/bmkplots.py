#!/usr/bin/env python3

import os, sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from subprocess import Popen

# Adapt plots to screen size
###plots_smallscreen=False # desktop
plots_smallscreen=True # laptop
if plots_smallscreen:
    plots_figsize=(6,3)
    plots_ftitlesize=9
    plots_txtsize=9
    plots_labelsize=9
    plots_legendsize=8
else:
    plots_figsize=(10,5)
    plots_ftitlesize=12
    plots_txtsize=15
    plots_legendsize=12
    plots_labelsize=12

#---------------------------------------

def loadOneRun( workdir, debug=False ):
    ###debug=True
    import json
    run_file = workdir + '/mg5amc-madgraph4gpu-2022_summary.json'
    print( 'Loading Run', run_file )
    run_dict = json.load( open( run_file ) )
    run_info = run_dict['run_info']
    run_info['app_version'] = run_dict['app']['version']
    if debug : print( run_info )
    run_scores = run_dict['report']['wl-scores']
    if debug : print( run_scores )
    return run_info, run_scores

#---------------------------------------

def loadCppRunSet( runsetdir, evtmatch='-e001', debug=False ):
    ###debug=True
    if not os.path.isdir( runsetdir ):
        print( 'Unknown directory', runsetdir )
        return
    # Go through all runs in runsetdir and fill cpprunset_scores[njob,nthr]
    cpprunset_scores = {}
    print( 'Loading runs in RunSetDir', runsetdir, 'for events', evtmatch )
    for d in sorted( os.listdir( runsetdir ) ) :
        if d.startswith( 'sa-cpp' ) and d.endswith( evtmatch ) and 'png' not in d : # e.g. sa-cpp-j004-t001-e001
            rundir = runsetdir + '/' + d
            if 'SKIP' in os.listdir( rundir ):
                print( 'WARNING! %s contains file SKIP and will be skipped'%rundir )
                continue
            dl = d.split( '-' )
            njob = int( dl[-3][-3:] )
            nthr = int( dl[-2][-3:] )
            nevt = int( dl[-1][-3:] )
            if debug : print( '\nRunDir=%30s %3i %3i %3i'%( rundir, njob, nthr, nevt ) )
            run_info, run_scores = loadOneRun( rundir )
            if len(run_scores) == 0 :
                print( 'WARNING! %s contains 0 scores and will be skipped'%rundir )
                continue
            njobkey, nthrkey, nevtkey = 'copies', 'threads_per_copy', 'events_per_thread'
            assert njob == run_info[njobkey], 'njob mismatch %i != %i'%( njob, run_info[njobkey] )
            assert nthr == run_info[nthrkey], 'nthr mismatch %i != %i'%( nthr, run_info[nthrkey] )
            assert nevt == run_info[nevtkey], 'nevt mismatch %i != %i'%( nevt, run_info[nevtkey] )
            cpprunset_scores[njob,nthr] = run_scores
    return cpprunset_scores

#---------------------------------------

def dumpCppScoresOneKey( cpprunset_scores, score_key, debug=False ):
    ###debug=True
    print( '\nSCORES[\'%s\']:'%score_key )
    score_key_none = score_key[:-4]+'none'
    njobs = set( [njobnthr[0] for njobnthr in cpprunset_scores] ) # use set(list) to get unique keys
    nthrs = set( [njobnthr[1] for njobnthr in cpprunset_scores] ) # use set(list) to get unique keys
    print( '%4s %4s %12s    %9s %12s %12s %16s'%( 'njob', 'nthr', 'Score', 'njob*nthr', 'S/S[1,1]', 'S/S-none', 'S/S-none[1,1]' ) )
    assert (1,1) in cpprunset_scores, 'no scores found for njob==1 and nthr==1?'
    tput1 = cpprunset_scores[1,1][score_key]
    tput1none = cpprunset_scores[1,1][score_key_none]
    for nthr in sorted(nthrs):
        for njob in sorted(njobs):
            if (njob,nthr) not in cpprunset_scores: continue
            tput = cpprunset_scores[njob,nthr][score_key]
            tputnone = cpprunset_scores[njob,nthr][score_key_none]
            print( '%4d %4d %12.6f    %9d %12.6f %12.6f %16.6f'%
                   ( njob, nthr, tput, njob*nthr, tput / tput1, tput / tputnone, tput / tput1none ) )

#---------------------------------------

def getSortedMatchingKeys( cpprunset_scores, keymatch=None, debug=False ):
    ###debug=True
    keys = []
    for njobnthr in cpprunset_scores : keys += list( cpprunset_scores[njobnthr].keys() )
    keys = set( keys ) # use set(list) to get unique keys
    if keymatch is not None: keys = [ key for key in keys if keymatch in key ]
    def sortableSimdKey( key ): # use keys sortable in this order: none, sse4, avx2, 512y, 512z, best
        key = key.replace( '-none', '-simd0' )
        key = key.replace( '-sse4', '-simd1' )
        key = key.replace( '-avx2', '-simd2' )
        key = key.replace( '-512y', '-simd3' )
        key = key.replace( '-512z', '-simd4' )
        key = key.replace( '-best', '-simd5' )
        return key
    keys2 = [ sortableSimdKey( key ) for key in keys ]
    keys = [ key for _, key in sorted( zip( keys2, keys ) ) ] # https://stackoverflow.com/a/6618543
    return keys

#---------------------------------------

def dumpCppScoresAllKeys( cpprunset_scores, keymatch=None, debug=False ):
    keys = getSortedMatchingKeys( cpprunset_scores, keymatch, debug )
    for key in keys : dumpCppScoresOneKey( cpprunset_scores, key )

#---------------------------------------

def loadCudaRunSet( runsetdir, evtmatch='-e001', gputhreads='gt00256', debug=False ):
    ###debug=True
    if not os.path.isdir( runsetdir ):
        print( 'Unknown directory', runsetdir )
        return
    # Go through all runs in runsetdir and fill cudarunset_scores[ngbl,ngth,njob,nthr]
    node, xht, ftitle = getNodeFeatures( runsetdir )
    cudarunset_scores = {}
    print( 'Loading runs in RunSetDir', runsetdir, 'for events', evtmatch )
    for d in sorted( os.listdir( runsetdir ) ) :
        if d.startswith( 'sa-cuda' ) and gputhreads in d and d.endswith( evtmatch )\
               and 'png' not in d : # e.g. sa-cuda-gb00064-gt00256-j001-t001-e0100
            rundir = runsetdir + '/' + d
            if 'SKIP' in os.listdir( rundir ):
                print( 'WARNING! %s contains file SKIP and will be skipped'%rundir )
                continue
            dl = d.split( '-' )
            ngbl = int( dl[-5][-5:] )
            ngth = int( dl[-4][-5:] )
            njob = int( dl[-3][-3:] )
            nthr = int( dl[-2][-3:] )
            nevt = int( dl[-1][-4:] )
            if debug : print( '\nRunDir=%30s %5i %5i %3i %3i %4i'%( rundir, ngbl, ngth, njob, nthr, nevt ) )
            run_info, run_scores = loadOneRun( rundir )
            if len(run_scores) == 0 :
                print( 'WARNING! %s contains 0 scores and will be skipped'%rundir )
                continue
            njobkey, nthrkey, nevtkey, nextkey = 'copies', 'threads_per_copy', 'events_per_thread', 'extra_arguments'
            assert '-p%d,%d,1'%(ngbl,ngth) in run_info[nextkey]
            assert njob == run_info[njobkey], 'njob mismatch %i != %i'%( njob, run_info[njobkey] )
            assert nthr == run_info[nthrkey], 'nthr mismatch %i != %i'%( nthr, run_info[nthrkey] )
            assert nevt == run_info[nevtkey], 'nevt mismatch %i != %i'%( nevt, run_info[nevtkey] )
            #if abs(xht) < njob:
            #    print( 'WARNING! %s contains overcommit data (njob=%d > nproc=%d) and will be skipped'%(rundir,njob,abs(xht)) )
            #    continue
            cudarunset_scores[ngbl,ngth,njob,nthr] = run_scores
    return cudarunset_scores

#---------------------------------------

def dumpCudaScoresOneKey( cudarunset_scores, score_key, debug=False ):
    ###debug=True
    print( '\nSCORES[\'%s\']:'%score_key )
    ngbls = set( [njobnthr[0] for njobnthr in cudarunset_scores] ) # use set(list) to get unique keys
    ngths = set( [njobnthr[1] for njobnthr in cudarunset_scores] ) # use set(list) to get unique keys
    njobs = set( [njobnthr[2] for njobnthr in cudarunset_scores] ) # use set(list) to get unique keys
    nthrs = set( [njobnthr[3] for njobnthr in cudarunset_scores] ) # use set(list) to get unique keys
    ngblmin = min( ngbls )
    ngthmin = min( ngths )
    assert (ngblmin,ngthmin,1,1) in cudarunset_scores, 'no scores found for (ngbl,ngth,njob,nthr)==(%d,%d,1,1)?'%(ngblmin,ngthmin)
    tputmin = cudarunset_scores[ngblmin,ngthmin,1,1][score_key]
    print( '%5s %5s %5s %5s %12s %20s %20s'%( 'nthr', 'njob', 'ngbl', 'ngth', 'Score', 'ngbl*ngth*njob*nthr', 'S/S[%d,%d,1,1]'%(ngblmin,ngthmin) ) )
    for nthr in sorted(nthrs):
        for njob in sorted(njobs):
            for ngbl in sorted(ngbls):
                for ngth in sorted(ngths):
                    if (ngbl,ngth,njob,nthr) not in cudarunset_scores: continue
                    tput = cudarunset_scores[ngbl,ngth,njob,nthr][score_key]
                    print( '%5d %5d %5d %5d %12.6f %20d %20.6f'%
                           ( nthr, njob, ngbl, ngth, tput, ngbl*ngth*njob*nthr, tput/tputmin ) )

#---------------------------------------

# Get node-dependent features
def getNodeFeatures( workdir ):
    if workdir == 'BMK-pmpe04' :
        node='pmpe04'
        xht=16
        ftitle='check.exe scalability on pmpe04 (2x 8-core 2.4GHz Xeon E5-2630 v3 with 2x HT)' # lscpu
    elif workdir == 'BMK-itscrd70' :
        node='itscrd70'
        xht=-4
        ftitle='check.exe scalability on itscrd70 (1x 4-core 2.1GHz Xeon Silver 4216 without HT)' # lscpu
    elif workdir == 'BMK-itscrd70-cuda' :
        node='itscrd70'
        xht=-4
        ftitle='gcheck.exe scalability on itscrd70 (V100 GPU on a 4-thread Xeon Silver 4216 CPU)' # lscpu
    elif workdir == 'BMK-jwlogin08' :
        node='jwlogin08'
        xht=40
        ftitle='check.exe scalability on jwlogin08 (2x 20-core 2.4GHz Xeon Gold 6148 with 2x HT)' # lscpu
    elif workdir == 'BMK-bmk6130' :
        node='bmk6130'
        xht=32
        ###ftitle='check.exe scalability on bmk-ironic-0731f1ce3b (2x 16-core 2.1GHz Xeon Gold 6130 with 2x HT)' # lscpu
        ftitle='check.exe scalability on "bmk6130" (2x 16-core 2.1GHz Xeon Gold 6130 with 2x HT)' # lscpu
    else:
        print( 'ERROR! Unknown workdir', workdir )
        sys.exit(-1)
    return node, xht, ftitle

#---------------------------------------

# Compare various curves in ST plots
# NB: xht>0 is the number of physical cores before 2xHT, xht<0 is the number of physical cores without HT
def axesCppST( ax, cpprunset_scores, keymatch=None, bestonly=False, abstput=True, ylog=False, xht=None, debug=False ):
    # Prepare axes labels
    ax.set_xlabel('Level of parallelism (number of ST jobs)', size=plots_labelsize )
    if abstput:
        ax.set_ylabel('Throughput (E6 events per second)', size=plots_labelsize )
    else:
        ax.set_ylabel('Throughput ratio to 1 no-SIMD job', size=plots_labelsize )
        ax.grid()
    if ylog: ax.set_yscale( 'log' )
    # Add one curve per matching score key
    xmax = 0
    ymax = 0
    keys = getSortedMatchingKeys( cpprunset_scores, keymatch, debug )
    if bestonly : keys = [ key for key in keys if 'best' in key ]
    elif 'best' not in keymatch: keys = [ key for key in keys if 'best' not in key ]
    for score_key in keys :
        score_key_none = score_key[:-4]+'none'
        njobs = set( [njobnthr[0] for njobnthr in cpprunset_scores] ) # use set(list) to get unique keys
        nthrs = set( [njobnthr[1] for njobnthr in cpprunset_scores] ) # use set(list) to get unique keys
        assert (1,1) in cpprunset_scores, 'no scores found for njob==1 and nthr==1?'
        ###tput1 = cpprunset_scores[1,1][score_key]
        tput1none = cpprunset_scores[1,1][score_key_none]
        # Prepare x-axis and y-axis lists
        xvals = []
        yvals = []
        for nthr in sorted(nthrs):
            for njob in sorted(njobs):
                if (njob,nthr) not in cpprunset_scores: continue
                xval = nthr*njob # 'npar' level of parallelism
                tput = cpprunset_scores[njob,nthr][score_key]
                ###tputnone = cpprunset_scores[njob,nthr][score_key_none]
                xvals.append( xval )
                if abstput: yvals.append( tput )
                else: yvals.append( tput / tput1none )
        xmax = max( xmax, max( xvals ) )
        ymax = max( ymax, max( yvals ) )
        # Add curve of y vs x
        p = ax.plot( xvals, yvals, marker='o', label=score_key )
    # Decorate axes
    loc = 'lower right'
    ax.legend( loc=loc, fontsize=plots_legendsize )
    xmin = 0
    if xht is None:
        hasht = False
        hasovercommit = False
    elif xht > 0:
        hasht = True
        hasovercommit = ( xmax > 2*xht )
    else:
        hasht = False
        xht = -xht
        hasovercommit = ( xmax > xht )
    xmax *= 1.8
    if ylog:
        ymin = 0.001
        ymax *= 12
        ytxt = 0.5 * ymax
    else:
        ymin = 0
        ymax *= 1.2
        ytxt = 0.92 * ymax
    ax.axis( [xmin, xmax, ymin, ymax] )
    if xht is not None :
        ax.axvline( xht, color='black', ls=':' )
        if hasht: ax.axvline( xht*2, color='black', ls='-.' )
        ax.text( xht/2, ytxt, 'No HT', ha='center', va='center', size=plots_txtsize )
        if hasht: ax.text( xht*3/2, ytxt, '2x HT', ha='center', va='center', size=plots_txtsize )
        if hasovercommit:
            xtxtover = ( xmax/2+xht if hasht else xmax/2+xht/2 )
            ax.text( xtxtover, ytxt, 'Overcommit', ha='center', va='center', size=plots_txtsize )

# Create a figure with a single plot
def plotCppST( workdir, keymatch=None, abstput=True, ylog=False, evtmatch='-e001', debug=False ):
    cpprunset_scores = loadCppRunSet( workdir, evtmatch=evtmatch )
    node, xht, ftitle = getNodeFeatures( workdir )
    pngpath = workdir + '/' + node + evtmatch + '-all-' + keymatch + '.png'
    # Create figure with one plot
    fig = plt.figure( figsize=plots_figsize )
    ax1 = fig.add_subplot( 111 )
    # Fill the plot in the figure
    axesCppST( ax1, cpprunset_scores, keymatch=keymatch, ylog=ylog, xht=xht, abstput=abstput, debug=debug )
    if ftitle is not None:
        if evtmatch.startswith( '-e' ) and evtmatch[2:].isdigit(): ftitle += ' for %d cycles'%int(evtmatch[2:])
        fig.suptitle( ftitle, size=plots_ftitlesize )
    # Save and show the figure
    fig.savefig( pngpath, format='png', bbox_inches="tight" )
    Popen( ['display', '-geometry', '+50+50', pngpath] )
    ###Popen( ['display', '-geometry', '+50+50', '-resize', '800', pngpath] )
    print( 'Plot successfully saved on', pngpath )

# Create a figure with two plots per process, absolute and normalized tput
def plotCppOneProcess2( workdir, oneprocess, keymatch, bestonly=False, evtmatch='-e001', debug=False ):
    cpprunset_scores = loadCppRunSet( workdir, evtmatch=evtmatch )
    node, xht, ftitle = getNodeFeatures( workdir )
    # One process or all processes?
    if oneprocess is not None: processes = [ oneprocess ]
    else: processes = [ 'eemumu', 'ggtt', 'ggttg', 'ggttgg' ]
    pngpath = workdir + '/' + node + evtmatch + '-' + \
              ( 'all' if oneprocess is None else oneprocess ) + '-' + keymatch + '.png'
    # Create figure with two plots per process
    fig = plt.figure( figsize = ( plots_figsize[0]*2, plots_figsize[1]*len(processes) ) )
    # Add two plots per process
    idx1 = len(processes)*100 + 21
    idx2 = len(processes)*100 + 22
    for process in processes:
        ###print( idx1, idx2 )
        fullkeymatch = process + '-' + keymatch
        ax1 = fig.add_subplot( idx1 )
        axesCppST( ax1, cpprunset_scores, keymatch=fullkeymatch, bestonly=bestonly, \
                   ylog=False, xht=xht, abstput=True, debug=debug )
        ax2 = fig.add_subplot( idx2 )
        axesCppST( ax2, cpprunset_scores, keymatch=fullkeymatch, bestonly=bestonly, \
                   ylog=False, xht=xht, abstput=False, debug=debug )
        idx1 += 2
        idx2 += 2
    # Save and show the figure
    if ftitle is not None:
        if oneprocess is not None: ftitle = oneprocess + ' ' + ftitle
        if evtmatch.startswith( '-e' ) and evtmatch[2:].isdigit(): ftitle += ' for %d cycles'%int(evtmatch[2:])
        fig.suptitle( ftitle, size=plots_ftitlesize )
    fig.set_tight_layout( True )
    fig.savefig( pngpath, format='png', bbox_inches="tight" )
    if oneprocess is not None:
        Popen( ['display', '-geometry', '+50+50', pngpath] )
    else:
        Popen( ['display', '-geometry', '+50+50', '-resize', '600', pngpath] )
    print( 'Plot successfully saved on', pngpath )

# Create a figure with one plots per process, with inl0 and inl1 best
def plotCppProcessesInl( workdir, keymatch, evtmatch='-e001', debug=False ):
    cpprunset_scores = loadCppRunSet( workdir, evtmatch=evtmatch )
    node, xht, ftitle = getNodeFeatures( workdir )
    # One process or all processes?
    processes = [ 'eemumu', 'ggtt', 'ggttg', 'ggttgg' ]
    pngpath = workdir + '/' + node + evtmatch + '-all-' + keymatch + '.png'
    # Create figure with two plots per process
    fig = plt.figure( figsize = ( plots_figsize[0], plots_figsize[1]*len(processes) ) )
    # Add two plots per process
    idx1 = len(processes)*100 + 11
    for process in processes:
        fullkeymatch = process + '-' + keymatch
        ax1 = fig.add_subplot( idx1 )
        axesCppST( ax1, cpprunset_scores, keymatch=fullkeymatch, bestonly=True, \
                   ylog=False, xht=xht, abstput=True, debug=debug )
        idx1 += 1
    # Save and show the figure
    if ftitle is not None:
        if evtmatch.startswith( '-e' ) and evtmatch[2:].isdigit(): ftitle += ' for %d cycles'%int(evtmatch[2:])
        fig.suptitle( ftitle, size=plots_ftitlesize )
    fig.set_tight_layout( True )
    fig.savefig( pngpath, format='png', bbox_inches="tight" )
    Popen( ['display', '-geometry', '+50+50', '-resize', '300', pngpath] )
    print( 'Plot successfully saved on', pngpath )

#---------------------------------------

def allCppPlots( workdir, evtmatch='-e001', debug=False ):
    plotCppST( workdir, keymatch='sa-cpp-d-inl0-best', ylog=True, evtmatch=evtmatch )
    plotCppST( workdir, keymatch='sa-cpp-f-inl0-best', ylog=True, evtmatch=evtmatch )
    plotCppOneProcess2( workdir, 'ggttgg', 'sa-cpp-d-inl0', evtmatch=evtmatch )
    plotCppOneProcess2( workdir, 'ggttgg', 'sa-cpp-f-inl0', evtmatch=evtmatch )
    if 'pmpe04' in workdir or 'itscrd70' in workdir or ( 'bmk6130' in workdir and '-e010' in evtmatch ) :
        plotCppOneProcess2( workdir, None, 'sa-cpp-d-inl0', evtmatch=evtmatch )
        plotCppOneProcess2( workdir, None, 'sa-cpp-f-inl0', evtmatch=evtmatch )
        plotCppProcessesInl( workdir, 'sa-cpp-d-inl', evtmatch=evtmatch )
        plotCppProcessesInl( workdir, 'sa-cpp-f-inl', evtmatch=evtmatch )

#---------------------------------------

# Compare various curves in cuda ST plots
def axesCudaST( ax, cudarunset_scores, score_key, xjob=False, xlog=True, abstput=True, ylog=False, debug=False ):
    ngbls = set( [njobnthr[0] for njobnthr in cudarunset_scores] ) # use set(list) to get unique keys
    ngths = set( [njobnthr[1] for njobnthr in cudarunset_scores] ) # use set(list) to get unique keys
    njobs = set( [njobnthr[2] for njobnthr in cudarunset_scores] ) # use set(list) to get unique keys
    nthrs = set( [njobnthr[3] for njobnthr in cudarunset_scores] ) # use set(list) to get unique keys
    assert len(ngths) == 1, 'assume all jobs have the same number of GPU threads'
    ngth0 = list(ngths)[0]
    # Prepare axes labels
    if xjob: ax.set_xlabel('nblocksGPU * nthreadsGPU * njobsCPU', size=plots_labelsize )
    else: ax.set_xlabel('nblocksGPU * nthreadsGPU', size=plots_labelsize )
    if abstput:
        ax.set_ylabel('Throughput (E6 events per second)', size=plots_labelsize )
    else:
        ax.set_ylabel('Throughput ratio to (1,%d,1)'%ngth0, size=plots_labelsize )
        ax.grid()
        assert (1,ngth0,1,1) in cudarunset_scores, 'no scores found for (ngbl,ngth,njob,nthr)==(1,%d,1,1)?'%ngth0
        tputmin = cudarunset_scores[1,ngth0,1,1][score_key]
    if xlog: ax.set_xscale( 'log' )
    if ylog: ax.set_yscale( 'log' )
    # Add one curve per njob
    nthr = 1
    assert nthrs == {nthr}, 'assume ST CPU jobs with nthr==1'
    xmax = 0
    ymax = 0
    for njob in sorted(njobs) :
        # Prepare x-axis and y-axis lists
        xvals = []
        yvals = []
        for ngbl in sorted(ngbls):
            for ngth in sorted(ngths):
                if (ngbl,ngth,njob,nthr) not in cudarunset_scores: continue
                if xjob: xval = ngbl*ngth*njob # 'gpugridsize*njob'
                else: xval = ngbl*ngth # 'gpugridsize'
                tput = cudarunset_scores[ngbl,ngth,njob,nthr][score_key]
                xvals.append( xval )
                if abstput: yvals.append( tput )
                else: yvals.append( tput / tputmin )
        if xjob: xmax = max( xmax, max(xvals) )
        else: xmax = max( xmax, max(xvals)*njob ) # use the same xmax both with/without xjob
        ymax = max( ymax, max(yvals) )
        # Add curve of y vs x
        p = ax.plot( xvals, yvals, marker='o', label=score_key+' (njobsCPU=%d)'%njob )
    # Decorate axes
    loc = 'lower right'
    ax.legend( loc=loc, fontsize=plots_legendsize )
    if xlog: xmin, xmax = 100, xmax*4
    else: xmin, xmax = 0, xmax*1.1
    if ylog: ymin, ymax = 0.001, ymax*12
    else: ymin, ymax = 0, ymax*1.1
    ax.axis( [xmin, xmax, ymin, ymax] )

# Create a figure with 2x2 plots (use gridsize or gridsize*njob as x-axis; use absolute or relative tput as y-axis)
def plotCudaST( workdir, score_key='ggttgg-sa-cuda-d-inl0', ylog=False, evtmatch='-e001', gputhreads='gt00256', debug=False ):
    cudarunset_scores = loadCudaRunSet( workdir, evtmatch=evtmatch, gputhreads=gputhreads )
    node, xht, ftitle = getNodeFeatures( workdir )
    pngpath = workdir + '/' + node + '-' + gputhreads + evtmatch + '-' + score_key + '.png'
    # Create figure with 2x2 plots
    fig = plt.figure( figsize = ( 1.2*plots_figsize[0]*2, 1.2*plots_figsize[1]*2 ) )
    ax1 = fig.add_subplot( 221 )
    ax2 = fig.add_subplot( 222 )
    ax3 = fig.add_subplot( 223 )
    ax4 = fig.add_subplot( 224 )
    # Fill the plots in the figure
    axesCudaST( ax1, cudarunset_scores, xjob=False, score_key=score_key, abstput=True, ylog=ylog, debug=debug )
    axesCudaST( ax2, cudarunset_scores, xjob=True, score_key=score_key, abstput=True, ylog=ylog, debug=debug )
    axesCudaST( ax3, cudarunset_scores, xjob=False, score_key=score_key, abstput=False, ylog=ylog, debug=debug )
    axesCudaST( ax4, cudarunset_scores, xjob=True, score_key=score_key, abstput=False, ylog=ylog, debug=debug )
    ax1.add_artist( matplotlib.offsetbox.AnchoredText( gputhreads, loc='upper left' ) )
    ax2.add_artist( matplotlib.offsetbox.AnchoredText( gputhreads, loc='upper left' ) )
    ax3.add_artist( matplotlib.offsetbox.AnchoredText( gputhreads, loc='upper left' ) )
    ax4.add_artist( matplotlib.offsetbox.AnchoredText( gputhreads, loc='upper left' ) )
    if ftitle is not None:
        if evtmatch.startswith( '-e' ) and evtmatch[2:].isdigit(): ftitle += ' for %d cycles'%int(evtmatch[2:])
        ###ftitle += ' - ' + ( 'DOUBLE' if '-d-' in score_key else 'FLOAT' )
        ###ftitle += ' - ' + score_key
        fig.suptitle( ftitle, size=plots_ftitlesize )
    # Save and show the figure
    fig.savefig( pngpath, format='png', bbox_inches="tight" )
    Popen( ['display', '-geometry', '+50+50', pngpath] )
    ###Popen( ['display', '-geometry', '+50+50', '-resize', '800', pngpath] )
    print( 'Plot successfully saved on', pngpath )

#---------------------------------------

def compareNodesCpp( ax=None, process=None, df='-d-', inl='inl0' ):
    plot = ax is not None
    nodes = ( 'itscrd70', 'pmpe04', 'bmk6130' )
    refnode = nodes[0] # refnode is the first one (itscrd70)
    npcores_node = {}
    tputs_node = {}
    keys = []
    for node in nodes:
        workdir = 'BMK-' + node
        cpprunset_scores = loadCppRunSet( workdir, evtmatch='-e010' )
        node1, xht, ftitle = getNodeFeatures( workdir )
        assert node == node1, 'node mismatch in getNodeFeatures?'
        npcores = abs( xht ) # use benchmarks for the case when njob equals the number of physical cores
        assert (npcores,1) in cpprunset_scores, 'no scores found for njob==%d and nthr==1?'%npcores
        tputs = cpprunset_scores[ npcores, 1 ]
        npcores_node[ node ] = npcores
        tputs_node[ node ] = tputs
        keys += tputs.keys()
    keys = set( keys )
    def sortableSimdKey( key ): # use keys sortable in this order: none, sse4, avx2, 512y, 512z, best
        if   key.endswith( '-none' ) : key = 'simd0' + key
        elif key.endswith( '-sse4' ) : key = 'simd1' + key
        elif key.endswith( '-avx2' ) : key = 'simd2' + key
        elif key.endswith( '-512y' ) : key = 'simd3' + key
        elif key.endswith( '-512z' ) : key = 'simd4' + key
        elif key.endswith( '-best' ) : key = 'simd5' + key
        return key
    keys2 = [ sortableSimdKey( key ) for key in keys ]
    keys = [ key for _, key in sorted( zip( keys2, keys ) ) ] # https://stackoverflow.com/a/6618543
    if not plot:
        print( '%-30s %-10s %10s %10s %12s %15s'%( 'ScoreName', 'Node', 'Tput', 'Tput/Ref', 'TputPerCore', 'TputPerCore/Ref' ) )
    else:
        xkey = 'ggttgg-sa-cpp-d-inl0-none' # Use ggttgg/d/none as the x axis (the hepspec06-like benchmark)
        xmin, ymin = 1, 1
        ###xmax, ymax = 1, 1
        ax.set_xlabel( 'Tput(node)/Tput(ref_node) for %s'%xkey )
        ax.set_ylabel( 'Tput(node)/Tput(ref_node)' )
    for key in keys :
        if process is not None and not key.startswith( process ): continue # look only at one specific process
        if not df in key: continue # look at -d- or -f- only
        if not '-%s-'%inl in key: continue # look at inl0 or inl1 only
        if plot: xvals, yvals = [], [] # prepare data to plot
        else: print()
        for node in nodes:
            if key not in tputs_node[node] : continue
            tput = tputs_node[node][key]
            reftput = tputs_node[refnode][key]
            tputpc = tput / npcores_node[node]
            reftputpc = tputs_node[refnode][key] / npcores_node[refnode]
            if not plot:
                print( '%-30s %-10s %10.5f %10.5f %12.5f %15.5f'%( key, node, tput, tput / reftput, tputpc, tputpc / reftputpc ) )
            else:
                xtput = tputs_node[node][xkey]
                xreftput = tputs_node[refnode][xkey]
                xvals.append( xtput / xreftput )
                yvals.append( tput / reftput )
        if plot:
            ###xmax = max( xmax, max(xvals) )
            ###ymax = max( ymax, max(yvals) )
            ax.plot( xvals, yvals, marker='o', label=key )
    if plot:
        ###xmax = max( xmax*1.4, ymax*1.2 )
        xmax = 14 # hardcoded default (same for all plots)
        ymax = xmax
        ax.axis( [xmin, xmax, ymin, ymax] )
        loc = 'lower right'
        ax.legend( loc=loc, fontsize=plots_legendsize )
        ax.plot( [xmin, xmax], [ymin, ymax], ls="--", c='black' )

def compareNodesCppPlot( df='-d-', inl='inl0' ):
    fig = plt.figure( figsize = ( plots_figsize[0]*1.4, plots_figsize[1]*2.8 ) )
    # Create figure with 2x2 plots
    ax1 = fig.add_subplot( 221 )
    ax2 = fig.add_subplot( 222 )
    ax3 = fig.add_subplot( 223 )
    ax4 = fig.add_subplot( 224 )
    compareNodesCpp( ax1, df=df, inl=inl, process='eemumu-' )
    compareNodesCpp( ax2, df=df, inl=inl, process='ggtt-' )
    compareNodesCpp( ax3, df=df, inl=inl, process='ggttg-' )
    compareNodesCpp( ax4, df=df, inl=inl, process='ggttgg-' )
    pngpath = 'BMK-COMPARE/all-sa-cpp%s%s.png'%(df,inl)
    fig.set_tight_layout( True )
    fig.savefig( pngpath, format='png', bbox_inches="tight" )
    Popen( ['display', '-geometry', '+50+50', pngpath] )
    print( 'Plot successfully saved on', pngpath )

#---------------------------------------

if __name__ == '__main__':

    # TESTS (CPP)
    #loadOneRun( 'BMK-pmpe04/sa-cpp-j032-t001-e001', debug=True )
    #loadCppRunSet( 'BMK-pmpe04', debug=True )
    #dumpCppScoresOneKey( loadCppRunSet( 'BMK-pmpe04' ), 'ggttgg-sa-cpp-d-inl0-best' )
    #dumpCppScoresAllKeys( loadCppRunSet( 'BMK-pmpe04' ) )
    #dumpCppScoresAllKeys( loadCppRunSet( 'BMK-pmpe04'), keymatch='best' )
    #dumpCppScoresAllKeys( loadCppRunSet( 'BMK-pmpe04'), keymatch='inl0-best' )
    #dumpCppScoresAllKeys( loadCppRunSet( 'BMK-pmpe04'), keymatch='ggttgg-sa-cpp-d-inl0' )

    # PRODUCTION PLOTS (CPP)
    #allCppPlots( 'BMK-pmpe04', '-e001' )
    #allCppPlots( 'BMK-pmpe04', '-e010' )
    #allCppPlots( 'BMK-itscrd70', '-e001' )
    #allCppPlots( 'BMK-itscrd70', '-e010' )
    #allCppPlots( 'BMK-jwlogin08', '-e001' )
    #allCppPlots( 'BMK-bmk6130', '-e001' )
    #allCppPlots( 'BMK-bmk6130', '-e010' )

    # TESTS (CUDA)
    #loadCudaRunSet( 'BMK-itscrd70-cuda', evtmatch='-e0100', debug=True )
    #dumpCudaScoresOneKey( loadCudaRunSet( 'BMK-itscrd70-cuda', evtmatch='-e0100' ), 'ggttgg-sa-cuda-d-inl0' )
    #dumpCudaScoresOneKey( loadCudaRunSet( 'BMK-itscrd70-cuda', evtmatch='-e0100' ), 'ggttgg-sa-cuda-f-inl0' )
    #dumpCudaScoresOneKey( loadCudaRunSet( 'BMK-itscrd70-cuda', evtmatch='-e0800' ), 'ggttgg-sa-cuda-d-inl0' )
    #dumpCudaScoresOneKey( loadCudaRunSet( 'BMK-itscrd70-cuda', evtmatch='-e0800' ), 'ggttgg-sa-cuda-f-inl0' )
    #dumpCudaScoresOneKey( loadCudaRunSet( 'BMK-itscrd70-cuda', evtmatch='-e0800', gputhreads='gt00032' ), 'ggttgg-sa-cuda-d-inl0' )
    #dumpCudaScoresOneKey( loadCudaRunSet( 'BMK-itscrd70-cuda', evtmatch='-e0800', gputhreads='gt00032' ), 'ggttgg-sa-cuda-f-inl0' )

    # PRODUCTION PLOTS (CUDA)
    #plotCudaST( 'BMK-itscrd70-cuda', score_key='ggttgg-sa-cuda-d-inl0', ylog=False, evtmatch='-e0100', debug=False )
    #plotCudaST( 'BMK-itscrd70-cuda', score_key='ggttgg-sa-cuda-f-inl0', ylog=False, evtmatch='-e0100', debug=False )
    #plotCudaST( 'BMK-itscrd70-cuda', score_key='ggttgg-sa-cuda-d-inl0', ylog=False, evtmatch='-e0800', debug=False )
    #plotCudaST( 'BMK-itscrd70-cuda', score_key='ggttgg-sa-cuda-f-inl0', ylog=False, evtmatch='-e0800', debug=False )
    #plotCudaST( 'BMK-itscrd70-cuda', score_key='ggttgg-sa-cuda-d-inl0', ylog=False, evtmatch='-e0800', gputhreads='gt00032', debug=False )
    #plotCudaST( 'BMK-itscrd70-cuda', score_key='ggttgg-sa-cuda-f-inl0', ylog=False, evtmatch='-e0800', gputhreads='gt00032', debug=False )

    # TEST COMPARISONS (CPP)
    #compareNodesCpp( process='ggtt-' )

    # PRODUCTION COMPARISONS (CPP)
    compareNodesCpp()
    #compareNodesCppPlot()
