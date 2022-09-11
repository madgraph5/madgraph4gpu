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

def loadRunSet( runsetdir, evtmatch='-e001', debug=False ):
    ###debug=True
    if not os.path.isdir( runsetdir ):
        print( 'Unknown directory', runsetdir )
        return
    # Go through all runs in runsetdir and fill runset_scores[njob,nthr]
    runset_scores = {}
    print( 'Loading runs in RunSetDir', runsetdir, 'for events', evtmatch )
    for d in sorted( os.listdir( runsetdir ) ) :
        if d.startswith( 'sa-cpp-' ) and d.endswith( evtmatch ) and 'png' not in d : # e.g. sa-cpp-j004-t001-e001
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
            runset_scores[njob,nthr] = run_scores
    return runset_scores

#---------------------------------------

def dumpScoresOneKey( runset_scores, score_key, debug=False ):
    ###debug=True
    print( '\nSCORES[\'%s\']:'%score_key )
    score_key_none = score_key[:-4]+'none'
    njobs = set( [njobnthr[0] for njobnthr in runset_scores] ) # use set(list) to get unique keys
    nthrs = set( [njobnthr[1] for njobnthr in runset_scores] ) # use set(list) to get unique keys
    print( '%4s %4s %12s    %9s %12s %12s %16s'%( 'njob', 'nthr', 'Score', 'njob*nthr', 'S/S[1,1]', 'S/S-none', 'S/S-none[1,1]' ) )
    assert (1,1) in runset_scores, 'no scores found for njob==1 and nthr==1?'
    tput1 = runset_scores[1,1][score_key]
    tput1none = runset_scores[1,1][score_key_none]
    for nthr in sorted(nthrs):
        for njob in sorted(njobs):
            if (njob,nthr) not in runset_scores: continue
            tput = runset_scores[njob,nthr][score_key]
            tputnone = runset_scores[njob,nthr][score_key_none]
            print( '%4d %4d %12.6f    %9d %12.6f %12.6f %16.6f'%
                   ( njob, nthr, tput, njob*nthr, tput / tput1, tput / tputnone, tput / tput1none ) )

#---------------------------------------

def getSortedMatchingKeys( runset_scores, keymatch=None, debug=False ):
    ###debug=True
    keys = []
    for njobnthr in runset_scores : keys += list( runset_scores[njobnthr].keys() )
    keys = set( keys ) # use set(list) to get unique keys
    if keymatch is not None: keys = [ key for key in keys if keymatch in key ]
    def sortableSimdKey( key ): # use keys sortable in this order: none, sse4, avx2, 512y, 512z, best
        key = key.replace( '-none', '-simd0' )
        key = key.replace( '-sse4', '-simd1' )
        key = key.replace( '-avx2', '-simd3' )
        key = key.replace( '-512y', '-simd4' )
        key = key.replace( '-512z', '-simd5' )
        key = key.replace( '-best', '-simd6' )
        return key
    keys2 = [ sortableSimdKey( key ) for key in keys ]
    keys = [ key for _, key in sorted( zip( keys2, keys ) ) ] # https://stackoverflow.com/a/6618543
    return keys

#---------------------------------------

def dumpScoresAllKeys( runset_scores, keymatch=None, debug=False ):
    keys = getSortedMatchingKeys( runset_scores, keymatch, debug )
    for key in keys : dumpScoresOneKey( runset_scores, key )

#---------------------------------------

# Compare various curves in ST plots
# NB: xht>0 is the number of physical cores before 2xHT, xht<0 is the number of physical cores without HT
def axesST( ax, runset_scores, keymatch=None, bestonly=False, abstput=True, ylog=False, xht=None, debug=False ):
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
    keys = getSortedMatchingKeys( runset_scores, keymatch, debug )
    if bestonly : keys = [ key for key in keys if 'best' in key ]
    elif 'best' not in keymatch: keys = [ key for key in keys if 'best' not in key ]
    for score_key in keys :
        score_key_none = score_key[:-4]+'none'
        njobs = set( [njobnthr[0] for njobnthr in runset_scores] ) # use set(list) to get unique keys
        nthrs = set( [njobnthr[1] for njobnthr in runset_scores] ) # use set(list) to get unique keys
        assert (1,1) in runset_scores, 'no scores found for njob==1 and nthr==1?'
        ###tput1 = runset_scores[1,1][score_key]
        tput1none = runset_scores[1,1][score_key_none]
        # Prepare x-axis and y-axis lists
        xvals = []
        yvals = []
        for nthr in sorted(nthrs):
            for njob in sorted(njobs):
                if (njob,nthr) not in runset_scores: continue
                xval = nthr*njob # 'npar' level of parallelism
                tput = runset_scores[njob,nthr][score_key]
                ###tputnone = runset_scores[njob,nthr][score_key_none]
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

# Create a figure with a single plot
def plotST( workdir, keymatch=None, abstput=True, ylog=False, evtmatch='-e001', debug=False ):
    runset_scores = loadRunSet( workdir, evtmatch=evtmatch )
    node, xht, ftitle = getNodeFeatures( workdir )
    pngpath = workdir + '/' + node + evtmatch + '-all-' + keymatch + '.png'
    # Create figure with one plot
    fig = plt.figure( figsize=plots_figsize )
    ax1 = fig.add_subplot( 111 )
    # Fill the plot in the figure
    axesST( ax1, runset_scores, keymatch=keymatch, ylog=ylog, xht=xht, abstput=abstput, debug=debug )
    if ftitle is not None:
        if evtmatch.startswith( '-e' ) and evtmatch[2:].isdigit(): ftitle += ' for %d cycles'%int(evtmatch[2:])
        fig.suptitle( ftitle, size=plots_ftitlesize )
    # Save and show the figure
    fig.savefig( pngpath, format='png', bbox_inches="tight" )
    Popen( ['display', '-geometry', '+50+50', pngpath] )
    ###Popen( ['display', '-geometry', '+50+50', '-resize', '800', pngpath] )
    print( 'Plot successfully saved on', pngpath )

# Create a figure with two plots per process, absolute and normalized tput
def plotOneProcess2( workdir, oneprocess, keymatch, bestonly=False, evtmatch='-e001', debug=False ):
    runset_scores = loadRunSet( workdir, evtmatch=evtmatch )
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
        axesST( ax1, runset_scores, keymatch=fullkeymatch, bestonly=bestonly, \
                ylog=False, xht=xht, abstput=True, debug=debug )
        ax2 = fig.add_subplot( idx2 )
        axesST( ax2, runset_scores, keymatch=fullkeymatch, bestonly=bestonly, \
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
def plotProcessesInl( workdir, keymatch, evtmatch='-e001', debug=False ):
    runset_scores = loadRunSet( workdir, evtmatch=evtmatch )
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
        axesST( ax1, runset_scores, keymatch=fullkeymatch, bestonly=True, \
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

if __name__ == '__main__':

    # TESTS
    #loadOneRun( 'BMK-pmpe04/sa-cpp-j032-t001-e001', debug=True )
    #loadRunSet( 'BMK-pmpe04', debug=True )
    #dumpScoresOneKey( loadRunSet( 'BMK-pmpe04' ), 'ggttgg-sa-cpp-d-inl0-best' )
    #dumpScoresAllKeys( loadRunSet( 'BMK-pmpe04' ) )
    #dumpScoresAllKeys( loadRunSet( 'BMK-pmpe04'), keymatch='best' )
    #dumpScoresAllKeys( loadRunSet( 'BMK-pmpe04'), keymatch='inl0-best' )
    #dumpScoresAllKeys( loadRunSet( 'BMK-pmpe04'), keymatch='ggttgg-sa-cpp-d-inl0' )

    # PRODUCTION PLOTS
    #workdir = 'BMK-pmpe04'
    #workdir = 'BMK-itscrd70'
    #workdir = 'BMK-jwlogin08'
    workdir = 'BMK-bmk6130'
    evtmatch='-e001'
    #evtmatch='-e010'
    plotST( workdir, keymatch='sa-cpp-d-inl0-best', ylog=True, evtmatch=evtmatch )
    plotST( workdir, keymatch='sa-cpp-f-inl0-best', ylog=True, evtmatch=evtmatch )
    plotOneProcess2( workdir, 'ggttgg', 'sa-cpp-d-inl0', evtmatch=evtmatch )
    plotOneProcess2( workdir, 'ggttgg', 'sa-cpp-f-inl0', evtmatch=evtmatch )
    #plotOneProcess2( workdir, None, 'sa-cpp-d-inl0', evtmatch=evtmatch )
    #plotOneProcess2( workdir, None, 'sa-cpp-f-inl0', evtmatch=evtmatch )
    #plotProcessesInl( workdir, 'sa-cpp-d-inl', evtmatch=evtmatch )
    #plotProcessesInl( workdir, 'sa-cpp-f-inl', evtmatch=evtmatch )

