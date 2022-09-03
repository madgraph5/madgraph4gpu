#!/usr/bin/env python3

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

def loadRunSet( runsetdir, debug=False ):
    ###debug=True
    import os
    import json
    if not os.path.isdir( runsetdir ):
        print( 'Unknown directory', runsetdir )
        return
    # Go through all runs in runsetdir and fill runset_scores[njob,nthr]
    runset_scores = {}
    print( 'Loading runs in RunSetDir', runsetdir )
    for d in sorted( os.listdir( runsetdir ) ) :
        if d.find( 'sa-cpp-' ) != -1: # e.g. sa-cpp-j004-t001-e001
            dl = d.split( '-' )
            njob = int( dl[-3][-3:] )
            nthr = int( dl[-2][-3:] )
            nevt = int( dl[-1][-3:] )
            rundir = runsetdir + '/' + d
            if debug : print( '\nRunDir=%30s %3i %3i %3i'%( rundir, njob, nthr, nevt ) )
            run_info, run_scores = loadOneRun( rundir ) 
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

def dumpScoresAllKeys( runset_scores, keymatch=None, debug=False ):
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
    for key in keys : dumpScoresOneKey( runset_scores, key )

#---------------------------------------

if __name__ == '__main__':

    #loadOneRun( 'BMK-pmpe04/sa-cpp-j032-t001-e001', debug=True )
    #loadRunSet( 'BMK-pmpe04', debug=True )
    #dumpScoresOneKey( loadRunSet( 'BMK-pmpe04' ), 'ggttgg-sa-cpp-d-inl0-best' )
    #dumpScoresAllKeys( loadRunSet( 'BMK-pmpe04' ) )
    #dumpScoresAllKeys( loadRunSet( 'BMK-pmpe04'), keymatch='best' )
    dumpScoresAllKeys( loadRunSet( 'BMK-pmpe04'), keymatch='inl0-best' )
    #dumpScoresAllKeys( loadRunSet( 'BMK-pmpe04'), keymatch='ggttgg-sa-cpp-d-inl0' )

