#!/usr/bin/env python3

import os

#------------------------------------------------------------------------------

# Format tuple [RSS, PSS, SWAP] where measurements can be None
def memstr(mem):
    ###print( mem )
    out= '[%8s, %8s, %8s]'%\
         tuple('%8.1f'%m if m is not None else 'None' for m in mem)
    ###print( out )
    return out

# Parse mem.xml for a job to extract max RSS/PSS/SWAP over the job time range
def parseMemXml(fileName, debug=False):
    f = open(fileName, 'r')
    memsum = [None, None, None] # rss, pss, swp (sum over all processes of a MP job)
    maxmemsum = [None, None, None] # rss, pss, swp (max sum over all processes of a MP job)
    for l in f.readlines():
        l = l.replace('<',' ').replace('>',' ')
        # New memprof line group (new time, new memory measurements)
        if l.startswith(' memprof'):
            if debug: print( memstr(memsum), '\n\n', l.split()[1] )
            for i in range(3):
                if maxmemsum[i] is None or memsum[i] > maxmemsum[i]:
                    maxmemsum[i] = memsum[i]
            memsum = [None, None, None] # rss, pss, swp (sum over all processes)
        # New process or thread line
        if not l.startswith(' process') and not l.startswith(' thread') :
            continue
        lsplit = l.split()
        ###print( lsplit )
        pid = lsplit[1]
        if lsplit[2].startswith('ppid='):
            if not lsplit[4].startswith('rss=') :
                raise Exception('Internal error - no rss in %s'%lsplit[4])
            rss = int(lsplit[4].split('"')[1])/1000.0 # MB (kB / 1000)
            if not lsplit[5].startswith('pss=') :
                raise Exception('Internal error - no pss in %s'%lsplit[5])
            pss = int(lsplit[5].split('"')[1])/1000.0 # MB (kB / 1000)
            if lsplit[6].startswith('swap='):
                # lbsmaps version "21": BC's second version, plus AV's swap
                swp = int(lsplit[6].split('"')[1])/1000.0 # MB (kB / 1000)
                cmd = lsplit[7]
            else:
                # lbsmaps version "20": BC's second version, without AV's swap
                swp = None
                cmd = lsplit[6]
        else:
            if not lsplit[3].startswith('rss=') :
                raise Exception('Internal error - no rss in %s'%lsplit[3])
            rss = int(lsplit[3].split('"')[1])/1000.0 # MB (kB / 1000)
            if not lsplit[4].startswith('pss=') :
                raise Exception('Internal error - no pss in %s'%lsplit[4])
            pss = int(lsplit[4].split('"')[1])/1000.0 # MB (kB / 1000)
            if lsplit[5].startswith('swap='):
                # lbsmaps version "11": BC's first version, plus AV's swap
                swp = int(lsplit[5].split('"')[1])/1000.0 # MB (kB / 1000)
                cmd = lsplit[6]
            else:
                # lbsmaps version "10": BC's first version, without AV's swap
                swp = None
                cmd = lsplit[5]
        mem = [rss, pss, swp]
        ###if debug: print( memstr(memsum), memstr(mem), cmd )
        if not cmd.startswith('{') and not cmd.endswith('}'):
            if debug: print( memstr(memsum), memstr(mem), cmd ) # NB check.exe has a single process!
            for i in range(3):
                if memsum[i] is None : memsum[i] = mem[i]
                else : memsum[i] += mem[i]
    f.close()
    if debug: print( '\nmaxmemsum\n', memstr(maxmemsum) )
    for i in range(3):
        if maxmemsum[i] is not None: maxmemsum[i] = int(maxmemsum[i]*1000)/1000.
    ###print( 'maxmemsum', maxmemsum )
    return maxmemsum

# Parse log.txt for a job to extract total times
def parseLogTxt(file):
    ###print( 'FILE:', file )
    tot, tot1, tot2, tot3, nevt = 0, 0, 0, 0, 0
    with open(file, 'r') as fh:
        for line in fh.readlines() :
            lline = line.split()
            if lline[0]=='TotalEventsComputed' : nevt=float(lline[2])
            elif lline[0]=='TOTAL' :
                if lline[1]==':' : tot=float(lline[2])
                elif lline[1]=='(1)' : tot1=float(lline[3])
                elif lline[1]=='(2)' : tot2=float(lline[3])
                elif lline[1]=='(3)' : tot3=float(lline[3])
    ###print( tot, tot1, tot2, tot3, nevt )
    return tot, tot1, tot2, tot3, nevt

# Append to timdata the job times extracted from log.txt
# Here timdata[ijob] = [tot, tot1, tot2, tot3, nevt]
def appendTimdata(timdata, ijob, file, debug=False):
    tot, tot1, tot2, tot3, nevt = parseLogTxt(file)
    if debug: print( "appendTimdata", ijob, tot, tot1, tot2, tot3, nevt )
    timdata[ijob] = ( tot, tot1, tot2, tot3, nevt )

# Append to memdata the rss, pss extracted from mem.txt
# Here memdata[ijob] = [max sumrss, max sumpss, max sumswap]
# where max is over all measurements in the job time range
# and sum is the sum over all processes existing at a given time
def appendMemdata(memdata, ijob, file, debug=False):
    rss, pss, swp = parseMemXml(file)
    if debug: print( "appendMemdata", ijob, rss, pss, swp )
    memdata[ijob] = ( rss, pss, swp )

# Browse the test directory for a given njob and nthr and parse files therein
def extractTimdataAndMemdata(workdir, debug=False):
    print( '-> Processing files in', workdir )
    # Extract timdata and memdata
    timdata = {} # x[ijob] = [tot, tot1, tot2, tot3, nevt]
    memdata = {} # x[ijob] = [max sumrss, max sumpss, max sumswap]
    for x in sorted(os.listdir(workdir)):
        workdir2 = workdir + x
        if os.path.isdir(workdir2) and x.isdigit():
            ijob = int(x)
            if debug : print( '---> Processing files in', workdir2 )
            for y in sorted(os.listdir(workdir2)):
                if y == 'log.txt':
                    logfile = workdir2 + os.sep + y
                    if os.path.isfile(logfile):
                        if debug : print( '----> Parsing', logfile )
                        appendTimdata(timdata, ijob, logfile, debug)
                elif y == 'mem.txt':
                    memfile = workdir2 + os.sep + y
                    if os.path.isfile(memfile):
                        if debug : print( '----> Parsing', memfile )
                        appendMemdata(memdata, ijob, memfile, debug)
    return timdata, memdata

#------------------------------------------------------------------------------

timsum_njob_nthr = {} # x[njob,nthr] = [avgtot, avgtot1, avgtot2, avgtot3, sumnevt] -> avg times and sum nevt over all jobs
memsum_njob_nthr = {} # x[njob,nthr] = [sumrss, sumpss, sumswap] -> sum over all jobs (and threads/processes of MT/MP jobs)

allnjob = []
allnthr = []

def processTimdata(njob, nthr, timdata, debug=False):
    ###debug=True
    ###if debug : print( 'processTimdata', njob, nthr, timdata )
    ###print( 'Running1', timsum_njob_nthr )
    sumtot = 0
    sumtot1 = 0
    sumtot2 = 0
    sumtot3 = 0
    sumnevt = 0
    for x in sorted(timdata): # one per job
        sumtot += timdata[x][0]
        sumtot1 += timdata[x][1]
        sumtot2 += timdata[x][2]
        sumtot3 += timdata[x][3]
        sumnevt += timdata[x][4]
    avgtot = sumtot/njob
    avgtot1 = sumtot1/njob
    avgtot2 = sumtot2/njob
    avgtot3 = sumtot3/njob
    if njob not in timsum_njob_nthr: timsum_njob_nthr[njob] = {}
    timsum_njob_nthr[njob][nthr] = [avgtot, avgtot1, avgtot2, avgtot3, sumnevt]
    ###print( 'Running2', timsum_njob_nthr )
    
def processMemdata(njob, nthr, memdata, debug=False):
    ###debug=True
    ###if debug : print( 'processMemdata', njob, nthr, memdata )
    ###print( 'Running1', memsum_njob_nthr )
    sumrss = 0
    sumpss = 0
    sumswp = None
    for x in sorted(memdata): # one per job
        sumrss += memdata[x][0]
        sumpss += memdata[x][1]
        if sumswp is None : sumswp = memdata[x][2]
        else : sumswp += memdata[x][2]
    if njob not in memsum_njob_nthr: memsum_njob_nthr[njob] = {}
    memsum_njob_nthr[njob][nthr] = [sumrss, sumpss, sumswp]
    ###print( 'Running2', memsum_njob_nthr )
    
def processFiles(rundir, avx='none', debug=False):
    ###debug=True
    if not os.path.isdir(rundir):
        print( 'Unknown directory', rundir )
        return
    # GO THROUGH ALL TESTS IN RUNDIR AND FILL TIMSUM_NJOB_NTHR AND MEMSUM_NJOB_NTHR
    global timsum_njob_nthr
    global memsum_njob_nthr
    timsum_njob_nthr = {}
    memsum_njob_nthr = {}
    print( 'Processing files in', rundir )
    curlist = os.listdir(rundir)
    curlist.sort()
    for d in curlist:
        if d.find('check-test.'+avx+'.') != -1: # e.g. check-test.none.j002.t004
            dl = d.split('.')
            njob = int(dl[-2][-3:])
            nthr = int(dl[-1][-3:])
            workdir = rundir + '/' + d + '/'
            print( 'Workdir', workdir, njob, nthr )
            data = extractTimdataAndMemdata(workdir, debug)
            if data is not None:
                timdata = data[0]
                memdata = data[1]
                processTimdata(njob, nthr, timdata, debug)
                processMemdata(njob, nthr, memdata, debug)
            else:
                print( '   Empty data found in', workdir )
    # DEBUG: DUMP TIMSUM_NJOB_NTHR AND MEMSUM_NJOB_NTHR
    print( 'Processed files' )
    ###print( 'timsum_njob_nthr', timsum_njob_nthr )
    #print( '%4s %4s %11s    %9s %11s'%( 'njob', 'nthr', 'tput[Mev/s]', 'njob*nthr', 'tput/tput1' ) )
    #tput1=0
    #for njob in sorted(timsum_njob_nthr) :
    #    for nthr in sorted(timsum_njob_nthr[njob]) :
    #        tput = timsum_njob_nthr[njob][nthr][4] / timsum_njob_nthr[njob][nthr][3] / 1E6
    #        if tput1 == 0 : tput1=tput
    #        print( '%4d %4d %11.5f    %9d %11.5f'%( njob, nthr, tput, njob*nthr, tput / tput1 ) )
    ###print( 'memsum_njob_nthr', memsum_njob_nthr )
    #print( '%4s %4s %11s'%( 'njob', 'nthr', 'maxPSS[GB]' ) )
    #for njob in sorted(memsum_njob_nthr) :
    #    for nthr in sorted(memsum_njob_nthr[njob]) :
    #        maxpss = memsum_njob_nthr[njob][nthr][1]/1000
    #        print( '%4d %4d %11.5f'%( njob, nthr, maxpss ) )
    # DUMP THEM TOGETHER (ASSUMING THE SAME ARRAY STRUCTURE)
    #print( '%4s %4s %11s    %9s %11s   %11s'%( 'njob', 'nthr', 'tput[Mev/s]', 'njob*nthr', 'tput/tput1', 'maxPSS[GB]' ) )
    #tput1=0
    #for njob in sorted(timsum_njob_nthr) :
    #    for nthr in sorted(timsum_njob_nthr[njob]) :
    #        tput = timsum_njob_nthr[njob][nthr][4] / timsum_njob_nthr[njob][nthr][3]
    #        if tput1 == 0 : tput1=tput
    #        maxpss = memsum_njob_nthr[njob][nthr][1]/1000
    #        print( '%4d %4d %11.5f    %9d %11.5f   %11.5f'%( njob, nthr, tput / 1E6, njob*nthr, tput / tput1, maxpss ) )
    # DUMP THEM TOGETHER (ASSUMING THE SAME ARRAY STRUCTURE)
    global allnjob
    global allnthr
    allnjob = sorted(timsum_njob_nthr)
    allnthr = []
    for njob in allnjob: allnthr += sorted(timsum_njob_nthr[njob])
    allnthr = list(set(allnthr)) # unique items
    print( '%4s %4s %11s    %9s %11s   %11s'%( 'njob', 'nthr', 'tput[Mev/s]', 'njob*nthr', 'tput/tput1', 'maxPSS[GB]' ) )
    tput1=0
    for nthr in sorted(allnthr):
        for njob in sorted(allnjob):
            if njob not in timsum_njob_nthr: continue
            if nthr not in timsum_njob_nthr[njob]: continue
            tput = timsum_njob_nthr[njob][nthr][4] / timsum_njob_nthr[njob][nthr][3]
            if tput1 == 0 : tput1=tput
            if njob in memsum_njob_nthr and nthr in memsum_njob_nthr[njob]: maxpss = memsum_njob_nthr[njob][nthr][1]/1000
            else: maxpss = 0
            print( '%4d %4d %11.5f    %9d %11.5f   %11.5f'%( njob, nthr, tput / 1E6, njob*nthr, tput / tput1, maxpss ) )
                
#------------------------------------------------------------------------------

# Compare MT and ST for many njobs (eg with no simd)
def plot1(rundir, debug=False):
    processFiles(rundir, debug=debug)
    # Create figure with two plots
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(12,6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    # First plot: throughput
    ax1.set_xlabel('Level of parallelism (#jobs times #threads per MT job)')
    ax1.set_ylabel('Node throughput (E6 events per second)')
    for nthr in sorted(allnthr):
        npars = []
        tputs = []
        for njob in sorted(allnjob):            
            if njob not in timsum_njob_nthr: continue
            if nthr not in timsum_njob_nthr[njob]: continue
            npar = nthr*njob
            tput = timsum_njob_nthr[njob][nthr][4] / timsum_njob_nthr[njob][nthr][3] / 1E6
            npars.append(npar)
            tputs.append(tput)
        if nthr == 1 : ax1.plot(npars, tputs, marker='*', label='1 (ST)', ms=20, mew=1)
        else: ax1.plot(npars, tputs, marker='o', label=str(nthr))
    # Second plot: memory
    ax2.set_xlabel('Level of parallelism (#jobs times #threads per MT job)')
    ax2.set_ylabel('Total PSS memory (GB)')
    for nthr in sorted(allnthr):
        npars = []
        mpsss = []
        for njob in sorted(allnjob):            
            if njob not in memsum_njob_nthr: continue
            if nthr not in memsum_njob_nthr[njob]: continue
            npar = nthr*njob
            mpss = memsum_njob_nthr[njob][nthr][1]/1000
            npars.append(npar)
            mpsss.append(mpss)
        if nthr == 1 : ax2.plot(npars, mpsss, marker='*', label='1 (ST)', ms=20, mew=1)
        else: ax2.plot(npars, mpsss, marker='o', label=str(nthr))
    # Decorate both plots 
    title='#threads\n per MT job'
    #loc = 'upper left'
    loc = 'lower right'
    ax1.legend(loc=loc, title=title)
    ax2.legend(loc=loc, title=title)
    node = 'pmpe'
    if node == 'pmpe' :
        ftitle = 'check.exe scalability on pmpe04 (2x 8-core 2.4GHz Haswell with 2x HT)'
        fig.suptitle(ftitle)
        nodet = 'WITHOUT SIMD'
        ax1.set_title(nodet)
        ax2.set_title(nodet)
        xmax=54
        ymax1=25
        ymax2=80
        ax1.axis([0,xmax,0,ymax1])
        ax2.axis([0,xmax,0,ymax2])
        xht=16
        ax1.axvline(xht, color='black', ls=':')
        ax1.axvline(xht*2, color='black', ls='-.')
        ax1.text(xht/2, 0.92*ymax1, 'No HT', ha='center', va='center', size=15)
        ax1.text(xht*3/2, 0.92*ymax1, '2x HT', ha='center', va='center', size=15)
        ax1.text(xmax/2+xht, 0.92*ymax1, 'Overcommit', ha='center', va='center', size=15)
        ax2.axhline(y=64, color='black', ls='-')
        ax2.text(xmax/2+xht/2, 64*1.05, 'MAXIMUM MEMORY: 64 GB', ha='center', va='center', size=12)
        ax2.axvline(xht, color='black', ls=':')
        ax2.axvline(xht*2, color='black', ls='-.')
        ax2.text(xht/2, 0.92*ymax2, 'No HT', ha='center', va='center', size=15)
        ax2.text(xht*3/2, 0.92*ymax2, '2x HT', ha='center', va='center', size=15)
        ax2.text(xmax/2+xht, 0.92*ymax2, 'Overcommit', ha='center', va='center', size=15)
    # Save and show the figure
    # NB: savefig may issue WARNING 'Unable to parse the pattern' (https://bugzilla.redhat.com/show_bug.cgi?id=1653300)
    print( 'Please ignore the warning "Unable to parse the pattern"' )
    save = rundir + '/' + node + '-nosimd.png'
    fig.savefig(save, format='png', bbox_inches="tight")
    from subprocess import Popen
    ###Popen(['eog', '-w', save])
    Popen(['display', save])
    print( 'Plot successfully completed' )

#---------------------------------------

# Compare various simd ST options for many njobs
def plot2(rundir, avxs, debug=False):
    # Loop through files for different avx
    global timsum_njob_nthr
    global memsum_njob_nthr
    timsum_avx_njob_nthr = {} # x[avx,njob,nthr] = [avgtot, avgtot1, avgtot2, avgtot3, sumnevt] -> avg times and sum nevt over all jobs
    memsum_avx_njob_nthr = {} # x[avx,njob,nthr] = [sumrss, sumpss, sumswap] -> sum over all jobs (and threads/processes of MT/MP jobs)
    allnjob2 = []
    allnthr2 = []
    for avx in avxs:
        processFiles(rundir+'.'+avx, avx, debug)
        timsum_avx_njob_nthr[avx] = timsum_njob_nthr
        memsum_avx_njob_nthr[avx] = memsum_njob_nthr
        allnjobtmp = sorted(timsum_njob_nthr)
        allnjob2 += allnjobtmp
        for njob in allnjobtmp: allnthr2 += sorted(timsum_njob_nthr[njob])
    allnjob2 = list(set(allnjob2)) # unique items
    allnthr2 = list(set(allnthr2)) # unique items
    # Create figure with two plots
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(12,6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    # First plot: throughput
    nthr = 1 # only ST
    tput1_avx = {} # throughput with 1 ST job for different AVX
    tcol1_avx = {} # colors of throughput plots for different AVX
    ax1.set_xlabel('Level of parallelism (number of ST jobs)')
    ax1.set_ylabel('Node throughput (E6 events per second)')
    for avx in avxs:
        timsum_njob_nthr = timsum_avx_njob_nthr[avx]
        npars = []
        tputs = []
        for njob in sorted(allnjob2):
            if njob not in timsum_njob_nthr: continue
            if nthr not in timsum_njob_nthr[njob]: continue
            npar = nthr*njob
            tput = timsum_njob_nthr[njob][nthr][4] / timsum_njob_nthr[njob][nthr][3] / 1E6
            npars.append(npar)
            tputs.append(tput)
            if njob == 1: tput1_avx[avx] = tput
        p = ax1.plot(npars, tputs, marker='o', label=avx)
        tcol1_avx[avx] = p[0].get_color()
    # Second plot: memory
    nthr = 1 # only ST
    ax2.set_xlabel('Level of parallelism (number of ST jobs)')
    ax2.set_ylabel('Total PSS memory (GB)')
    for avx in avxs:
        memsum_njob_nthr = memsum_avx_njob_nthr[avx]
        npars = []
        mpsss = []
        for njob in sorted(allnjob2):
            if njob not in memsum_njob_nthr: continue
            if nthr not in memsum_njob_nthr[njob]: continue
            npar = nthr*njob
            mpss = memsum_njob_nthr[njob][nthr][1]/1000
            npars.append(npar)
            mpsss.append(mpss)
        ax2.plot(npars, mpsss, marker='o', label=avx)
    # Decorate both plots 
    title='SIMD mode'
    loc = 'right'
    ax1.legend(loc=loc, title=title)
    loc = 'lower right'
    ax2.legend(loc=loc, title=title)
    node = 'pmpe'
    if node == 'pmpe' :
        ftitle = 'check.exe scalability on pmpe04 (2x 8-core 2.4GHz Haswell with 2x HT)'
        fig.suptitle(ftitle)
        nodet = 'VARIOUS SIMD MODES'
        ax1.set_title(nodet)
        ax2.set_title(nodet)
        xmax=54
        ymax1=110
        ymax2=80
        ax1.axis([0,xmax,0,ymax1])
        ax2.axis([0,xmax,0,ymax2])
        xht=16
        ax1.axvline(xht, color='black', ls=':')
        ax1.axvline(xht*2, color='black', ls='-.')
        ax1.text(xht/2, 0.92*ymax1, 'No HT', ha='center', va='center', size=15)
        ax1.text(xht*3/2, 0.92*ymax1, '2x HT', ha='center', va='center', size=15)
        ax1.text(xmax/2+xht, 0.92*ymax1, 'Overcommit', ha='center', va='center', size=15)
        for avx in avxs:
            tput1 = tput1_avx[avx]
            tcol1 = tcol1_avx[avx]
            ax1.plot( [0,xht], [0,xht*tput1], marker='', ls=':', lw=2, color=tcol1 )
        ax2.axhline(y=64, color='black', ls='-')
        ax2.text(xmax/2+xht/2, 64*1.05, 'MAXIMUM MEMORY: 64 GB', ha='center', va='center', size=12)
        ax2.axvline(xht, color='black', ls=':')
        ax2.axvline(xht*2, color='black', ls='-.')
        ax2.text(xht/2, 0.92*ymax2, 'No HT', ha='center', va='center', size=15)
        ax2.text(xht*3/2, 0.92*ymax2, '2x HT', ha='center', va='center', size=15)
        ax2.text(xmax/2+xht, 0.92*ymax2, 'Overcommit', ha='center', va='center', size=15)
    # Save and show the figure
    # NB: savefig may issue WARNING 'Unable to parse the pattern' (https://bugzilla.redhat.com/show_bug.cgi?id=1653300)
    print( 'Please ignore the warning "Unable to parse the pattern"' )
    save = rundir + '.none/' + node + '-simd.png'
    fig.savefig(save, format='png', bbox_inches="tight")
    from subprocess import Popen
    ###Popen(['eog', '-w', save])
    Popen(['display', save])
    print( 'Plot successfully completed' )

#---------------------------------------

if __name__ == '__main__':

    ###parseMemXml('BMKTST/check-test.none.j016.t001/1/mem.txt', debug=False)

    ###parseLogTxt('BMKTST/check-test.none.j016.t001/1/log.txt')

    ###processFiles('BMKTST', debug=True)
    ###processFiles('BMKTST', debug=False)
    ###processFiles('BMKTST.sse4', avx='sse4', debug=False)
    ###processFiles('BMKTST.avx2', avx='avx2', debug=True)

    #plot1('BMKTST', debug=False)

    plot2('BMKTST', ['none', 'sse4', 'avx2'], debug=False)
