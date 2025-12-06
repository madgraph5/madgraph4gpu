#!/bin/env python3
import sys
if len(sys.argv) != 2:
    print('Usage:', sys.argv[0], '<filename>')
    sys.exit(1)
filename=sys.argv[1]
tjcm_bld_fp={}
with open(filename) as file:
    for line in file:
        ###print(line.rstrip())
        lsplit=line.rstrip().split()
        if line.startswith('PROC'):
            fp=lsplit[1].replace('FPTYPE=','')
            bld=lsplit[2].replace('BLD=','')
            if fp not in tjcm_bld_fp: tjcm_bld_fp[fp]={}
            if bld not in tjcm_bld_fp[fp]: tjcm_bld_fp[fp][bld]=[]
        elif len(lsplit)>1 and lsplit[1]=='Jamps':
            tjcm_bld_fp[fp][bld].append(lsplit[7])
            tjcm_bld_fp[fp][bld].append(lsplit[5])
        elif len(lsplit)>1 and lsplit[1]=='ColorSum':
            tjcm_bld_fp[fp][bld].append(lsplit[5])
        elif len(lsplit)>1 and lsplit[1]=='MeanMatrixElemValue':
            tjcm_bld_fp[fp][bld].append(lsplit[3].replace('e-07',''))
###for fp in tjcm_bld_fp:
for fp in ('d','m','f'): # reorder
    bs,ts,js,cs,ms=[],[],[],[],[]
    for bld in tjcm_bld_fp[fp]:
        #print(fp, bld, tjcm_bld_fp[fp][bld])
        bs.append(bld)
        ts.append(tjcm_bld_fp[fp][bld][0])
        js.append(tjcm_bld_fp[fp][bld][1])
        cs.append(tjcm_bld_fp[fp][bld][2])
        ms.append(tjcm_bld_fp[fp][bld][3])
    print('FPTYPE=%s'%fp)
    print('BLD   ', '  '.join('%-8s'%v for v in bs))
    print('Total ', '  '.join('%-8s'%v for v in ts))
    print('Jamps ', '  '.join('%-8s'%v for v in js))
    print('ColSum', '  '.join('%-8s'%v for v in cs))
    print('MeanME', '  '.join('%-8s'%v for v in ms))
    print()
