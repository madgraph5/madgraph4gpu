import re
from datetime import datetime


class ReadData():

    def __init__(self, filename):
        self.fh = open(filename)
        self.matchdate = r'([A-Z,a-z]{3})\s([A-Z,a-z]{3})\s{1,2}([0-9]{1,2})\s([0-9]{2}\:[0-9]{2}\:[0-9]{2})\s([CEST]{3,4})\s([0-9]{4})'
        self.matchscie = r'([0-9]{1}\.[0-9]+e[+-][0-9]{2})'

    def genData(self):
        starttime = None
        endtime = None
        values = []
        localdict = {}
        command = None
        for line in self.fh.readlines():
            if line.find('check.exe') != -1:
                command = line
            if line.find('***********************************') != -1:
                if localdict:
                    values.append(localdict)
                localdict = {}
                localdict['CommandExecuted'] = command
            if line.find(' = ') != -1:
                ls = line.split()
                val = ls[-1]
                match = re.search(self.matchscie, val)
                if val.isdigit():
                    val = int(val)
                elif match:
                    val = float(val)
                localdict[ls[0]] = val
            match = re.search(self.matchdate, line)
            if match:
                g = match.groups()
                dstr = '%s-%s-%s %s' % (g[5], g[1], g[2].rjust(2, '0'), g[3])
                if starttime:
                    pass
                    # print endtime - starttime
                if endtime:
                    starttime = endtime
                endtime = datetime.strptime(dstr, '%Y-%b-%d %H:%M:%S')
                if starttime and endtime:
                    localdict['WallStartTime'] = starttime
                    localdict['WallEndTime'] = endtime
        values.append(localdict)
        return values


if __name__ == '__main__':
    r = ReadData('20200402-perf-test-run2.txt')
    v = r.genData()
    print v
    print len(v)
    print v[-1].keys()
