#!/usr/bin/env python3

from optparse import OptionParser
from datetime import datetime
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import ScalarFormatter
import numpy as np
import copy
import sys
import json
from operator import itemgetter


class Perf():

    def __init__(self, date, run, x, y, z, xrem, yrem, loc):
        perffile = '%s/%s-perf-test-run%s.json' % (loc, date, run)
        data = open(perffile, 'r')
        readJson = json.loads(data.read())
        data.close()
        self.axesn = [x, y, z]
        self.axesr = [xrem, yrem]  # remove outer bands from axes
        self.axesv = [[], [], []]
        self.data = self.prepData(readJson)

    def prepData(self, jsonData):
        for data in jsonData:
            for i in data:
                if isinstance(data[i], type('test')):
                    idx = -1
                    if data[i].find("sec") != -1:
                        idx = data[i].find("sec")
                    elif data[i].find("GEV") != -1:
                        idx = data[i].find("GeV")

                    if idx != -1:
                        data[i] = float(data[i][:idx - 1])
        return jsonData

    def prepAxes3D(self):
        for d in self.data:
            ks = list(d.keys())
            for ax in self.axesn:
                idx = self.axesn.index(ax)
                axlist = self.axesv[idx]
                if ax in ks:
                    axval = d[ax]
                    if axval not in axlist:
                        axlist.append(axval)
                else:
                    print('Error: cannot find axes name %s in %s' % (ax, d))
        if len(self.axesv[0]) * len(self.axesv[1]) != len(self.axesv[2]):
            print('Error: axes don\'t match x * y != z (%d * %d != %d' %
                  (len(self.axesv[0]), len(self.axesv[1]), len(self.axesv[2])))
        self.axesv[0].sort()
        self.axesv[1].sort()
        self.axesv[0] = self.axesv[0][self.axesr[0]:]  # sr
        self.axesv[1] = self.axesv[1][self.axesr[1]:]  # sr

    def prepData3D(self):
        xlen = len(self.axesv[0])
        ylen = len(self.axesv[1])
        self.data2d = []
        ylist = [0] * ylen
        for i in range(xlen):
            self.data2d.append(copy.deepcopy(ylist))
        for d in self.data:
            xpos = -1
            ypos = -1
            if d[self.axesn[0]] in self.axesv[0]:
                xpos = self.axesv[0].index(d[self.axesn[0]])
            if d[self.axesn[1]] in self.axesv[1]:
                ypos = self.axesv[1].index(d[self.axesn[1]])
            if xpos != -1 and ypos != -1:
                zval = d[self.axesn[2]]
                self.data2d[xpos][ypos] = zval

    def plot3D(self):
        self.prepAxes3D()
        self.prepData3D()

        data_array = np.array(self.data2d)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x_data, y_data = np.meshgrid(np.arange(data_array.shape[1]),
                                     np.arange(data_array.shape[0]))
        xticks = x_data[0]
        yticks = np.array(list(range(len(y_data))))
        x_data = x_data.flatten()
        y_data = y_data.flatten()
        z_data = data_array.flatten()
        ax.set_xlabel(self.axesn[1], {'fontsize': 'small'})
        ax.set_xticks(xticks)
        # consider 'fontsize': 'small' for dict also yticklabels
        ax.set_xticklabels(self.axesv[1], {'rotation': 45, 'fontsize': 'small'})
        ax.set_ylabel(self.axesn[0], {'fontsize': 'small'})
        ax.set_yticks(yticks)
        # consider 'fontsize': 'small' for dict
        ax.set_yticklabels(self.axesv[0], {'rotation': 45, 'fontsize': 'small'})
        ax.set_zlabel(self.axesn[2], {'fontsize': 'small'})
        # ax.set_zscale('log')
        # z_data = np.log10(z_data)
        ax.bar3d(x_data, y_data, np.zeros(len(z_data)), 1, 1, z_data)
        plt.show()

    def prepData2D(self):
        self.dataDict2D = {}
        xname = self.axesn[0]
        yname = self.axesn[1]
        zname = self.axesn[2]

        for d in self.data:
            xval = d[xname]
            yval = d[yname]
            zval = d[zname]
            dim = xval * yval
            tick = '%s/%s' % (str(xval), str(yval))
            vallist = [float(str(zval).split()[0]), tick]
            if dim not in self.dataDict2D:
                self.dataDict2D[dim] = [vallist]
            else:
                self.dataDict2D[dim].append(vallist)

    def plot2D(self):
        self.prepData2D()

        cmap = {'32': 'red', '64': 'orange', '128': 'blue', '256': 'green'}
        smap = {'32': 20, '64': 40, '128': 80, '256': 160}

        dims = list(self.dataDict2D.keys())
        dims.sort()
        xlist = list(range(1, len(dims) + 1))
        ylist = []
        clist = []
        slist = []
        ylabels = []
        for d in dims:
            ysublist = []
            for y in self.dataDict2D[d]:
                ysublist.append(y)  # y[0]
            ysublist = sorted(ysublist, key=itemgetter(0), reverse=True)
            clist.append([cmap[x[1].split('/')[0]] for x in ysublist])
            slist.append([smap[x[1].split('/')[0]] for x in ysublist])
            # Temporary conversion for total time for events -> events per sec
            # ysublist[0][0] = d / ysublist[0][0]
            ylabels.append([x[1] for x in ysublist])
            ylist.append([x[0] for x in ysublist])

        fig, ax = plt.subplots()
        print(xlist)
        print(ylist)
        for xe, ye, ce, se in zip(xlist, ylist, clist, slist):
            ax.scatter([xe] * len(ye), ye, s=se, facecolors='none', edgecolors=ce)

        ax.set_xticks(xlist)
        ax.set_xlabel('%s * %s' % (self.axesn[0], self.axesn[1]))
        ax.set_ylabel('%s' % (self.axesn[2]))
        # ax.set_yscale('log')
        ax.set_xticklabels(dims, rotation=45)
        ax.yaxis.set_major_formatter(ScalarFormatter())
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        # Commenting only for the current example due to an overlap of the
        # product labels
        # xpos = 1
        # for y in ylabels:
        #     xstr = ''
        #     for x in y:
        #         # xstr += x.replace('/', '\n')
        #         xstr += x
        #         xstr += '\n'
        #     ax.text(xpos, 1, xstr, {'fontsize': 'xx-small',
        #                             'ha': 'center',
        #                             'va': 'bottom'})
        #     xpos += 1

        import matplotlib.patches as mpatches
        from matplotlib.lines import Line2D

        handlelist = []
        for k in cmap:
            handlelist.append(plt.scatter([],[], s=smap[k], marker='o',
                                          color=cmap[k], facecolor='none'))

        print(handlelist)
        plt.legend(handlelist, [str(x) for x in cmap.keys()], title="# threads / block")


        plt.show()


def print_keys(loc, date, run):
    perffile = '%s/%s-perf-test-run%s.json' % (loc, date, run)
    data = open(perffile, 'r')
    readJson = json.loads(data.read())
    data.close()
    for k in list(readJson[0].keys()):
        print(k)


if __name__ == '__main__':

    n = datetime.now()
    today = str(n.year) + str(n.month).rjust(2, '0') + str(n.day).rjust(2, '0')
    parser = OptionParser()
    parser.add_option('-l', '--location', dest='dir', default='data',
                      help='directory with data (default: data)')
    parser.add_option('-d', '--date', dest='date', default=today,
                      help='date of data files YYYYMMDD (default: today)')
    parser.add_option('-r', '--run', default='1', dest='run',
                      help='run number (default: 1)')
    parser.add_option('-x', dest='xax', default='NumThreadsPerBlock',
                      help='variable name for x axis \
                            (default: NumThreadsPerBlock)')
    parser.add_option('-y', dest='yax', default='NumBlocksPerGrid',
                      help='variable name for y axis \
                            (default: NumBlocksPerGrid)')
    parser.add_option('-z', dest='zax', default='TotalTimeInWaveFuncs',
                      help='variable name for z axis \
                            (default: TotalTimeInWaveFuncs)')
    parser.add_option('--xrm', dest='xrm', default=0,
                      help='# of outer x dimensions to remove')
    parser.add_option('--yrm', dest='yrm', default=0,
                      help='# of outer y dimensions to remove')
    parser.add_option('-k', '--keys', dest='keys', action='store_true',
                      help='print available keys from data')

    (op, ar) = parser.parse_args()

    plotnames = ['2D', '3D']
    plot = '2D'

    xrm = 0
    yrm = 0
    if op.xrm:
        xrm = int(op.xrm)
    if op.yrm:
        yrm = int(op.yrm)

    if op.keys:
        print_keys(op.dir, op.date, op.run)
        sys.exit(0)

    if (len(ar) == 1 and ar[0].upper() not in plotnames) or len(ar) > 1:
        print(parser.print_help())
        sys.exit(1)
    elif len(ar) == 1:
        plot = ar[0].upper()

    p = Perf(op.date, op.run, op.xax, op.yax, op.zax, xrm, yrm, op.dir)
    if plot == '3D':
        p.plot3D()
    if plot == '2D':
        p.plot2D()
