from optparse import OptionParser
from datetime import datetime
from readdata import ReadData
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys


class Perf():

    def __init__(self, date, run, x, y, z, xrem, yrem, loc):
        perffile = '%s/%s-perf-test-run%s.txt' % (loc, date, run)
        rd = ReadData(perffile)
        self.axesn = [x, y, z]
        self.axesr = [xrem, yrem]  # remove outer bands from axes
        self.axesv = [[], [], []]
        self.data = rd.genData()
        # print self.data[0]

    def prepAxes(self):
        for d in self.data:
            ks = d.keys()
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

    def prepData(self):
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

    def plot(self):
        self.prepAxes()
        self.prepData()

        data_array = np.array(self.data2d)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x_data, y_data = np.meshgrid(np.arange(data_array.shape[1]),
                                     np.arange(data_array.shape[0]))
        xticks = x_data[0]
        yticks = np.array(range(len(y_data)))
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

    (op, ar) = parser.parse_args()

    xrm = 0
    yrm = 0
    if op.xrm:
        xrm = int(op.xrm)
    if op.yrm:
        yrm = int(op.yrm)

    if (ar):
        print parser.print_help()
        sys.exit(1)

    Perf(op.date, op.run, op.xax, op.yax, op.zax, xrm, yrm, op.dir).plot()
