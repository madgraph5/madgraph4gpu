from readdata import ReadData
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import copy


class Perf():

    def __init__(self, date, run, x, y, z, xrem, yrem):
        perffile = '%s-perf-test-%s.txt' % (date, run)
        rd = ReadData(perffile)
        self.axesn = [x, y, z]
        self.axesr = [xrem, yrem]  # remove outer bands from axes
        self.axesv = [[], [], []]
        self.data = rd.genData()
        print self.data[0]

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
            xpos = 0
            ypos = 0
            if d[self.axesn[0]] in self.axesv[0]:
                xpos = self.axesv[0].index(d[self.axesn[0]])
            if d[self.axesn[1]] in self.axesv[1]:
                ypos = self.axesv[1].index(d[self.axesn[1]])
            if xpos and ypos:
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
        x_data = x_data.flatten()
        y_data = y_data.flatten()
        z_data = data_array.flatten()
        # ax.set_zscale('log')
        # z_data = np.log10(z_data)
        ax.bar3d(x_data, y_data, np.zeros(len(z_data)), 1, 1, z_data)
        plt.show()


if __name__ == '__main__':

    Perf('20200402', 'run2', 'NumBlocksPerGrid', 'NumThreadsPerBlock',
         'MeanTimeinWaveFuncs', 0, 0).plot()
#         'TotalTimeInWaveFuncs').plot()
