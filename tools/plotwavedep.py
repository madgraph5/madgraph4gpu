#!/usr/bin/env python3

import re
import sys
from optparse import OptionParser


class Parse():

    def __init__(self, filename, title):
        self.title = title
        fh = open(filename)
        self.filedat = fh.readlines()
        fh.close()
        self.data = []
        self.labels = {}

    def genSig(self, str):
        repl = ('[', ']', '\n')
        for r in repl:
            str = str.replace(r, '')
        label = str
        if label[0] == 'a':
            label = label.replace('a', 'amplitude ')
        elif label[0] == 'm':
            label = label.replace('m', 'momentum ')
        repl = (',', ')')
        for r in repl:
            str = str.replace(r, '')
        repl2 = ('(')
        for r in repl2:
            str = str.replace(r, '_')
        if str not in self.labels:
            self.labels[str] = label
        return str

    def parse(self):
        for l in self.filedat:
            ls = l.split('(')
            args = ls[-1]
            while args[-1] in (')', '\n'):
                args = args[:-1]
            argl = args.split(',')
            argl = [x.strip() for x in argl]
            self.data.append([ls[0], argl[:-1], argl[-1], self.genSig(l), []])

        # rets = []
        # for x in self.data:
        #     ret = x[2]
        #     if ret in rets:
        #         print "double return value", ret
        #     else:
        #         rets.append(ret)

        self.data.reverse()

        ipos = 0
        for x in self.data:
            arg = x[1]
            ret = x[2]
            sig = x[3]

            if ret[:1] == 'a':
                self.data.append([ret, [], '', self.genSig(ret), [sig]])

            if len(arg) == 1 and arg[0][:1] == 'm':
                x[4].append(arg[0])
            else:
                for a in arg:
                    for y in self.data[ipos:]:
                        if a == y[2]:
                            x[4].append(y[3])
                            break
            ipos += 1

        deps = []
        for x in self.data:
            dep = x[4]
            sig = x[3]
            for d in dep:
                deps.append(self.genSig(d) + ' -- ' + sig)

        print ('graph gg_ttxgg{')

        if self.title:
            print('    label="%s";' % self.title)
            print('    labelloc="t";')
            print('    fontsize=30;')

        for l in self.labels:
            print('    %s [label="%s"];' % (l, self.labels[l]))

        for d in deps:
            print('    %s;' % d)

        print ('}')

    def clean(self):
        localdat = []
        tmpline = ''
        for l in self.filedat:
            l = l.strip()
            if l[:2] == '//':
                continue
            elif l[-1] != ';':
                tmpline = l
            elif tmpline != '' and l[-1] == ';':
                tmpline = tmpline + l.strip()
                localdat.append(tmpline)
                tmpline = ''
            else:
                localdat.append(l)

        for i in range(len(localdat)):
            repl = [' ', '&', '0.,', '+1,', '-1,', ';', '\n']
            for r in repl:
                localdat[i] = localdat[i].replace(r, '')
            repl2 = [['local_mom', 'm'], ['amp', 'a']]
            for r in repl2:
                localdat[i] = localdat[i].replace(r[0], r[1])
            res = ['thrust::complex<double>\(cIPC\[[0-9]+\],cIPC\[[0-9]+\]\),',
                   'cIPD\[[0-9]+\],',
                   'cHel\[ihel\]\[[0-9]+\],']
            for r in res:
                localdat[i] = re.sub(r, '', localdat[i])

        self.filedat = localdat

        # for x in self.filedat:
        #     print x

    def run(self):
        self.clean()
        self.parse()


class MyParser(OptionParser):
    def format_epilog(self, formatter):
        return self.epilog


if __name__ == '__main__':

    epilog = """
Argument:
  The script accepts one argument which is the filename containing the
  verbatim C++ functions generated usually in calculate_wavefunctions

  e.g.
  [...]
  vxxxxx(local_mom[4], 0., cHel[ihel][4], +1, w[4]);
  VVV1P0_1(w[0], w[1], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[5]);
  FFV1P0_3(w[3], w[2], thrust::complex<double> (cIPC[2], cIPC[3]), 0., 0.,
      w[6]);
  // Amplitude(s) for diagram number 1
  [...]

How to proceed:
   The script generates graphviz/dot format text. Pipe the output in a file
   and run dot on it to generate your prefered output. E.g.

   dot -Tpng -O <filename.dot>"""
    usage = "usage: %prog [options] filename"
    parser = MyParser(usage=usage, epilog=epilog)
    parser.add_option("-t", "--title", dest="ttl",
                      help="diagram title")

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        sys.exit(1)

    title = ''
    if options.ttl is not None:
        title = options.ttl

    Parse(args[0], title).run()
