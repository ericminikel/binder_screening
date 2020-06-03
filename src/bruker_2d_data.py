#!/usr/bin/env python

import sys
import numpy
import argparse

# suggested usage: src/bruker_2d_data.py -p data/2016-02-11-HuPrP23-230-15N-rprp00007-trosy.txt > data/2016-02-11-HuPrP23-230-15N-rprp00007-trosy.matrix

def parse_2d_spectrum(path, dest):
    with open(path, 'r') as f:
        for i in range(0,12):
            line = f.readline()
            if i + 1 <= 8: # first 8 lines will be kept as a header in output
                dest.write(line)
            if i + 1 == 4:
                info = line.strip().split()
                f1left = float(info[3])
                f1right = float(info[7])
            if i + 1 == 5:
                info = line.strip().split()
                f2left = float(info[3])
                f2right = float(info[7])
            if i + 1 == 7:
                info = line.strip().split()
                nrows = float(info[3])
            if i + 1 == 8:
                info = line.strip().split()
                ncols = float(info[3])
                cols = numpy.arange(f2left, f2right, (f2right - f2left)/ncols)
                rows = numpy.arange(f1left, f1right, (f1right - f1left)/nrows)
        dest.write('y\t' + '\t'.join(['x' + x for x in map(str, cols)]))
        row = 0
        for line in f.readlines():
            if line[0] == '#':
                row = int(line.strip().split()[3])
                dest.write('\n' + 'y' + str(rows[row]))
            else:
                dest.write('\t' + line.strip())

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a 2D NMR spectrum exported from Bruker as a text file to a matrix you can manipulate in R')
    parser.add_argument('-p','--path', dest='path',
                       type=str, help='Path to Bruker-exported text file')
    parser.add_argument('-o', '--out', nargs='?', type=argparse.FileType('w'),
                       default=sys.stdout)
    args = parser.parse_args()
    parse_2d_spectrum(path=args.path,dest=args.out)

