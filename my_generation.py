#!/usr/bin/env python

# my_generation.py: generates a heatmap for reading by chreval.py.
# This is intended to be used to build test structures for chreval.py
# to find. It is particularly needed for pseudoknot like contraptions
# where all sorts of things hairy things could go wrong.

from Vienna import Vstruct
from MatrixTools import MatrixTools
import sys
import random
import argparse
import os

def add_noise(N, hv):
    for j in range(0, N):
        for i in range(0,j):
            u1 = int(random.uniform(0,1) + 0.5)
            u2 = int(random.uniform(0,1) + 0.5)
            v = u1*u2*int(random.uniform(0,2) + 0.5)
            hv[i][j] = v
            hv[j][i] = v
        #
    #
    return hv
#

def read_SeqFile(iflnm):
    print iflnm
    if not os.path.isfile(iflnm):
        print "ERROR: %s not found" % iflnm
        sys.exit(1)
    #
    fp = open(iflnm, 'r')
    lfp = fp.readlines()
    fp.close()
    ss_seq = []
    for k in range(0, len(lfp)):
        s = lfp[k].strip().split()
        if len(s) > 0:
            if s[0][0] == '#':
                continue
            ss_seq += [s[0]]
        #
    #
    for ss in ss_seq:
        print ss
    return ss_seq
#
        
def generate(ss_seq, w_noise, oflnm):
    vs = []
    k = 0
    N = len(ss_seq[0])
    
    # setup objects and check input sequence (or sequences)
    for k in range(0, len(ss_seq)):
        vs += [Vstruct()]
        n = len(ss_seq[k])
        if not n == N:
            print "ERROR: input sequence is not the same length as the others"
            m = 1
            for m in range(0, len(ss_seq)):
                print "seq(%2d): %s" % (m, ss_seq[m])
            print "new seq: %s" % ss_seq[k]
            sys.exit(1)
    #
    
    for k in range(0, len(ss_seq)):
        vs[k].set_Vstruct(ss_seq[k])
        vs[k].scan_ctcf(0)
        vs[k].print_vstr()
    
        # display for information purposes
        if len(vs[k].BPlist) > 0:
            print "secondary structure"
            vs[k].print_pairlist(vs[k].BPlist)
        #
        if len(vs[k].PKlist) > 0:
            print "pk connects"
            vs[k].print_pairlist(vs[k].PKlist)
        #
        if len(vs[k].CTCFlist) > 0:
            print "CTCF connects"
            vs[k].print_pairlist(vs[k].CTCFlist)
        #
    #
    # sys.exit(0)
    mtools = MatrixTools()
    hv = []
    hv = mtools.initialize_matrix(hv, N, 0)
    
    # add noise if requested
    if w_noise:
        hv = add_noise(N, hv)
    #
    for k in range(0, len(vs)):
        for pair in vs[k].BPlist:
            i = pair.i
            j = pair.j
            hv[i][j] = 3
            hv[j][i] = hv[i][j]
        #
        for pair in vs[k].PKlist:
            i = pair.i
            j = pair.j
            hv[i][j] = 5
            hv[j][i] = hv[i][j]
        #
        for pair in vs[k].CTCFlist:
            i = pair.i
            j = pair.j
            hv[i][j] = 100
            hv[j][i] = hv[i][j]
        #
    #
    
    heatmap = mtools.make_heatmap(hv)
    print heatmap
    fp = open(oflnm, 'w')
    fp.write(heatmap)
    fp.close()
#
    



def main(cl):
    
    
    # parse the command line
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-o', nargs=1, default=["test.heat"],
                        dest='f_heatmap',
                        help="name of output heatmap file [default 'test.heat']")
    
    parser.add_argument('-add_noise', action='store_true', default=False,
                        dest='w_noise',
                        help='Request to include noise in the output heatmap.')
    
    input_opt = parser.add_mutually_exclusive_group()    
    input_opt.add_argument('-f', nargs=1, default=["test.ss"],
                        dest='f_ss',
                        help="name of input 1D structure file [default 'test.ss']")
    
    input_opt.add_argument('-seq',     action='store', nargs='+', default=None,
                        dest='seq', 
                        help='The input structure (in Fontana-like format).')
    
    
    #
    # assign arguments
    args = parser.parse_args()
    # print args
    ss_seq = args.seq
    oflnm = args.f_heatmap[0]
    iflnm = args.f_ss[0]
    w_noise = args.w_noise
    print 'iflnm ', iflnm
    print 'oflnm ', oflnm
    
    ss_seq = read_SeqFile(iflnm)
    generate(ss_seq, w_noise, oflnm)
#


# Main
if __name__ == '__main__':
    main(sys.argv)

