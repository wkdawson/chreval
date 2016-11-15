#!/usr/bin/python

# anal_CTCF.py 
#
#    This is supposed to be a self-contained data object. Not sure if
#    it really functions as so, but I was able to decouple it from
#    Anal at least.

import sys
import os
import chreval
from difflib import SequenceMatcher
from GetOpts import GetOpts
from ChromatinData import Data

#
PROGRAM = "anal_CTCFs.py"


class Anal:
    def __init__(self, data):
        self.debug              = False
        self.data               = data
        self.cdata              = data.cdata
        self.keylist            = data.keylist
        self.datatype           = data.datatype
        self.res                = "res5kb"
        self.range_TddS         = 2.0
        self.range_dG           = 2.0
    #
    
    def get_PET_wt(self, ky, wt):
        ctcf1 = self.cdata[ky].ctcf1
        ctcf2 = self.cdata[ky].ctcf2
        wt_ctcf = wt
        if (ctcf1 == 'R' and ctcf2 == 'R') or (ctcf1 == 'L' and ctcf2 == 'L'):
            wt_ctcf *= 0.25
        elif (ctcf1 == 'L' and ctcf2 == 'R'):
            wt_ctcf *= 0.05
        #
        return wt_ctcf
    #
    
    def build_CTCF_map(self, cl):
        print "build_CTCF_map()"
        keylist = self.keylist
        if len(self.keylist) == 0:
            print "setting keylist"
            print "length of the dataset: ", len(self.cdata)
            keylist = self.data.ordered_keylist()
        print keylist
        print self.cdata.keys()
        print self.data.chrmsm_grp.keys()
        for kk in self.data.chrmsm_grp.keys():
            for cgk in self.data.chrmsm_grp[kk].chrsegment:
                #
                v = self.data.makekey([kk, cgk[0], cgk[1]])
                print v, self.cdata[v].ctcf1, self.cdata[v].ctcf2 
            # print self.data.chrmsm_grp[kk].show_chr_fragments()
        #
        sys.exit(0)

    


def main(cl):
    
    CL = GetOpts(PROGRAM)
    
    if len(cl) < 2:
        emsg = "ERROR: too few arguments"
        CL.error(emsg)
    #
    
    dt = Data()
    files = CL.f_activity
    print files
    
    
    for f in files:
        if os.path.exists(f):
            dt.readfile(f)
        else:
            print "ERROR: cannot open %s; the file missing" % f
            sys.exit(1)
        #
    #
    dt.ordered_keylist() # make an ordered list in terms of chromosome
                         # number and region
    dt.display_Data(False)
    
    an = Anal(dt)
    an.build_CTCF_map(CL)
    
    
#

if __name__ == '__main__':
    # running the program
    main(sys.argv)

