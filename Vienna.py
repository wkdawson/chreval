#!/usr/bin/python

# Tools for handles Vienna package data with SimRNA

# This can be tested in the following way
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# Python 2.7.6 (default, Mar 22 2014, 22:59:56) 
# [GCC 4.8.2] on linux2
# Type "help", "copyright", "credits" or "license" for more information.
# >>> import Vienna
# >>> vs = Vienna.Vstruct()
# >>> vs.convert_Vstruct("(((((....)))))")
# (((((....)))))
# (    0,   13)
# (    1,   12)
# (    2,   11)
# (    3,   10)
# (    4,    9)
# 0
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

import sys
debug = False
# debug = True

TEST = 2

# used for PKs and parallel stems
# this notation at least works with VARNA
pairing = {  0 : ['(', ')'],
             1 : ['[', ']'],
             2 : ['{', '}'],
             3 : ['<', '>'],
             4 : ['A', 'a'],
             5 : ['B', 'b'],
             6 : ['C', 'c'],
             7 : ['D', 'd'],
             8 : ['E', 'e'],
             9 : ['F', 'f'],
            10 : ['G', 'g'],
            11 : ['H', 'h'],
            12 : ['I', 'i'],
            13 : ['J', 'j'],
            14 : ['K', 'k'],
            15 : ['L', 'l'],
            16 : ['M', 'm'],
            17 : ['N', 'n'],
            #18 : ['O', 'o'], # doesn't work beyond N
            #19 : ['P', 'p'],
            #20 : ['Q', 'q'],
            #21 : ['R', 'r'],
            #22 : ['S', 's'],
            #23 : ['T', 't'],
            #24 : ['U', 'u'],
            #25 : ['V', 'v'],
            #26 : ['W', 'w'],
            #27 : ['X', 'x'],
            #28 : ['Y', 'y'],
            #29 : ['Z', 'z'],
}


# error handling 
class MyException(Exception):
   def __init__(self, st):
      self.stmnt = st
   #
   def __str__(self):
      return self.stmnt
   # Example:
   # if not self.set_var:
   #    mssg = "\nERROR, value of var undefined\n"
   #    ex=MyException(mssg)
   #    raise ex
#



class BP:
   def __init__(self):
      self.i = -1
      self.j = -1
   #
   def put_BP(self, i, j):
      self.i = i
      self.j = j
      return 0
   #
   def put_BPi(self, i):
      self.i = i
      return 0
   #
   def put_BPj(self, j):
      self.j = j
      return 0
   #
   def print_BP(self):
      print '(%5d,%5d)' % (self.i, self.j)
      return 0
   #       
#        

class Vstruct:
   def __init__(self):
      self.vstr=''
      self.N = -1
      # standard secondary structure (ss)
      self.BPlist=[]
      self.stack = 0
      self.point = 0
      self.jcount = 0
      # a both ss and pseudoknot (PK) in same line
      self.PKlist=[]
      self.hstack = 0
      self.hpoint = 0
      self.jhcount = 0
      # CTCF islands
      self.CTCFlist = []
      self.CTCFpoint = []
   #
   def set_Vstruct(self, s):
      self.vstr = s
      self.N = len(s)
      return 0
   #
   #
   def reset_Vstruct(self):
      self.vstr=''
      self.N = -1
      # standard secondary structure (ss)
      self.BPlist=[]
      self.stack = 0
      self.point = 0
      self.jcount = 0
      # both ss and pseudoknot (PK) in same line
      self.PKlist=[]
      self.hstack = 0
      self.hpoint = 0
      self.jhcount = 0
      # CTCF islands
      self.CTCFlist = []
      self.CTCFpoint = []
      return 0
   # 
   def scan_vs(self, i):
      if debug:
         print "scan_vs(%d):" % i
      #
      mm = -1
      if i < self.N:
         s = self.vstr[i]
         if not (s == '(' or s == ')' or s == '[' or s == ']' or s == '.'):
            print 'ERROR: improperly defined Fontana formated sequence'
            print '       offending character \'%c\' located at position %d' % (s, i+1)
            sys.exit(1)
         #
         if s == '(':
            # add to the stack
            b = BP()
            b.put_BPi(i)
            self.BPlist += [b]
            if debug:
               b.print_BP()
               # print 'i = %5d, %s' % (i, s)
            mm = i + 1
            self.stack += 1 
            self.point = self.stack
            self.jcount +=1
            self.scan_vs(mm)
         elif s == ')':
            self.point = self.find_next_BPpoint(self.point)
            self.BPlist[self.point-1].put_BPj(i)
            if debug:
               self.BPlist[self.point-1].print_BP()
               # print 'j = %5d, %s' % (i, s)
            self.point -= 1
            self.jcount -= 1
            mm = i + 1
            self.scan_vs(mm)
         else:
            mm = i
            # ignore pseudoknot brackets [[[...]]], if present.
            # Converts these to unpaired strand data.
            if (s == '[') or (s == ']'):
               s = '.'
            # Either look for next closing bracket in sequence or
            # terminate at the end of the sequence if nothing is
            # found.
            while (s == '.') and (mm < self.N):
               if debug:
                  print 'i = %5d, %s' % (mm, s)
               mm += 1
               if mm == self.N:
                  break
               s = self.vstr[mm]
               
            self.scan_vs(mm)
            #
         #
      #
      else:
         if debug:
            print 'point = %d' % self.point
         if not self.jcount == 0:
            case = 0
            if self.jcount > 0:
               case = 0
            else:
               case = 1
            print 'jcount = %d' %  self.jcount
            print 'ERROR!!! Fontana notation is not correct!'
            if case == 0:
               print '         %d too many \'(\' brackets!' % self.jcount
            else:
               print '         %d too many \')\' brackets!' % (-self.jcount)
            #
            print 'input structure:'
            print self.vstr
            sys.exit(1)
         if not self.check_BP():
            print 'ERROR!!! Fontana notation is not correct!'
            print '         At least one structure of type \')...(\' was found'
            print 'input structure:'
            print self.vstr
            sys.exit(1)
                
      return mm
   #
   
   
   def scan_ctcf(self, i):
      if debug:
         print "scan_ctcf(%d):" % i
      #
      
      mm = -1
      if i < self.N:
         s = self.vstr[i]
         if not (s == '(' or s == ')' or s == '[' or s == ']' or s == '.' or s == '|'):
            print 'ERROR: improperly defined Fontana formated sequence'
            print '       offending character \'%c\' located at position %d' % (s, i+1)
            sys.exit(1)
         #
         if s == '(':
            # add to the stack
            b = BP()
            b.put_BPi(i)
            self.BPlist += [b]
            if debug:
               b.print_BP()
               # print 'i = %5d, %s' % (i, s)
            mm = i + 1
            self.stack += 1 
            self.point = self.stack
            self.jcount +=1
            self.scan_ctcf(mm)
         elif s == ')':
            self.point = self.find_next_BPpoint(self.point)
            self.BPlist[self.point-1].put_BPj(i)
            if debug:
               self.BPlist[self.point-1].print_BP()
               # print 'j = %5d, %s' % (i, s)
            self.point -= 1
            self.jcount -= 1
            mm = i + 1
            self.scan_ctcf(mm)
         elif s == '[':
            if debug:
               print "pk ["
            #
            mm = i
            b = BP()
            b.put_BPi(i)
            self.PKlist += [b]
            if debug:
               b.print_BP()
               print 'i = %5d, %s' % (i, s)
            mm = i + 1
            self.hstack += 1 
            self.hpoint = self.hstack
            self.jhcount +=1
            self.scan_ctcf(mm)
         elif s == ']':
            if debug:
               print "pk ]"
            #
            self.hpoint = self.find_next_PKpoint(self.hpoint)
            self.PKlist[self.hpoint-1].put_BPj(i)
            if debug:
               self.PKlist[self.hpoint-1].print_BP()
               print 'j = %5d, %s' % (i, s)
            self.hpoint -= 1
            self.jhcount -= 1
            mm = i + 1
            self.scan_ctcf(mm)
            #
         elif s == '|':
            if debug:
               print "ctcf |"
            #
            self.CTCFpoint += [i]
            mm = i + 1
            if debug:
               print 'j = %5d, %s' % (i, s)
            #
            self.scan_ctcf(mm)
            #
         else:
            # exit out the non-interacting point
            mm = i
            while (s == '.') and (mm < self.N):
               if debug:
                  print 'i = %5d, %s' % (mm, s)
               mm += 1
               if mm == self.N:
                  break
               s = self.vstr[mm]
               
            self.scan_ctcf(mm)
            #
         #
      #
      else:
         if debug:
            print 'point = %d' % self.point
         if not self.jcount == 0:
            case = 0
            if self.jcount > 0:
               case = 0
            else:
               case = 1
            print 'jcount = %d' %  self.jcount
            print 'ERROR!!! Fontana notation is not correct!'
            if case == 0:
               print '         %d too many \'(\' brackets!' % self.jcount
            else:
               print '         %d too many \')\' brackets!' % (-self.jcount)
            #
            print 'input structure:'
            print self.vstr
            sys.exit(1)
         if not self.jhcount == 0:
            case = 0
            if self.jhcount > 0:
               case = 0
            else:
               case = 1
            print 'jhcount = %d' %  self.jhcount
            print 'ERROR!!! Fontana notation is not correct!'
            if case == 0:
               print '         %d too many \'[\' brackets!' % self.jhcount
            else:
               print '         %d too many \']\' brackets!' % (-self.jhcount)
            #
            print 'input structure:'
            print self.vstr
            sys.exit(1)
         if not self.check_BP():
            print 'ERROR!!! Fontana notation is not correct!'
            print '         At least one structure of type \')...(\' was found'
            print 'input structure:'
            print self.vstr
            sys.exit(1)
         if len(self.CTCFpoint) > 0:
            if debug:
               print "found CTCF references:"
            #
            if len(self.CTCFpoint) < 2:
               print "there must be more than one CTCF marker to assign this structure"
               sys.exit(1)
            self.CTCFpoint.sort()
            for l in range(1, len(self.CTCFpoint)):
               for k in range(0,l):
                  b = BP()
                  b.put_BPi(self.CTCFpoint[k])
                  b.put_BPj(self.CTCFpoint[l])
                  self.CTCFlist += [b]
               #
            #
         #
      #
      return mm
   #
   
   
   def find_next_BPpoint(self, p):
      pp = p
      j = self.BPlist[p-1].j
      # print 'find_next_BPpoint: j = %d, pp = %d' % (j, pp)
      while j > -1 and pp >= 0:
         pp -= 1
         j = self.BPlist[pp-1].j
      self.point = pp
      return pp
   #
   
   def find_next_PKpoint(self, p):
      pp = p
      j = self.PKlist[p-1].j
      # print 'find_next_PKpoint: j = %d, pp = %d' % (j, pp)
      while j > -1 and pp >= 0:
         pp -= 1
         j = self.PKlist[pp-1].j
      self.hpoint = pp
      return pp
   #
   
   def check_BP(self):
      pass_BP = True
      for b in self.BPlist:
         if b.j < 0:
            pass_BP = False
         #
      #
      return pass_BP
   #
   def print_vstr(self):
      print self.vstr
   
   def print_BPlist(self):
      # left for backward compatibility
      for b in self.BPlist:
         b.print_BP()
      #
   #
    
   def print_pairlist(self, plist):
      for b in plist:
         b.print_BP()
      #
   #
    
   def convert_Vstruct(self, s):
      self.set_Vstruct(s)
      self.scan_vs(0)
      self.print_vstr()
      self.print_pairlist(self.BPlist)
      return 0
   #

   def convert_CTCFstruct(self, s):
      self.set_Vstruct(s)
      self.scan_ctcf(0)
      self.print_vstr()
      print "secondary structure"
      self.print_pairlist(self.BPlist)
      if len(self.PKlist) > 0:
         print "pk connects"
         self.print_pairlist(self.PKlist)
      if len(self.CTCFlist) > 0:
         print "CTCF connects"
         self.print_pairlist(self.CTCFlist)
      return 0
   #

#


# reads in Vienna package data created by the SimRNA_trafl2pdbs
# processing function
class ViennaData:
   #
   def __init__(self):
      # secondary structure created from SimRNA_trafl2pdbs processing
      self.ss_flnm = ''
      self.ss_data = []
      self.n_structs = -1  
   #
   
   ################################
   ########   Functions    ########  
   ################################
   
   def reset_ViennaData(self):
      self.ss_flnm = ''
      self.ss_data = []
      self.n_structs = -1  
      return 0
   #
       
   # reading in trajectory data (*.trafl) from SimRNA run
   def get_ViennaData(self, fl):
      self.ss_flnm = fl
      try: 
         input = open(self.ss_flnm, 'r')
      except:
         print 'ERROR: %s does not exist!' % self.ss_flnm
         sys.exit(1)
      #
      while True:
         string = input.readline()        # ((((.....)))) data
         string = string[:len(string)-1]  # remove the '\n'
         if not string: break
         self.ss_data += [string]
      #
      input.close()
      self.n_structs   = len(self.ss_data)
      if self.n_structs == 0:
         mssg = '\nERROR: no data read\n'
         raise MyException(mssg)
      return 0
   #

   # reading in trajectory data (*.trafl) from SimRNA run
   def read_SimRNA_ss(self, data):

      # print "read_SimRNA_ss(): \n%s" % data

      self.ss_data = data
      self.n_structs   = len(self.ss_data)
      if self.n_structs == 0:
         mssg = '\nERROR: no data read\n'
         raise MyException(mssg)
      # print self.n_structs, self.ss_data
      return 0
   #


#

def test1(cl):
   vs = Vstruct()
   ss_seq = "(((((....)))))"
   if len(cl) > 1:
      ss_seq = cl[1]
   vs.convert_Vstruct(ss_seq)
   # output 
   # (((((....)))))
   # (    0,   13)
   # (    1,   12)
   # (    2,   11)
   # (    3,   10)
   # (    4,    9)
   # 0
#   

def test2(cl):
   vs = Vstruct()
   ss_seq = "|(((((.[[...)))))..]]...|..([)]([)]...|"
   if len(cl) > 1:
      ss_seq = cl[1]
   vs.convert_CTCFstruct(ss_seq)
   # output 
   # |(((((.[[...)))))..]]...|..([)]([)]...|
   # secondary structure
   # (    1,   16)
   # (    2,   15)
   # (    3,   14)
   # (    4,   13)
   # (    5,   12)
   # (   27,   29)
   # (   31,   33)
   # pk connects
   # (    7,   20)
   # (    8,   19)
   # (   28,   30)
   # (   32,   34)
   # CTCF connects
   # (    0,   24)
   # (    0,   38)
   # (   24,   38)
#


def main(cl):
   if TEST == 0 or TEST == 1:
      test1(cl)
   elif TEST == 2:
      test2(cl)
   #
#
   
# Main
if __name__ == '__main__':
   main(sys.argv)

