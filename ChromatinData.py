#!/usr/bin/python

# ChromatinData.py 
#
#    This is supposed to be a self-contained data object. Not sure if
#    it really functions as so, but I was able to decouple it from
#    Anal at least.

#    The input file is a list of loops formatted according to
#    Przemek's approach (with extension 'bed') Teresa's approach
#    (currently with extension 'txt').

#    command line example:
#    > chreval.py -ff loops.CTCF.annotated.bed -dG

#    where the file "loops.CTCF.annotated.bed" contains a list of data
#    formatted according to Przemek. This has the format

#    file format1 (loops.CTCF.annotated.bed):
#
#                       CTCF    CTCF    PET   cmplx   cmplx   cmplx     
# chr    bgn     end     1       2      cnt    ity      a       b       open            active
# chr1	838908	912011	R	L	11	0	4	6	0.0686291944243	0.426918183932
# chr1	838908	922335	R	R	5	1	6	8	0.0601364066789	0.397329401752
# chr1	838908	1000167	R	L	7	3	15	17	0.0910398799447	0.5333717808
# chr1	918286	969271	R	R	5	0	4	4	0.119015396685	0.543963910954
# chr1	918286	1000167	R	L	75	1	8	8	0.11802493863	0.65697780926
# chr1	918286	1059032	R	R	4	2	12	12	0.09421226891	0.648785755901

#    file format2 (structure.loops.GSE6352.annotated.bed):
# chr1	1050000	1190000	241	-1	9	9	0.0529785714286	0.418914285714
# chr1	1585000	1650000	80	-2	0	0	0.0553384615385	0.737907692308
# chr1	1710000	1840000	154	-2	22	22	0.0699230769231	0.874938461538

#    file format2 (structure.TAD_Rao.annotated.bed):
# chr1	915000	1005000	0.72251	1	10	10	0.116533333333	0.643166666667
# chr1	1030000	1235000	1.1954	-1	12	12	0.0460487804878	0.507936585366
# chr1	1255000	1450000	0.9312	0	27	27	0.0741435897436	0.653230769231


#    In general, bed type files should only be one instance.

#    An alternative file input would involve the following command
#    line

#    command line example:
#    > chreval.py -ff wayn_active_inactive.txt wayn_compartments.txt wayn_open.txt -dG

#    where the various list of "txt" files (wayn_active_inactive.txt,
#    etc.) are formatted according to Teresa.


# used as a identifier                 |used in the program ----------------------------------
# region of chromosome                 |specific location in region             state   length
# chr   bgn             end            |chr     bgn             end              
# chr1	10436116	10556545	chr1	10437400	10438000	active	600
# chr1	10436116	10556545	chr1	10438000	10438466	active	466
# chr1	10436116	10556545	chr1	10438558	10440200	active	1642
        

#    These files often have different data related to open, active,
#    inactive, compartment A or B, etc., so the results typically
#    require more than one file.  Therefore, multiple file entries are
#    allowed for both file types; however, in general, for "bed"
#    files, there should only be one input file.

#    The chromatin loop is specified by the region chr1 10436116 to
#    10556545 and the zones listed above as "active" are bands between
#    (10437400-10438000), (10438000-10438466) and
#    (10438558-10440200). The percentage of active region is the sum
#    of the lengths of the small segments divide by the length of the
#    loop.



import sys
import os
from CUtils import String
from GetOpts import GetOpts

PROGRAM = "anal_loops.py"
string = String(10)

allowed_tags = ['chr', 'bgn', 'end',
                'lCTCF', 'rCTCF',
                'PETcnt', 'cmplx1', 'cmplx2', 'cmplx3',
                'open', 'active', 'A', 'B', 'repressed']


class Segment:
    def __init__(self, bgn, end, state):
        self.bgn    = bgn
        self.end    = end
        self.state  = state
    #
    
    def disp(self):
        return "%10d  %10d  %10d  %s" % (self.bgn, self.end, (self.end - self.bgn), self.state)
    #
    
    def get_length(self):
        length = self.end - self.bgn
        return length
    #
#


# Serves as a container for Stores data on a particular chromosome, or
# part of the chromosome.
class Chromatin:
    """stores records of chromatin segments of a specific type"""
    def __init__(self, name, bgn, end, ctcf1 = '*', ctcf2 = '*', nPET = 0, cmplx = []):
        self.name  = name       # name of the chromosome
        self.bgn   = int(bgn)   # beginning position
        self.end   = int(end)   # ending position
        self.ctcf1 = ctcf1      # CTCF direction at the beginning point
        self.ctcf2 = ctcf2      # CTCF direction at the ending point
        
        # number of PET counts
        if string.isFloat(nPET):
            self.nPET = float(nPET)
        else:
            self.nPET  = int(nPET)
        #
        
        self.cmplx = cmplx      # complexity
        self.state = { "active"    : 0.0,\
                       "open"      : 0.0,\
                       "repressed" : 0.0,\
                       "A"         : 0.0,\
                       "B"         : 0.0 }
        #
        self.segments = [] # mainly for type1 data
    #
    
    
    def add_segments(self, bgn, end, state):
        # print "add_segments"
        self.segments += [Segment(bgn, end, state)]
    #
    
    def disp_head(self):
        return "%5s  %10d  %10d" % (self.name, self.bgn, self.end)
    #
    
    def disp_data(self, fill_region = 0):
        s = ''
        s += "%5s  %10d  %10d      %s      %s" \
             % (self.name, self.bgn, self.end, self.ctcf1, self.ctcf2)
        s += "       p_%s" % string.field(self.nPET, 3)
        k = 1
        for c in self.cmplx:
            s += "    c%d_%s" % (k, string.field(c, 4))
            k += 1
        #
        
        # in case we have to add additional spaces for non existing
        # data and fit it within the same table as some other
        # data. This normally shouldn't be necessary, but there are a
        # variety of data types, and they could be drawn from a
        # variety of sources and mixed together for some reason.
        if fill_region > 0:
            for i in range(0, fill_region):
                s += "  c-_0     "
            #
        #
        return s
    #
#
        
class Chromosome:
    """stores information on the length and contents of each chromosome unit"""
    def __init__(self, name, bgn, end):
        self.name = name
        self.chrsegment = [(int(bgn), int(end))]
    #
    def add_chrsegment(self, bgn, end):
        self.chrsegment += [(int(bgn), int(end))]
    #
    def show_chr_fragments(self):
        s =  self.name + '\n'
        for chrsx in self.chrsegment:
            s += "%10d   %10d\n" % (chrsx[0], chrsx[1])
        return s
    #
#


class Data:
    def __init__(self):
        self.cdata       = {}
        self.chrmsm_grp  = {}
        self.ext         = []
        self.annotation  = []
        self.interaction = []
        self.taglist     = []
        self.description = []
        self.flhd        = []
        self.keylist     = []
        self.datatype    = ''
        self.set_Data    = False
        
        self.version    = '0'
        self.debug_Data = False
        
    #
    
    
    
    # makes a key out of the chromosome region specified
    def makekey(self, s):
        if not(len(s) == 3):
            print "ERROR(makekey): input data does not contain the right number of fields"
            print "                %s" % s
            usage()
            sys.exit(1)
        #
        bgn = int(s[1])
        end = int(s[2])
        
        sp = s[0] + '_' + string.hollerith(bgn) + '_' + string.hollerith(end)
        
        return sp
    #
    
    def parse_basic_file_name_info(self, flnm):
        # First check the extensions to vet the input
        xtns = flnm.split('.')
        ext_k = xtns[len(xtns)-1]
        flhd_k = flnm[:len(flnm)-len(ext_k)-1]
        self.ext  += [ext_k]
        self.flhd += [flhd_k]
        
        
        saved_datatype = self.datatype
        
        if ext_k == "txt": 
            # the extension should be txt
            self.datatype = "type1"
        elif ext_k == "bed":
            # the extension must be txt or bed
            self.datatype = "type2"
        else:
            # shouldn't actually have to call this, but it is here
            # anyway. Never know what bright ideas we have when it is
            # way past midnight.
            print "ERROR(DATA): unrecognized file type '%s'" % flnm
            print "accepted file types have extension 'txt' or 'bed'"
            sys.exit(1)
        #
        
        # There may be a way to combine these _very_ different data
        # types, but it is not easy, so I have to insist on only one
        # type. The program will read multiple files of the _same_
        # type of data (though without any critical assessment of the
        # source presently). Anyway, multiple bed _XOR_ txt files can
        # be read in principle and were read for the txt files.
        if len(self.ext) > 1:
            if not self.datatype == saved_datatype:
                print "ERROR: program does not accept mixed file types."
                print "       Files must have an extension 'txt' or 'bed'"
                print "       but not both!"
                print "       current requested extensions: ", self.ext
                sys.exit(1)
            #
        #
        
        # Next, we make sure the data types match for all the input
        # files. Perhaps in principle, I could mix Teresa's data with
        # Przemek's, but presently, I am not sure how this would work.
        # Therefore, I have forced the situation where only one
        # datatype is allowed.
        
        annotation     = False
        interaction    = ''
        flag_fail      = True
        err_msg        = ''
        
        for xtnsnk in xtns:
            if xtnsnk == "annotated":
                # file must have a tag "annotated"
                annotation = True
            elif xtnsnk == "txt" or xtnsnk == "bed":
                continue
            else:
                interaction += (xtnsnk + " ")
            #
        #
        
        
        
        # Now example specific fatures of each file
        if self.datatype == 'type2':
            # bed type file
            if annotation:
                self.annotation += [True]
            else:
                err_msg += "       + %s does not contain annotated data.\n" % flnm
            #
                
            # file must have a tag indicated the type of
            # interaction; e.g., CTCF, CCD, RNAPII, etc.
            if len(interaction) > 0 and annotation:
                self.interaction += [interaction]
                flag_fail = False
            else:
                err_msg += "       + %s does not have sufficient information.\n" % flnm
            #
            if flag_fail: 
                print "ERROR: unrecognized file type in '%s'" % flnm
                print "cause: ----"
                print err_msg
                print "       'bed' files should look as the following:"
                print ""
                print "       examples:"
                print "       > loops.CTCF_RNAPII.annotated.bed"
                print "       > structure.loops.GSE6352.annotated.bed"
                print ""
                print "       * The 'annotated' label is REQUIRED and should probably"
                print "         be placed just before the file extension type, but that"
                print "         is not essential."
                print ""
                print "       * The 'CTCF_RNAPII' or 'GSE6352' tag indicates the source"
                print "         of the data and should probably be there for information,"
                print "         though again, it is not required."
                print ""
                print "       please check your file for accuracy"
                sys.exit(1)
            #
        else:
            self.interaction += ["type1"] # nothing else to call it
        #
        
        
        #
        return True
    #
    
    def readfile(self, flnm):
        
        # There is a fair amount of information in the ChIA-PET data
        # files, so parse_basic_file_name_info() is used to obtain
        # some other clues about the file before the calcuations
        # begin.
        self.parse_basic_file_name_info(flnm)
        
        
        
        # open the file
        try:
            fp = open(flnm, 'r')
        except IOError:
            print "ERROR: cannot open '%s'" % flnm
            usage()
            sys.exit(1)
        #
        lfp = fp.readlines()
        fp.close()
        
        
        if self.datatype == "type1":
            self.read_ext_txt(lfp, flnm)
        elif self.datatype == "type2":
            self.read_ext_bed(lfp, flnm)
        else:
            # shouldn't be here!!!
            print "ERROR(DATA): problems with readfile(%s)" % flnm
            sys.exit(1)
        #
        return 0
        
    #
    
    def read_ext_bed(self, lfp, flnm):
        debug_read_ext_bed = False
        if debug_read_ext_bed:
            print "read_ext_bed()"
        #
        
        # I think this is getting too crazy, so I have to do some file
        # specification at the very start, so I have introduced two
        # indicators that can help me adapt to all the different file
        # formats.
        
        
        # the first line must contain the version number (if v1) or
        # just the data (if v2).
        s = lfp[0].strip().split()
        if s[0] == 'version:':
            # read the version number
            self.version = s[1]
            self.read_ext_bed_data_v1(lfp, flnm)
        else:
            self.read_ext_bed_data_v0(lfp, flnm)
        #
    #
    
    def read_bed_taglist(self, lfp, flnm):
        debug_read_ext_bed = False
        if debug_read_ext_bed:
            print "read_bed_taglist()"
        #
        
        #01 version: 1
        #02 # CTCF analysis from ChIA-PET data and epigenetic info
        #03 # active
        #04 # open
        #05 # compartment A
        #06 # compartment B
        #07 
        #08 
        #09 chromatin_tags:
        #10 chr      bgn      end    lCTCF rCTCF  PETcnt  cmplx1 cmplx2   cmplx3  open             active           A       B
        #11 chromatin_data:
        #12 chr15 20558294 22455914  R      R       9      1       2       2      0.00473698633025 0.0317128824528 1.0     0.0
        #13 chr15 21925132 22455914  R      R       17     0       1       0      0.0136930039074  0.0404572875493 1.0     0.0
        #14 chr15 23887615 25747260  R      L       7      0       7       9      0.00475413318133 0.164474402372  1.0     0.0
        #15 chr15 25752736 26179721  R      L       17     2       -2      2      0.0429710645573  0.061978758036  1.0     0.0
        #16 chr15 25939090 26179721  R      L       7      1       -2      1      0.0750069608654  0.0585419168769 1.0     0.0
        
        # line01 tells the version index (so the program can adjust to
        # the new formats people will create)
        
        # line02 to line06 A brief description. Not absolutely
        # required, but it becomes useful when you have hundreds of
        # files and you want to figure out what you were doing 6
        # months ago. Any line that starts with '#' is automatically
        # assumed to be a comment.
        
        # line07 to line08 this end of description. The empty lines or
        # lines with some only spaces are ignored.
        
        # line09 specifies that you are reading in the names of the
        # fields in the chromatin file.
        
        # line10 The line that follows contains the list of all tags
        # that will be used. The number of tags must match the number
        # of columns in the data file. Moreover, the tags must adhere
        # to particular names.
        
        # line11 specifies the the data should follow. anything that
        # follows this line should either be the data or a comment
        # "#". Anything else will probably lead to an error.
        
        
        #01 version: 1
        #02 # CTCF analysis from ChIA-PET data and epigenetic info
        #03 # active
        #04 # open
        #05 # compartment A
        #06 # compartment B
        #07 
        #08 
        #09 chromatin_tags:
        #10 chr      bgn      end    lCTCF rCTCF  PETcnt  cmplx1 cmplx2   cmplx3  open             active           A       B
        #11 chromatin_data:
        
        k = 1
        taglist = []
        flag_continue = True
        description = []
        while flag_continue:
            s = lfp[k].strip().split()
            # print s
            if len(s) == 0:
                k += 1
                continue
            #
            if s[0][0] == '#': # first element on the line is a comment
                description += [lfp[k].strip("^#").strip()]
                k += 1
                continue
            #
            if s[0] == "chromatin_tags:":
                k += 1
                taglist = lfp[k].strip().split()
                # print "taglist: ", taglist
                k += 1
            #
            if s[0] == "chromatin_data:":
                k += 1
                break
            #
            if s[0][:3] == "chr" and not s[0] == "chromatin_tags:":
                print "ERROR: reading chromosome data before finished reading the tags!"
                print "       Something is wrong with the file %s" % flnm
                sys.exit(1)
            #
            
        #
        
        return k, taglist, description
    #
    
    # read file formats that Przemek came up with. These were later
    # versions.
    def read_ext_bed_data_v1(self, lfp, flnm):
        
        debug_read_ext_bed = False
        if debug_read_ext_bed:
            print "read_ext_bed_data_v1()"
        #
        
        print "scanning %s...." % flnm
        
        dtbgn, taglist, description = self.read_bed_taglist(lfp, flnm)
        self.taglist += [taglist]
        self.description += [description]
        
        # print dtbgn
        # print self.taglist
        # print self.description
        
        # this obtains the tags below
        
        #09 chromatin_tags:
        #10 chr      bgn      end    lCTCF rCTCF  PETcnt  cmplx1 cmplx2   cmplx3  open             active           A       B
        #11 chromatin_data:
        #12 chr15 20558294 22455914  R      R       9      1       2       2      0.00473698633025 0.0317128824528 1.0     0.0
        #13 chr15 21925132 22455914  R      R       17     0       1       0      0.0136930039074  0.0404572875493 1.0     0.0
        #14 chr15 23887615 25747260  R      L       7      0       7       9      0.00475413318133 0.164474402372  1.0     0.0
        #15 chr15 25752736 26179721  R      L       17     2       -2      2      0.0429710645573  0.061978758036  1.0     0.0
        #16 chr15 25939090 26179721  R      L       7      1       -2      1      0.0750069608654  0.0585419168769 1.0     0.0
        
        kk = 0
        
        #  
        name = ''    # name of the chromosome
        bgn = -1     # beginning locus on the chromosome
        end = -1     # ending locus on the chromosome
        ctcf1 = ''   # right side CTCF
        ctcf2 = ''   # right side CTCF
        nPET = 0     # number of PET counts
        cmplx = []   # complexity
        popen  = 0.0 # open
        pactv  = 0.0 # active
        cmprtA = 0.0 # compartment A
        cmprtB = 0.0 # compartment B
        #
        
        
        for k in range(dtbgn, len(lfp)):
            s = lfp[k].strip().split()
            # print s
            if s[0][0] == '#': # first element on the line is a comment
                continue
            kky = self.makekey(s[0:3])
            
            cmplx = []
            for tlk in range(0, len(taglist)):
                try:
                    if taglist[tlk] == 'chr':
                        name  = s[tlk]
                    elif taglist[tlk] == 'bgn':
                        bgn = int(s[tlk])
                    elif taglist[tlk] == 'end':
                        end = int(s[tlk])
                    elif taglist[tlk] == 'lCTCF':
                        ctcf1 = s[tlk]
                    elif taglist[tlk] == 'rCTCF':
                        ctcf2 = s[tlk]
                    elif taglist[tlk] == 'PETcnt':
                        nPET = int(s[tlk])
                    elif taglist[tlk] == 'cmplx1':
                        cmplx += [int(s[tlk])]
                    elif taglist[tlk] == 'cmplx2':
                        cmplx += [int(s[tlk])]
                    elif taglist[tlk] == 'cmplx3':
                        cmplx += [int(s[tlk])]
                    elif taglist[tlk] == 'open':
                        popen = float(s[tlk])
                    elif taglist[tlk] == 'active':
                        pactv = float(s[tlk])
                    elif taglist[tlk] == 'A':
                        cmprtA = float(s[tlk])
                    elif taglist[tlk] == 'B':
                        cmprtB = float(s[tlk])
                    else:
                        print "ERROR: undefined tag. "
                        print "       allowed tags are the following:"
                        for aa in allowed_tags:
                            print aa
                        sys.exit(1)
                        
                except (ValueError):
                    print "ERROR: incompatable data found for entry"
                    print s
                    sys.exit(1)
                #
            #
            
            
            # create a new dictionary entry if the key is new
            if not self.cdata.has_key(kky):
                self.cdata.update({ kky : Chromatin(name, bgn, end, ctcf1, ctcf2, nPET, cmplx) })
                
                if not self.chrmsm_grp.has_key(name):
                    self.chrmsm_grp.update( { name: Chromosome(name, bgn, end) } )
                else:
                    self.chrmsm_grp[name].add_chrsegment(bgn, end)
            #
            
            
            self.cdata[kky].add_segments(bgn, end, "-")
            self.cdata[kky].state["active"] = pactv
            self.cdata[kky].state["open"]   = popen
            self.cdata[kky].state["A"]      = cmprtA
            self.cdata[kky].state["B"]      = cmprtB
            kk += 1
            if debug_read_ext_bed:
                # tests a small section of the input file
                if kk > 10:
                    self.display_Data(True)
                    sys.exit(0)
                #
            #
        #
        if debug_read_ext_bed:
            self.display_Data(False)
            print "finished read_ext_bed"
            sys.exit(0)
        #
    #
    
    
    
    # read file formats that Przemek came up with. These were later
    # versions.
    def read_ext_bed_data_v0(self, lfp, flnm):
        # ##############################
        # prior versions
        # ##############################
        
        debug_read_ext_bed = False
        if debug_read_ext_bed:
            print "read_ext_bed_data_v0()"
        #
        
        print "scanning %s...." % flnm
        
        # version 0 formats
        # file format1 (loops.CTCF.annotated.bed):
        #
        #                       CTCF    CTCF    PET   cmplx   cmplx   cmplx     
        # chr    bgn     end     1       2      cnt    ity      a       b       open            active
        # chr1	838908	912011	R	L	11	0	4	6	0.0686291944243	0.426918183932
        # chr1	838908	922335	R	R	5	1	6	8	0.0601364066789	0.397329401752
        # chr1	838908	1000167	R	L	7	3	15	17	0.0910398799447	0.5333717808
        # chr1	918286	969271	R	R	5	0	4	4	0.119015396685	0.543963910954
        # chr1	918286	1000167	R	L	75	1	8	8	0.11802493863	0.65697780926
        # chr1	918286	1059032	R	R	4	2	12	12	0.09421226891	0.648785755901
        
        # file format2 (structure.loops.GSE6352.annotated.bed):
        # chr1	1050000	1190000	241	-1	9	9	0.0529785714286	0.418914285714
        # chr1	1585000	1650000	80	-2	0	0	0.0553384615385	0.737907692308
        # chr1	1710000	1840000	154	-2	22	22	0.0699230769231	0.874938461538
        
        # file format2 (structure.TAD_Rao.annotated.bed):
        # chr1	915000	1005000	0.72251	1	10	10	0.116533333333	0.643166666667
        # chr1	1030000	1235000	1.1954	-1	12	12	0.0460487804878	0.507936585366
        # chr1	1255000	1450000	0.9312	0	27	27	0.0741435897436	0.653230769231
        
        # There are other formats, but it appears that if the tag
        # "annotated" is not present, we can discriminate between
        # them.
        
        
        if len(self.ext) > 1:
            print "SORRY: The file '%s' is version 0 type 'bed' file." % flnm
            print ""
            print "        Unfortunately, these version 0 files cannot be input "
            print "        as a group. You can only read them one at a time."
            sys.exit(1)
        #
        nn = len(lfp[0].strip().split())
        if nn == 9:
            self.taglist = [['chr', 'bgn', 'end', 'lCTCF', 'rCTCF', 'PETcnt', 'cmplx1', 'open', 'active']]
        elif nn == 11:
            self.taglist = [['chr', 'bgn', 'end', 'lCTCF', 'rCTCF', 'PETcnt', 'cmplx1', 'cmplx2', 'cmplx3', 'open', 'active']]
        #
        
            
        kk = 0
        for lfpk in lfp:
            s = lfpk.strip().split()
            # print s
            if s[0][0] == '#': # first element on the line is a comment
                continue
            kky = self.makekey(s[0:3])
            
            
            name  = s[0]
            bgn   = int(s[1])
            end   = int(s[2])
            ctcf1 = 'N'
            ctcf2 = 'N'
            nPET  = 0
            cmplx = []
            popen  = 0.0  # "open"
            pactv  = 0.0  # "active"
            if len(s) == 11:
                ctcf1 = s[3]
                ctcf2 = s[4]
                nPET  = int(s[5])
                cmplx = [int(s[6]), int(s[7]), int(s[8])]
                popen  = float(s[9])
                pactv  = float(s[10])
            elif len(s) == 9:
                # not sure why some files have floats, others have
                # integer. So what _is_ this variable?
                if string.isFloat(s[3]):
                    nPET  = float(s[3])
                    # print "set float"
                else:
                    nPET  = int(s[3])
                    # print "set int"
                #
                cmplx = [int(s[4]), int(s[5]), int(s[6])]
                popen  = float(s[7])
                pactv  = float(s[8])
            else:
                print "ERROR: unrecognized structure of 'bed' file '%s'" % flnm
                sys.exit(1)
            #
            #
            
            
            
            if not self.cdata.has_key(kky):
                self.cdata.update({ kky : Chromatin(name, bgn, end, ctcf1, ctcf2, nPET, cmplx) })
                
                if not self.chrmsm_grp.has_key(name):
                    self.chrmsm_grp.update( { name: Chromosome(name, bgn, end) } )
                else:
                    self.chrmsm_grp[name].add_chrsegment(bgn, end)
            #
            
            
            self.cdata[kky].add_segments(bgn, end, "-")
            self.cdata[kky].state["active"] = pactv
            self.cdata[kky].state["open"]   = popen
            kk += 1
            if debug_read_ext_bed:
                # tests a small section of the input file
                if kk > 10:
                    self.display_Data(True)
                    sys.exit(0)
                #
            #
        #
        if debug_read_ext_bed:
            self.display_Data(False)
            print "finished read_ext_bed"
            sys.exit(0)
        #
    #
    
    
    
    # read file formats that Teresa set up. This was the initial set
    # of data.
    def read_ext_txt(self, lfp, flnm):
        # file format:
        
        self.taglist += [['chr', 'bgn', 'end', 'lCTCF', 'open', 'active', 'repressed', 'A', 'B']]
        # used as a identifier                 |used in the program ----------------------------------                           
        # region of chromosome                 |specific location in region             state   length
        # chr   bgn             end            |chr     bgn             end              
        # chr1	10436116	10556545	chr1	10437400	10438000	active	600
        # chr1	10436116	10556545	chr1	10438000	10438466	active	466
        # chr1	10436116	10556545	chr1	10438558	10440200	active	1642
        # 0     1               2               3       4               5               6       7
        
        
        # 160617wkd: I must admit that building a hash in this way
        # made some very simple code that is far more flexible
        # compared to the original where data was a vector of
        # objects. In this version, I can read in more than one file
        # and update according to my needs.
        for lfpk in lfp:
            s = lfpk.strip().split()
            if s[0][0] == '#': # first element on the line is a comment
                continue
            kky = self.makekey(s[0:3])
            name = s[0]
            
            if not self.cdata.has_key(kky):
                c_bgn = int(s[1])
                c_end = int(s[2])
                self.cdata.update({ kky : Chromatin(s[0], c_bgn, c_end) })
                if not self.chrmsm_grp.has_key(name):
                    self.chrmsm_grp.update( { name: Chromosome(name, c_bgn, c_end) } )
                else:
                    self.chrmsm_grp[name].add_chrsegment(c_bgn, c_end)
            #
            
            
            l_bgn   = int(s[4])
            l_end   = int(s[5])
            l_state = s[6]
            self.cdata[kky].add_segments(l_bgn, l_end, l_state)
        #
        # self.display_Data(False)
        # sys.exit(0)
    #
    
    
    
    
    
    def display_Data(self, show_all):
        # print "keys: ", self.cdata.keys()
        keylist = self.cdata.keys()
        if len(self.keylist) > 0:
            print " using ordered keylist"
            keylist = self.keylist
        #
            
        for ky in keylist:
            
            # this now shows that I can extract whatever datapoint I
            # desire from the set of experimental data. Therefore, I
            # can read in all the data and extract it according to my
            # needs.
            
            if not show_all:
                print "total number of segments: ", len(self.cdata[ky].segments)
            #
            for seg in self.cdata[ky].segments:
                if show_all:
                    print "%4d  %s | %s  %s" \
                        % (len(self.cdata[ky].segments), \
                           self.cdata[ky].disp_head(), \
                           seg.disp(), self.cdata[ky].state)
                else:
                    print "%4d  %s | %s" \
                        % (len(self.cdata[ky].segments), \
                           self.cdata[ky].disp_head(), \
                           seg.disp())
                    break # just show the first segment
                #
            #
        #
        print "length of the dataset: ", len(self.cdata)
    #
    
    def ordered_keylist(self):
        print "length of the dataset: ", len(self.cdata)
        self.keylist = self.ins_sort_keys(self.cdata.keys())
        return self.keylist
    #
    
    def get_type1_prob(self, ky):
        # files with extension 'txt'
        s_active    = []
        s_open      = []
        s_repressed = []
        s_A         = []
        s_B         = []
        
        chr_len = float(self.cdata[ky].end - self.cdata[ky].bgn)
        
        for seg in self.cdata[ky].segments:
            #
            if seg.state == 'open':
                seg = self.check_range(ky, seg)
                s_open += [seg]
                if self.debug_Data:
                    print "%4d  %s | %s" %(len(self.cdata[ky].segments), self.cdata[ky].disp_head(), seg.disp())
                #
            #
            
            if seg.state == 'active':
                seg = self.check_range(ky, seg)
                s_active += [seg]
                if self.debug_Data:
                    print "%4d  %s | %s" %(len(self.cdata[ky].segments), self.cdata[ky].disp_head(), seg.disp())
                #
            #
            if seg.state == 'repressed':
                seg = self.check_range(ky, seg)
                s_repressed += [seg]
                if self.debug_Data:
                    print "%4d  %s | %s" %(len(self.cdata[ky].segments), self.cdata[ky].disp_head(), seg.disp())
                #
            #
            if seg.state == 'A':
                seg = self.check_range(ky, seg)
                s_A += [seg]
                if self.debug_Data:
                    print "%4d  %s | %s" %(len(self.cdata[ky].segments), self.cdata[ky].disp_head(), seg.disp())
                #
            #
            if seg.state == 'B':
                seg = self.check_range(ky, seg)
                s_B += [seg]
                if self.debug_Data:
                    print "%4d  %s | %s" %(len(self.cdata[ky].segments), self.cdata[ky].disp_head(), seg.disp())
                #
            #
        #
        
        
        p_active    = self.process_segs(s_active, chr_len)
        if p_active > 1.0:
            print "%s     %10d" % (self.cdata[ky].disp_head(), int(chr_len))
            sys.exit(1)
        #
        
        p_open      = self.process_segs(s_open, chr_len)
        if p_open > 1.0:
            print "%s     %10d" % (self.cdata[ky].disp_head(), int(chr_len))
            sys.exit(1)
        #
        
        p_repressed = self.process_segs(s_repressed, chr_len)
        if p_repressed > 1.0:
            print "%s     %10d" % (self.cdata[ky].disp_head(), int(chr_len))
            sys.exit(1)
        #
        
        p_A         = self.process_segs(s_A, chr_len)
        if p_A > 1.0:
            print "%s     %10d" % (self.cdata[ky].disp_head(), int(chr_len))
            sys.exit(1)
        #
        
        p_B         = self.process_segs(s_B, chr_len)
        if p_B > 1.0:
            print "%s     %10d" % (self.cdata[ky].disp_head(), int(chr_len))
            sys.exit(1)
        #
        
        return p_active, p_open, p_repressed, p_A, p_B
    #
    
    def get_type2_prob(self, ky):
        # files with extension 'bed'
        p_active    = self.cdata[ky].state["active"]    
        p_open      = self.cdata[ky].state["open"]      
        p_repressed = self.cdata[ky].state["repressed"] 
        p_A         = self.cdata[ky].state["A"]         
        p_B         = self.cdata[ky].state["B"]         
        
        return p_active, p_open, p_repressed, p_A, p_B
    #
    
    
    def check_range(self, ky, seg):
        flag_clash11 = False
        flag_clash12 = False
        flag_clash21 = False
        flag_clash22 = False
        chrxb = self.cdata[ky].bgn
        chrxe = self.cdata[ky].end
        
        #
        sbgn = seg.bgn; sbgno = seg.bgn
        send = seg.end; sendo = seg.end
        if seg.bgn < chrxb:
            flag_clash11 = True
            seg.bgn = chrxb
        #
        
        if seg.bgn > chrxe:
            flag_clash21 = True
            seg.bgn = chrxe
        #
        
        if seg.end > chrxe:
            flag_clash12 = True
            seg.end = chrxe
        #
        if seg.end < chrxb:
            flag_clash22 = True
        #
        
        chr_len = self.cdata[ky].end - self.cdata[ky].bgn
        if flag_clash11 or flag_clash12:
            print "WARNING! %s   range(%10d)  seg: %10d  %10d len(%10d)  state(%s)" \
                % (self.cdata[ky].disp_head(), chr_len, seg.bgn, seg.end, seg.end - seg.bgn, seg.state), sbgn, send
            s = ''
            if flag_clash11:
                s += "clash at beginning: " 
            if flag_clash12:
                s += "clash at ending:    "
            files = ''
            for fk in range(0, len(self.flhd)):
                files += "%s.%s " % (self.flhd[fk], self.ext[fk])
            s += "segment(%10d to %10d) vs chr(%10d to %10d), in files [ %s]" % (sbgno, sendo, chrxb, chrxe, files)
            print s
        #
        if flag_clash21 or flag_clash22:
            print "ERROR!   %s   range(%10d)  seg: %10d  %10d len(%10d)  state(%s)" \
                % (self.cdata[ky].disp_head(), chr_len, seg.bgn, seg.end, seg.end - seg.bgn, seg.state), sbgn, send
            print "incompatibility of data"
            s = ''
            if flag_clash21:
                s += "segment begins after the end:     "
            if flag_clash22:
                s += "segment ends after the beginning: "
            files = ''
            for fk in range(0, len(self.flhd)):
                files += "%s.%s " % (self.flhd[fk], self.ext[fk])
            s += "segment(%10d to %10d) vs chr(%10d to %10d), in file [ %s]" % (sbgno, sendo, chrxb, chrxe, files)
            print s
            sys.exit(1)
        #
        return seg
    #
    
    
    
    
    def ins_sort_segs(self, seg):
        for i in range(1,len(seg)):    
            j = i                    
            while j > 0 and seg[j].bgn < seg[j-1].bgn: 
                seg[j], seg[j-1] = seg[j-1], seg[j] 
                j=j-1 
        return seg
    #
    
    
    
    
    
    def ins_sort_keys(self, kky):
        for i in range(1,len(kky)):    
            j = i                    
            while j > 0 and kky[j] < kky[j-1]: 
                kky[j], kky[j-1] = kky[j-1], kky[j] 
                j=j-1 
        return kky
    #
    
    
    
    
    
    def process_segs(self, segs, chr_len):
        segs = self.ins_sort_segs(segs)
        if self.debug_Data:
            for ss in segs:
                print "%10d  %10d" % (ss.bgn, ss.end)
            # sys.exit(0)
        #
        
        # make sure that there are no overlaps between the segments
        # described in the experimental data
        
        for j in range(0, len(segs)-1):
            if segs[j+1].bgn < segs[j].end:
                print "WARNING resetting: (%10d, [%10d)|(%10d], %10d)"\
                    % (segs[j].bgn, segs[j].end, segs[j+1].bgn, segs[j+1].end)
                segs[j+1].bgn = segs[j].end + 1
                sys.exit(0) # presently, stop if there are problems!!!
            #
        #
        
        # calculate the fractional percentage of overlap 
        n = 0
        for j in range(0, len(segs)):
            n += segs[j].end - segs[j].bgn
        p = float(n)/chr_len
        
        return p 
    #
#    

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
    
    
#

if __name__ == '__main__':
    # running the program
    main(sys.argv)

