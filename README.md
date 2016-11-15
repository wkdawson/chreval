# chreval

# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

chreval.py (CHRomatin EVALuation program for heatmaps)

anal_loops.py

make_heatmap.py

my_generation.py

heatmaps.py (contributed progam by Michal Kadlof and/or Przemek Szalaj

_Chreval_ is the main driver for calculating the free energy of
observed heatmaps.  Chreval is a dynamic programming method for
determining the most probably chromatin structure arrangement as well
as the distribution of chromatin structure arrangements as a function
of free energy of various motifs found in the structure.

_Anal_loops_ is for handling a group of heatmap files with the
properly. The format is Przemek's bed format. The program calls
Chreval. You should make sure that the files listed in the bed file
also exist in your directory.

_Make_heatmap_ is a tool to generate heatmaps either from *.heat files
or the output from Chreval *.clust

_My_generation_ is a program to generate heatmaps. Entries must be
given in dot bracket format and pseudoknots must be added separately
on a different line.


* how to run

To run, the minimum information you need is a heat map data file.
This file contains experimental data from a source (particularly
ChIA-PET but possibly Hi-C). This data (particularly ChIA-PET) is
highly correlated with positions of the CTCF binding proteins and
cohesin. The CTCF sites largely correspond to regions that form loops
and typically regulate the chromatin.

The heat maps consists of a symmetric matrix where the indices (i,j),
corresponding to row and column positions, indicate points on the
chromatin chain where segment i and j interact with each other. For
simplicity, we assume i < j, and look only at the upper triangle. A
reflection of matrix is found across the diagonal. The intensity of a
given point (i,j) is proportional to the frequency that this
particular interaction was encountered in the experimental data, so
small numbers mean only weak interactions, and large numbers mean
strong interactions.

Once you obtain a heat make, or have created one using
my_generation.py, you can run this program with the following command
line example:

> chreval.py -f chr17_55674644_55751596_res5kb.heat  

Here, we call the file chr17_55674644_55751596_res5kb.heat
"chrN_x_y_res5kb.heat" for short. The program requires that
"chrN_x_y_res5kb.heat" consist of a matrix of integer values that
indicate the instances of CTCF and other interactions. The input file
"chrN_x_y_res5kb.heat" must have the extension "heat". The output
consists of a directory with the same name as the input file minus the
extension ".heat"; i.e., the directory name "chrN_x_y_res5kb". This
directory contains (in separate files) the top layer of suboptimal
structures within some specifable energy range from the mimumum free
energy (default is 10 kcal/mol). These files have the extension "DBN"
and can be read by the 3rd party program VARNA. Additionally, Chevral
deposits two additional files: chrN_x_y_res5kb_BDwt.clust contains a
matrix with the Boltzmann probabilities for different interactions and
chrN_x_y_res5kb_summary.txt contains a shorthand list of the secondary
structures.

There are a variety of additional options. Please run

> chreval.py -h 

to obtain additional information on additional command line options.

Version 1.0
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up

The package can run as is, but some applications require installation
of the following python packages

matplotlib and numpy


* Configuration

Standard

* Dependencies

matplot and numpy

* Database configuration

None

* How to run tests

some examples are provided in the directory test

* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

to consult about bugs and other issues of the code, please contact 

Wayne Dawson

Laboratory of Functional and Structural Genomeics

Center of New Technologies

University of Warsaw

Banacha 2C, 02-098 Warsaw

email w.dawson@cent.uw.edu.pl

* Repo owner or admin

Wayne Dawson
* Other community or team contact
