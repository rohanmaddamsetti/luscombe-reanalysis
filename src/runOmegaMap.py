#!/usr/bin/python
##runOmegaMap.py by Rohan Maddamsetti
##This short script is a wrapper for omegaMap.
##A second script will grep through the result files,
##printing out a text file with two headings: locus_tag, and theta_s.

from subprocess import Popen
from os import listdir
from os.path import splitext
import sys

def getLocusTag(path):
    '''This short helper function pulls out a locus_tag from a path.
    '''
    rest, extension = splitext(path)
    ##pull out the locus_tag from the path.
    try:
        locus_tag_start = path.index('ECB')
        locus_tag_end = path.index(extension)
        locus_tag = path[locus_tag_start:locus_tag_end]
    except ValueError:
        locus_tag = None
    return locus_tag
    
def createConfigFile (fasta_path):
    '''This function takes the path to an alignment file,
    and writes a *.ini config file to a directory in the cwd,
named 'configs'. The * wildcard is replaced by the name of the
fasta file, *.fasta. 

parameter settings for omegaMap (see Luscombe supplement section 2.1.7):
 norders = 10
 niter = 150000
 thinning = 150
 muPrior = improper_inverse
 kappaPrior = improper_inverse
 indelPrior = improper_inverse
 omegaPrior = inverse
 omegaParam = 0.01, 100
 rhoPrior = inverse
 rhoParam = 0.01, 100
 muStart = 0.1
 kappaStart = 3.0
 indelStart = 0.1
 omega_model = variable
 oBlock = 15
 rho_model = variable
 rBlock = 15
    '''

    locus_tag = getLocusTag(fasta_path)
    #print locus_tag

    ##open the init file for writing. The path is hardcoded.
    config_outfile = open('./configs/' + locus_tag + '.ini', 'w')

    config_outfile.write("FASTA = " + fasta_path + "\n")
    config_outfile.write("pi = .016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393" + "\n")
    config_outfile.write("norders = 10" + "\n")
    config_outfile.write("niter = 150000" + "\n")
    config_outfile.write("thinning = 150" + "\n")
    config_outfile.write("muPrior = improper_inverse" + "\n")
    config_outfile.write("kappaPrior = improper_inverse" + "\n")
    config_outfile.write("indelPrior = improper_inverse" + "\n")
    config_outfile.write("omegaPrior = inverse" + "\n")
    config_outfile.write("omegaParam = 0.01, 100" + "\n")
    config_outfile.write("rhoPrior = inverse" + "\n")
    config_outfile.write("rhoParam = 0.01, 100" + "\n")
    config_outfile.write("muStart = 0.1" + "\n")
    config_outfile.write("kappaStart = 3.0" + "\n")
    config_outfile.write("indelStart = 0.1" + "\n")
    config_outfile.write("omega_model = variable" + "\n")
    config_outfile.write("oBlock = 15" + "\n")
    config_outfile.write("rho_model = variable" + "\n")
    config_outfile.write("rBlock = 15" + "\n")

    config_outfile.close()

def createAllConfigs():
    '''This function creates all config files for all input files to
    OmegaMap.'''
    datafiles = [x for x in listdir("./data") if x.find(".fasta")]
    for basename in datafiles:
        path = "./data/" + basename
        createConfigFile(path)
    return None
    
def runOmegaMap(init_path):
    '''This function runs OmegaMap with the given init file.'''
    name = getLocusTag(init_path)
#    print name
    out_filename = "./results/omegamap_results/" + name + ".txt" 
    process_args = ['./omegamap/omegaMap.dat', init_path, '-outfile', out_filename]
    proc =  Popen(process_args, shell=False,
                  stdin=None, stdout=None, stderr=None, close_fds=True)
    return None

def main(argv=None):
    ##test_path = "./data/btrans_ECB_00002.fasta"
    #createConfigFile(test_path)
    createAllConfigs()
    ##runOmegaMap("configs/ECB_00002.ini")
    if argv is None:
        argv = sys.argv
    ##get the seed value from the command line, when run on the cluster.
    seed = argv[-1]
    print seed
    configs = [x for x in listdir("./configs") if getLocusTag(x)]
    print configs
    
if __name__ == "__main__":
    main()
