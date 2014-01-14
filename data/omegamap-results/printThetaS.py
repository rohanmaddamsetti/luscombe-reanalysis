#!/usr/bin/python
##printThetaS.py by Rohan Maddamsetti
##This script writes a csv file containing locus_tag, thetaS
##from the files in results/omegamap_results.
##This output file is called thetaS.csv.
##This also prints out the maximum omega (dN/dS) for the locus,
##to an output file called max_omega.csv.

from os import listdir
import subprocess
import csv 

def makeSummaries():
    '''This function runs the summarize.dat function on all the omegamap
    result files, with a burn-in of 50000. '''
    inputdir = "/Users/Rohan/Desktop/Projects/luscombe_reanalysis/results/omegamap_results/"
    inputfiles = [x for x in listdir(inputdir) if x.startswith('ECB')]
    for filename in inputfiles:
        locus_tag,ext = filename.split('.')
        thecall = ['../omegamap/summarize.dat', '50000',
                         'omegamap_results/'+filename,
                         '>omegamap_summaries/'+locus_tag+'_summary.txt']
        call_string = ' '.join(thecall)
#        print call_string
        subprocess.call(call_string, shell=True)
    return None

def print_thetaS():
    '''This function goes through the summary files, and prints the csv file
    that I want for the next bit of the analysis.'''

    myWriter = csv.writer(open('thetaS.csv', 'w'))
    myWriter.writerow(['locus_tag','thetaS'])

    summarydir = "/Users/Rohan/Desktop/Projects/luscombe_reanalysis/results/omegamap_summaries/"
    summary_files = [x for x in listdir(summarydir) if x.startswith('ECB')]
    for filename in summary_files:
        #print filename
        locus_tag,sep,rest = filename.rpartition('_')
        #print locus_tag
        summary_handle = open('omegamap_summaries/' + filename, 'r')
        line1 = summary_handle.readline()
        line2 = summary_handle.readline()
        line3 = summary_handle.readline()
        if len(line1):
            thetaS = summary_handle.readline().split('\t')[10]
        else:
            thetaS = "NA"
        #print thetaS
        myWriter.writerow([locus_tag,thetaS])
        summary_handle.close()

def print_max_omega():
    '''This function goes through the summary files, and prints the csv file
    that I want for the next bit of the analysis. It gets the maximum omega
    value for the current gene that it is going through.'''

    myWriter = csv.writer(open('max_omega.csv', 'w'))
    myWriter.writerow(['locus_tag','max_omega'])

    summarydir = "/Users/Rohan/Desktop/Projects/luscombe_reanalysis/results/omegamap_summaries/"
    summary_files = [x for x in listdir(summarydir) if x.startswith('ECB')]
    for filename in summary_files:
        #print filename
        locus_tag,sep,rest = filename.rpartition('_')
        #print locus_tag
        summary_handle = open('omegamap_summaries/' + filename, 'r')
        line1 = summary_handle.readline()
        line2 = summary_handle.readline()
        line3 = summary_handle.readline()
        line = summary_handle.readline()
        max_omega = 0
        if len(line1):
            while (line): 
                cur_omega = line.split('\t')[2]
                if cur_omega > max_omega:
                    max_omega = cur_omega
                line = summary_handle.readline()
        else:
            max_omega = "NA"
        print max_omega
        myWriter.writerow([locus_tag,max_omega])
        summary_handle.close()


##makeSummaries()
##print_thetaS()
##print_max_omega()
