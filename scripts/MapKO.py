from Bio import SeqIO
import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
import argparse
import math
import csv
from collections import defaultdict

def main(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument("module_map_file", help="mapping of kos to modules")
    
    parser.add_argument("ko_coverage_file", help="csv file of ko gene coverages")

    args = parser.parse_args()

#    import ipdb; ipdb.set_trace()
    
    koModuleMap = defaultdict(list)
   
    with open(args.module_map_file, 'r') as f:
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        for row in reader:
            (koId, moduleId) = row
            koModuleMap[koId].append(moduleId)
    
    moduleCov = defaultdict(float)
    with open(args.ko_coverage_file, 'r') as f:
        reader = csv.reader(f, dialect='excel', delimiter=',')
        
        for row in reader:
            (geneId, totalCov, normCov) = row
            totalCov = float(totalCov)
            NHits = float(len(koModuleMap[geneId]))
            for module in koModuleMap[geneId]:
                moduleCov[module] += totalCov/NHits
                

    for module, sumMap in moduleCov.iteritems():
        print module + "," + str(sumMap) 
    
   

if __name__ == "__main__":
    main(sys.argv[1:])
