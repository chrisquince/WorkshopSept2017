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

    parser.add_argument("blast_file", help="m8 blast file format")
    
    parser.add_argument("ko_length_file", help="csv file of mean ko gene lengths")

    parser.add_argument("gene_list_file", help="gene to ko mappings")

    args = parser.parse_args()

    import ipdb; ipdb.set_trace()
    
    #hardset reads lengths at 150 nt or 50aa
    readLength = 50
    minEValue = 1.0e-10
    
    ko_lengths = {}
    with open(args.ko_length_file, 'r') as f:
        reader = csv.reader(f, dialect='excel', delimiter=',')
        for row in reader:
            (koId, nGenes, meanLength, medLength, sdLength) = row
            ko_lengths[koId] = meanLength
    
    geneKOMap = {}
    with open(args.gene_list_file, 'r') as f:
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        for row in reader:
            (koId, geneId) = row
            geneKOMap[geneId] = koId
    
    koCov = defaultdict(float)
    with open(args.blast_file, 'r') as f:
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        currId = "None"
        hits = []
        
        for row in reader:
            (queryId, subjectId, percIdentity, alnLength, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore) = row
            if currId == "None":
                currId = queryId
            if currId != queryId:
                #process hits
                nHits = len(hits)
                if nHits > 0:
                    fHit = readLength/float(nHits)
                    for geneId in hits:
                        if geneId in geneKOMap:
                            koCov[geneKOMap[geneId]] += fHit
                hits = []
            
            currId = queryId
            if float(eVal) < minEValue:
                hits.append(subjectId)
            
        #process final hits
        nHits = len(hits)
        if nHits > 0:
            fHit = readLength/float(nHits)
            for geneId in hits:
                if geneId in geneKOMap:
			koCov[geneKOMap[geneId]] += fHit
                

    for koId, sumMap in koCov.iteritems():
        koNCov = sumMap/float(ko_lengths[koId])
        print koId + "," + str(sumMap) + "," + str(koNCov)  
    #import ipdb; ipdb.set_trace()
   

if __name__ == "__main__":
    main(sys.argv[1:])
