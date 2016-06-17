#!/usr/bin/env python

""" Initial implementation Diego Alvarado, heavily modified by Isaac Overcast """

import os, copy, numpy as np
import argparse

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument("-f", dest="VCFname", required=True, help="name of the VCF input file being converted")
parser.add_argument("-p", dest="populations", required=True, help="Input file containing population assignments per individual")
parser.add_argument("-o", dest="OUT", default='output', help="name for dadi output file")
parser.add_argument("--GQ", dest="GQual", help="minimum genotype quality tolerated", default=20)
parser.add_argument("-v", dest="verbose", help="Set verbosity. Dump tons of info to the screen", default=False)

options = parser.parse_args()

infilename = options.VCFname
outfilename = '%s.dadi' %(options.OUT)
minqual = int(options.GQual)
missfilename = '%s.miss' %(options.OUT)

# Read in the vcf file and grab the line with the individual names
# Here we populate the 'indnames' array from the vcf file.
# Would it be preferrable to just read in from the pop assignment
# file and only use individuals that are specified there? Could give
# us more flexibility.
# Add the 'U' to handle opening files in universal mode, squashes the
# windows/mac/linux newline issue.
with open(infilename, 'rU') as infile:
    for line in infile:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                row = line.strip().split()
                # VCF file format spec says that if the file contains genotype
                # data then "FORMAT" will be the last column header before
                # sample IDs start
                startcol = row.index('FORMAT')
                # This is trying to strip something off individual names 
                # that isn't there in my data. Can't hurt to leave it but it 
                # shouldn't be necessary. We only need indnames to be iterated here.
                # 3/9/15 iao
                indnames = [x.replace('.sorted.bam', '') for x in row[startcol+1:]]
            else:
                pass
        else:
            break

# Here we need to read in the individual population
# assignments file and do this:
# - populate the locs dictionary with each incoming population name
# - populate another dictionary with individuals assigned to populations
# Add the 'U' to handle opening files in universal mode, squashes the
# windows/mac/linux newline issue.
with open(options.populations, 'rU') as popsfile:
    pop_assignments = {}
    locs = {}

    # Skip the header line
    next( popsfile )
    for line in popsfile:
        tmp = line.split()
        pop_assignments[tmp[0]] = tmp[1]
        if(options.verbose): print tmp[0] , pop_assignments[tmp[0]]

        # Populate the 'locs' dict with just the location names
        # This is a not an ideal way of doing this but it'll work
        locs.setdefault( tmp[1], {} )

    print "Specifiying", len( locs ), "populations:", locs.keys() 

header = ['Gene', 'Pos', 'Ref', 'Alt']
for x in indnames:
    header.append(x)
geno = np.empty([1,len(header)])
geno = np.vstack((geno, header))
geno = np.delete(geno,0,0)
with open(missfilename, 'w') as missfile:
    print >> missfile, 'Tag\tP0/0\tP0/1\tP1/1\tM0/0\tM0/1\tM1/1'
    with open(infilename, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                pass
            else:
                countP = {'00':0, '01':0, '11':0}    #stores the number of passing homozygous reference, heterozygous, and homozygous alternative genotypes
                countF = {'00':0, '01':0, '11':0}    #stores the number of filtered out homozygous reference, heterozygous, and homozygous alternative genotypes
                row = line.strip().split()
                if(options.verbose): print( row )
                genorow = row[:2]
                genorow.append(row[3])
                genorow.append(row[4])
                if len(row[4])==1:    #this removes on-biallelic SNPs
                    dic = {'0':row[3],'1':row[4]}
                    for IND in row[startcol+1:]:
                        indinfo = IND.split(':')
                        #Print out each individuals genotype
                        #if(options.verbose): print( IND )
                        # Add the call to replace pipe with forward slash
                        # to handle the situation where vcf uses the pipe to
                        # separate alleles.
                        gen = indinfo[0].replace('|','/').split('/')
                        # Gqual is getting hard coded for now because my vcf
                        # is prefiltered, doesn't have qual scores. Could be
                        # useful to fix this some day. Here's the old way:
                        # int(indinfo[2])
                        Gqual = -1                     
                        if (Gqual >= minqual) or (Gqual == -1):
                            indgeno = '%s/%s' %(dic.get(gen[0]), dic.get(gen[1]))
                            if indinfo[0] == '0/0':
                                countP['00'] += 1
                            elif indinfo[0] == '0/1':
                                countP['01'] += 1
                            elif indinfo[0] == '1/1':
                                countP['11'] += 1
                        else:
                            if indinfo[0] == '0/0':
                                countF['00'] += 1
                            elif indinfo[0] == '0/1':
                                countF['01'] += 1
                            elif indinfo[0] == '1/1':
                                countF['11'] += 1
                            indgeno = '-/-'
                        genorow.append(indgeno)
                    geno = np.vstack((geno, genorow))
                # Print out each modified allele row
                if(options.verbose): print( genorow )
                print >> missfile, '%s\t%d\t%d\t%d\t%d\t%d\t%d' %(row[0], countP.get('00'), countP.get('01'), countP.get('11'), countF.get('00'), countF.get('01'), countF.get('11'))
with open(outfilename,'w') as dadifile:
    print( "************************" )
    print >> dadifile, 'InOne\tInTwo\tAllele1\t%s\tAllele2\t%s\tGene\tPosition' %('\t'.join(locs.keys()), '\t'.join(locs.keys()))
    for row in xrange(1,np.shape(geno)[0]):
        # Create a deep copy of the empty array of each population for tracking
        # allele counts
        locsRow = copy.deepcopy(locs)
        if(options.verbose): print( locsRow )
        alleRef = geno[row,2]
        alleDev = geno[row,3]
        for col in xrange(4,np.shape(geno)[1]):
            # Get the population for this individual
            # Accessing the dict with .get() we can pass
            # in a default in case the key doesn't exist
            # Then we test for this default key and if we find it
            # this means the individual doesn't have a population assignment
            # so we pass on including them in the output file.
            # <TODO> Would be nice to keep track of all the individuals
            # we _don't_ use and report this, perhaps as part of the miss file.
            loc = pop_assignments.get( geno[0,col], {} )
            if( not loc ):
                # print >> missfile, 'Skipping %s', geno[0,col]
                print "Skipping", geno[0,col]
                continue
            if(options.verbose): print geno[row,0]
            locsRow[loc].setdefault(alleRef,[])
            locsRow[loc].setdefault(alleDev,[])
            if( geno[row,col] != '-/-' and geno[row,col] != './.' and geno[row,col] != 'None/None' ):
                for a in geno[row,col].split('/'):
                    if(options.verbose): print( loc + " " + a )
                    locsRow[loc][a].append(1)
                    exit
        reffreq = []
        devfreq = []
        for l in locsRow.keys():
            reffreq.append(str(sum(locsRow[l].get(alleRef))))
            if alleDev != alleRef:
## This is making a stupid assumption about something.
                devfreq.append(str(sum(locsRow[l].get(alleDev))))
        print >> dadifile, '-%s-\t-%s-\t%s\t%s\t%s\t%s\t%s\t%s' %(alleRef, alleDev, alleRef, '\t'.join(reffreq), alleDev, '\t'.join(devfreq), geno[row,0], geno[row,1])

#Clean up the .miss file bcz it's not meaningful right now.
os.remove( missfilename )

print '{:*^60s}'.format(" Conversion succesfully finished ")
print
