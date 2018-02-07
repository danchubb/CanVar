#! /usr/bin/python

#Adapted and expanded from https://github.com/pontikos/uclex_browser/blob/master/load_vcfs/uclex_to_vcf.py 

from __future__ import print_function
import argparse
import os
import os.path
import sys
import re
#import tabix
import pysam
from pysam import VCF
import csv
import numpy
import cPickle as pickle
import pandas
import pandas as pd
from collections import defaultdict

metrics = ['DP', 'GQ']

def convert_to_int(val):
    """
    Converts string to int if possible, otherwise returns initial string
    """
	
    #print("WHAT?",val)
	
    try:
        return int(val)
    except ValueError:
        #print("DDDDDDDDDDD",val)
        pass
    try:
        return float(val)
    except ValueError:
        return val


def get_histogram_for_variant(distr):
    #hist=numpy.array([1,2])
    hist=numpy.histogram(distr,[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])
    return hist[0]

usage_example = ' '
parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, epilog = usage_example) 
#compulsory arguments
parser.add_argument('--infile', dest='infile', help = "VCF file required", required=False, default=None)
#parser.add_argument('--chrom', dest='chrom', help = "chromosome", required=False, default=None)
parser.add_argument('--phenotypes', dest='phenotypes', help = "file with phenotype data", default=None)
parser.add_argument('--pops', dest='populations', help = "file with population data", default=None)
args = parser.parse_args()

vcf_file=args.infile
popfile=args.pops
#chrom=args.chrom
phenotypes=args.phenotypes


#phenotypes='/data/cancgene-guard2/dchubb/canvar/NSCCG.csv'


#POPS=pandas.read_csv('/data/cancgene-guard2/dchubb/canvar/pop.csv',dtype={'phenotype':str})
POPS=pandas.read_csv(popfile,dtype={'phenotype':str})
POPS.index=POPS.phenotype


samplePheno=dict()
phenoSample=dict()
sampleSex=dict()
SexSample=dict()

for l in csv.DictReader(file(phenotypes,'r'),delimiter=','):
    samplePheno[l['sample']]=l['phenotype']
    phenoSample[l['phenotype']]=phenoSample.get(l['phenotype'],[])+[l['sample']]
    sampleSex[l['sample']]=l['sex']
    SexSample[l['sex']]=SexSample.get(l['sex'],[])+[l['sample']]



H=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']



vcf=pysam.VariantFile(vcf_file,'r')



desc=vcf.header.info['CSQ'].record['Description']
INF=re.sub(r"#CHROM.*","",str(vcf.header))
INF = os.linesep.join([s for s in INF.splitlines() if s])
print(INF)
print('##INFO=<ID=DP_HIST,Number=A,Type=String,Description="Histogram for DP; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">')
print('##INFO=<ID=GQ_HIST,Number=A,Type=String,Description="Histogram for GQ; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5">')

print('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">')
print('##INFO=<ID=AC_2,Number=A,Type=Integer,Description="CRC AF">')
print('##INFO=<ID=AC_3,Number=A,Type=Integer,Description="CRC AF">')
print('##INFO=<ID=AC_Adj,Number=A,Type=Integer,Description="CRC Heterozygous Counts">')
print('##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description="CRC Hemizygous Counts">')
print('##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="CRC AF">')
print('##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="CRC AF">')
print('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">')
print('##INFO=<ID=AF_2,Number=A,Type=Integer,Description="CRC AF">')
print('##INFO=<ID=AF_3,Number=A,Type=Integer,Description="CRC AF">')
print('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">')
print('##INFO=<ID=AN_2,Number=A,Type=Integer,Description="CRC AF">')
print('##INFO=<ID=AN_3,Number=A,Type=Integer,Description="CRC AF">')
print('##INFO=<ID=AN_Adj,Number=A,Type=Integer,Description="CRC HOMO counts">')
print('##INFO=<ID=Hemi_2,Number=A,Type=Integer,Description="CRC Hemizygous Counts">')
print('##INFO=<ID=Hemi_3,Number=A,Type=Integer,Description="CRC Hemizygous Counts">')
print('##INFO=<ID=Het_2,Number=A,Type=Integer,Description="CRC Heterozygous Counts">')
print('##INFO=<ID=Het_3,Number=A,Type=Integer,Description="CRC Heterozygous Counts">')
print('##INFO=<ID=Hom_2,Number=A,Type=Integer,Description="CRC HOMO counts">')
print('##INFO=<ID=Hom_3,Number=A,Type=Integer,Description="CRC HOMO counts">')




#print(vcf.header.replace("#^CHROM.*",''))
CSQ_HEAD=desc.replace('Consequence type as predicted by VEP. Format: ','').replace(' ','').replace('"','').strip().split('|') 

#print(CSQ_HEAD)

#print("XX",CSQ,"XX")
H=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
print('\t'.join(H))
d=0
#vcf_chrom=vcf.fetch(chrom)

for v in vcf:

#for v in vcf_chrom:

#for v in vcf:
    d+=1
    #print(d)
#iterate through each variant line in vcf/chrom

    # coding variants
    if 'CSQ' not in v.info.keys(): continue
    ac=defaultdict(list)
    idict=(dict(v.info.items()))
    genotypes={}
    genotypes_m={}
    GQs={}
    DPs={}
    ADs={}
    GQs_m={}
    DPs_m={}
    ADs_m={}

    num_males=len(SexSample['M'])
    num_females=len(SexSample['F'])
	#Dictionaries with GT GQ and DP per sample
    if v.chrom == 'Y':
        genotypes=dict([(s,v.samples[s]['GT'],) for s in SexSample['M']])

        deletelist=[]
        
        for male_gts in genotypes:
            #print(genotypes[male_gts])
            if genotypes[male_gts][0]!=genotypes[male_gts][1]:
                #print(genotypes[male_gts][0],genotypes[male_gts][1])
                deletelist.append(male_gts)
        for delthese in deletelist:
            del genotypes[delthese]

        GQs=dict([(s,v.samples[s]['GQ'],) for s in genotypes.keys()])
        DPs=dict([(s,v.samples[s]['DP'],) for s in genotypes.keys()])
        ADs=dict([(s,v.samples[s]['AD'],) for s in genotypes.keys()])

    elif v.chrom =='X':
        genotypes=dict([(s,v.samples[s]['GT'],) for s in SexSample['F']])
        GQs=dict([(s,v.samples[s]['GQ'],) for s in genotypes.keys()])
        DPs=dict([(s,v.samples[s]['DP'],) for s in genotypes.keys()])
        ADs=dict([(s,v.samples[s]['AD'],) for s in genotypes.keys()])
        genotypes_m=dict([(s,v.samples[s]['GT'],) for s in SexSample['M']])
        deletelist=[]
        for male_gts in genotypes_m:
            if genotypes_m[male_gts][0]!=genotypes_m[male_gts][1]:
                deletelist.append(male_gts)
        for delthese in deletelist:
            del genotypes_m[delthese]
        GQs_m=dict([(s,v.samples[s]['GQ'],) for s in genotypes_m.keys()])
        DPs_m=dict([(s,v.samples[s]['DP'],) for s in genotypes_m.keys()])   
        ADs_m=dict([(s,v.samples[s]['AD'],) for s in genotypes_m.keys()])


    else:
        genotypes=dict([(s,v.samples[s]['GT'],) for s in samplePheno])
        GQs=dict([(s,v.samples[s]['GQ'],) for s in genotypes])
        DPs=dict([(s,v.samples[s]['DP'],) for s in genotypes])
        ADs=dict([(s,v.samples[s]['AD'],) for s in genotypes.keys()])

    deletelist=[]
    for aGQ in GQs:
        #print(aGQ)
        
        if(GQs[aGQ]<20 or not isinstance(GQs[aGQ], int)):
            deletelist.append(aGQ)
    for delthese in deletelist:
        del genotypes[delthese]	
        del GQs[delthese]
        #print(ADs)
        del ADs[delthese]
        del DPs[delthese]
    deletelist=[]
    for GTs in genotypes:
        GT1=genotypes[GTs][0]
        GT2=genotypes[GTs][1]
        if(ADs[GTs][GT1]<=2 or ADs[GTs][GT1]<=2):
            for delthese in deletelist:
                del genotypes[delthese]
                del GQs[delthese]
                del ADs[delthese]
                del DPs[delthese]

    deletelist=[]
    #print("GQs_M",GQs_m)
    for aGQ in GQs_m:
        #print(aGQ)
        if(GQs_m[aGQ]<20  or not isinstance(GQs_m[aGQ], int)):
            deletelist.append(aGQ)
    for delthese in deletelist:
        del genotypes_m[delthese]	
        del GQs_m[delthese]
        del ADs_m[delthese]
        del DPs_m[delthese]
    deletelist=[]
    for GTs in genotypes_m:
        GT1=genotypes_m[GTs][0]
        GT2=genotypes_m[GTs][1]
        if(ADs_m[GTs][GT1]<=2 or ADs_m[GTs][GT1]<=2):
            for delthese in deletelist:
                del genotypes_m[delthese]
                del GQs_m[delthese]
                del ADs_m[delthese]
                del DPs_m[delthese]
    #num_males=len(SexSample['M'])
    #num_females=len(SexSample['F'])
    callrate=0
    if v.chrom == 'Y':
        if(num_males==0):
            callrate=0
        else:
            callrate=float(len(genotypes))/float(num_males)
    else:
        called_genotypes=len(genotypes)+len(genotypes_m)
        callrate=float(called_genotypes)/float(num_males+num_females)
    if callrate<0.5:
        #print(len(genotypes),len(genotypes_m))
        #print("callrate < 0.5")
        continue
    
    #print(GQs)
    #print(GQs_m)
    #print("COntinue?",d)
    
    GENOPOS={}
    GENOPOS_m={}    
    if(len(genotypes)!=0):
	 
        geno=[genotypes[s] for s in genotypes.keys()]

        gq=pd.Series([GQs[s] for s in genotypes.keys()])
        dp=pd.Series([DPs[s] for s in genotypes.keys()])
        gq_m=pd.Series()
        dp_m=pd.Series()
        GENO=pd.Series(geno)

        AC=0

	    #set of samples that are reference
        #print(GENO[GENO==(0,0)].index)
        GENOPOS[0]=set(GENO[GENO==(0,0)].index)
        GENOPOS_m.setdefault(0,set())

        
        GQ_hist=''
        GQ_mid='' 
        DP_hist=''
        DP_mid=''
	    
	    #aref llele count 
        if(v.chrom == 'Y'):
            ac['AC_Adj'].append(geno.count((0,0)))
        else:
            ac['AC_Adj'].append(geno.count((0,0))*2)
        for i,alt_allele, in enumerate(v.alts):
            ac['AC_Adj'].append(0)
            ac['AF'].append(0)
            ac['AC_Hom'].append(0)
            ac['AC_Het'].append(0)
            if(v.chrom == 'Y'):
	    	    ac['AC_Hemi'].append(0)
            GENOPOS.setdefault(i+1,set())
            GENOPOS_m.setdefault(i+1,set())
        for i,alt_allele, in enumerate(v.alts):
	    	#set of samples that have 1 reference 
            GENOPOS[0].update(set(GENO[GENO==(0,i+1)].index))
	    	#set of samples that have 1 alt
            GENOPOS[i+1].update(set(GENO[GENO==(0,i+1)].index))
	    	#set of samples that have 2 alt
            GENOPOS[i+1].update(set(GENO[GENO==(i+1,i+1)].index))
	    	
	    	#homozygous alt counts
            if(v.chrom == 'Y'):
                #print("hemilen",str(i+1),len(ac['AC_Hemi']))
                ac['AC_Hemi'][i]+=(geno.count((i+1,i+1)))
                ac['AC_Adj'][i+1]+=(geno.count((i+1,i+1)))
        
            else:
                ac['AC_Hom'][i]+=(geno.count((i+1,i+1)))
                ac['AC_Adj'][i+1]+=(geno.count((i+1,i+1))*2)
	    	# het counts
            
            ac['AC_Het'][i]+=(geno.count((0,i+1)))
	    	# alt counts ref
            ac['AC_Adj'][0]+=(geno.count((0,i+1)))
	    	#alt ocounts alt
            ac['AC_Adj'][i+1]+=(geno.count((0,i+1)))
            if(v.chrom != 'Y'):
                for j,alt_allele_2, in enumerate(v.alts):
                        if i==j:
                            continue
                        else:
	    	        		#set of samples with 1 of the alts
                            GENOPOS[i+1].update(set(GENO[GENO==(i+1,j+1)].index))
                            GENOPOS[j+1].update(set(GENO[GENO==(i+1,j+1)].index))

                            ac['AC_Adj'][i+1]+=(geno.count((i+1,j+1)))
                            ac['AC_Adj'][j+1]+=(geno.count((i+1,j+1)))
                            ac['AC_Het'][i]+=(geno.count((i+1,j+1)))
                            ac['AC_Het'][j]+=(geno.count((i+1,j+1)))

	    #get study specific numbers
        for p in phenoSample:
            #print ('p',p)
            pheno_samples=set(genotypes.keys())&set(phenoSample[p])
            geno=[genotypes[s] for s in pheno_samples]

            if(v.chrom == 'Y'):
                ac['AC_%s'%POPS.loc[p]['code']].append(geno.count((0,0)))

            else:
                ac['AC_%s'%POPS.loc[p]['code']].append(geno.count((0,0))*2)
            for i,alt_allele, in enumerate(v.alts):
                ac['AC_%s'%POPS.loc[p]['code']].append(0)
                ac['AF_%s'%POPS.loc[p]['code']].append(0)
                ac['Het_%s'%POPS.loc[p]['code']].append(0)
                ac['Hom_%s'%POPS.loc[p]['code']].append(0)
                if(v.chrom == 'Y'):
	    	        ac['Hemi_%s'%POPS.loc[p]['code']].append(0)
            for i,alt_allele, in enumerate(v.alts):
                if(v.chrom == 'Y'):
                    ac['Hemi_%s'%POPS.loc[p]['code']][i]+=(geno.count((i+1,i+1)))
                    ac['AC_%s'%POPS.loc[p]['code']][i+1]+=(geno.count((i+1,i+1)))
                else:
                    ac['Hom_%s'%POPS.loc[p]['code']][i]+=(geno.count((i+1,i+1)))
                    ac['AC_%s'%POPS.loc[p]['code']][i+1]+=(geno.count((i+1,i+1))*2)
                ac['Het_%s'%POPS.loc[p]['code']][i]+=(geno.count((0,i+1)))
                ac['AC_%s'%POPS.loc[p]['code']][0]+=(geno.count((0,i+1)))
                ac['AC_%s'%POPS.loc[p]['code']][i+1]+=(geno.count((0,i+1)))
                if(v.chrom != 'Y'):
                    for j,alt_allele_2, in enumerate(v.alts):
                        if i==j:
                            continue
                        else:
                            ac['AC_%s'%POPS.loc[p]['code']][i+1]+=(geno.count((i+1,j+1)))
                            ac['AC_%s'%POPS.loc[p]['code']][j+1]+=(geno.count((i+1,j+1)))
                            ac['Het_%s'%POPS.loc[p]['code']][i]+=(geno.count((i+1,j+1)))
                            ac['Het_%s'%POPS.loc[p]['code']][j]+=(geno.count((i+1,j+1)))


      
    if(len(genotypes_m)!=0):
        #print("Some male Xs")
        if v.chrom == 'X':
            geno=[genotypes_m[s] for s in genotypes_m.keys()]
            if(geno.count((0,1))==0 and geno.count((1,1))==0):
                continue
            gq_m=pd.Series([GQs_m[s] for s in genotypes_m.keys()])
            dp_m=pd.Series([DPs_m[s] for s in genotypes_m.keys()])
            GENO=pd.Series(geno)
            
            GENOPOS_m.setdefault(0,set())
            
            GENOPOS_m[0].update(set(GENO[GENO==(0,0)].index))

            if(len(ac['AC_Adj'])==0):
                ac['AC_Adj'].append(geno.count((0,0)))
            else:
	    	    ac['AC_Adj'][0]+=(geno.count((0,0)))
            for i,alt_allele, in enumerate(v.alts):
                if(len(ac['AC_Adj'])<=i+1):
                    ac['AC_Adj'].append(0)
                if(len(ac['AF'])<i+1):
                    ac['AF'].append(0)
                if(len(ac['AC_Het'])<i+1):
                    ac['AC_Het'].append(0)
                ac['AC_Hemi'].append(0)
                GENOPOS_m.setdefault(i+1,set())
                #GENOPOS[i+1]=set()
            for i,alt_allele, in enumerate(v.alts):
        		#set of samples that have 1 reference 
                GENOPOS_m[0].update(set(GENO[GENO==(0,i+1)].index))
        		#set of samples that have 1 alt
                GENOPOS_m[i+1].update(set(GENO[GENO==(0,i+1)].index))
        		#set of samples that have 2 alt
                GENOPOS_m[i+1].update(set(GENO[GENO==(i+1,i+1)].index))
                #print("iter",i)
        		#homozygous alt counts
                ac['AC_Hemi'][i]+=(geno.count((i+1,i+1)))
                #print(ac['AC_Adj'])
                ac['AC_Adj'][i+1]+=(geno.count((i+1,i+1)))
	    	
            for p in phenoSample:
                pheno_samples=set(genotypes_m.keys())&set(phenoSample[p])
                geno=[genotypes_m[s] for s in pheno_samples]
                if(len(ac['AC_%s'%POPS.loc[p]['code']])==0):
                    ac['AC_%s'%POPS.loc[p]['code']].append(geno.count((0,0)))
                else:
                    ac['AC_%s'%POPS.loc[p]['code']][0]+=(geno.count((0,0)))
                for i,alt_allele, in enumerate(v.alts):
                    if(len(ac['AC_%s'%POPS.loc[p]['code']])<=i+1):
                        ac['AC_%s'%POPS.loc[p]['code']].append(0)
                    if(len(ac['AF_%s'%POPS.loc[p]['code']])<i+1):
                        ac['AF_%s'%POPS.loc[p]['code']].append(0)
                    if(len(ac['Het_%s'%POPS.loc[p]['code']])<i+1):
                        ac['Het_%s'%POPS.loc[p]['code']].append(0)
                    ac['Hemi_%s'%POPS.loc[p]['code']].append(0)
                for i,alt_allele, in enumerate(v.alts):
                    ac['AC_%s'%POPS.loc[p]['code']][i+1]+=(geno.count((i+1,i+1)))
                    ac['Hemi_%s'%POPS.loc[p]['code']][i]+=(geno.count((i+1,i+1)))
	    
     
         					
    for p in phenoSample:
        #pheno_samples=set(genotypes.keys())&set(phenoSample[p])
        ac['AN_%s'%POPS.loc[p]['code']].append(0)
        for i in ac['AC_%s'%POPS.loc[p]['code']]:
            ac['AN_%s'%POPS.loc[p]['code']][0]+=i
        ac['AC_%s'%POPS.loc[p]['code']].pop(0)
        for i,an, in enumerate(ac['AC_%s'%POPS.loc[p]['code']]):
    #for i in ac['AC_adj']:
        #print(i,an)
            if(float(ac['AN_%s'%POPS.loc[p]['code']][0])>0):
                ac['AF_%s'%POPS.loc[p]['code']][i]=(float(an))/(float(ac['AN_%s'%POPS.loc[p]['code']][0]))
            else:
                ac['AF_%s'%POPS.loc[p]['code']][i]=0

            



		
    
    GENOPOS.setdefault(0,set())
    GENOPOS_m.setdefault(0,set())
    FULL_GENOPOS=GENOPOS[0]
    FULL_GENOPOS_m=GENOPOS_m[0]

    for i,alt_allele, in enumerate(v.alts):
        GENOPOS.setdefault(i+1,set())		
        GENOPOS_m.setdefault(i+1,set())		
        FULL_GENOPOS.update(GENOPOS[i+1])
        FULL_GENOPOS_m.update(GENOPOS_m[i+1])
        gqtest=list(gq[GENOPOS[i+1]])
        gqtest.extend(list(gq_m[GENOPOS_m[i+1]]))
        
        #print("-",v.chrom,v.pos,gqtest)
        
        dptest=list(dp[GENOPOS[i+1]].fillna(0))
        dptest.extend(list(dp_m[GENOPOS_m[i+1]].fillna(0)))

        gq_hist = get_histogram_for_variant(gqtest)
        dp_hist = get_histogram_for_variant(dptest)
        GQ_hist+=','+'|'.join(map(str, gq_hist))
        #GQ_mid+=','+'|'.join(map(str, gq_midpoints))
        DP_hist+=','+'|'.join(map(str, dp_hist))
        #DP_mid+=','+'|'.join(map(str, dp_midpoints))
    
    gqtest=list(gq[FULL_GENOPOS])
    #print(gqtest)
    gqtest.extend(list(gq_m[FULL_GENOPOS_m]))
    dptest=list(dp[FULL_GENOPOS].fillna(0))
    dptest.extend(list(dp_m[FULL_GENOPOS_m].fillna(0)))
    #print('_',v.chrom,v.pos,gqtest)
	#make the histogram data for all variants
    gq_hist = get_histogram_for_variant(gqtest)
    dp_hist = get_histogram_for_variant(dptest)   

    new_info=";"+'GQ_HIST='+'|'.join(map(str, gq_hist))+GQ_hist+";"+'DP_HIST='+'|'.join(map(str, dp_hist))+DP_hist
    new_info=";"+'GQ_HIST='+'|'.join(map(str, gq_hist))+GQ_hist+";"+'DP_HIST='+'|'.join(map(str, dp_hist))+DP_hist

    
    ac['AN_Adj'].append(0)
    
    
    for i in ac['AC_Adj']:
        ac['AN_Adj'][0]+=i
   
	#REMOVING REFERENCE CALLS FROM ALLELE COUNTS
    ac['AC_Adj'].pop(0)

    
    for i,an, in enumerate(ac['AC_Adj']):
        if(float(ac['AN_Adj'][0])>0):
            ac['AF'][i]=(float(an))/(float(ac['AN_Adj'][0]))
        else:
            ac['AF'][i]=0    
    #print(ac)
    anything_there=0
    for AF in ac['AF']:
        anything_there=AF
    if(anything_there==0):
        #print("NO ALTS!")
        continue
    ACs=[]

    #print([','.join(map(lambda x: str(x), ac[k])) for k in ac])
    for key, value in ac.iteritems():
        ACs.append(key+"="+','.join([str(i) for i in value]))
    
    
    infs=[]
    for key, value in idict.iteritems():
        #print(key,value,type(value),len(value))
        if(type(value)==tuple):
            value=','.join(map(str,value))
            #print(key,value,type(value),len(value))
        infs.append(key+"="+str(value))
    ACs.extend(infs)
    INFO=';'.join(ACs)
    INFO+=new_info
    FILTER=','.join([v.filter[k].name for k in v.filter]) or '.'
    print('\t'.join(map(lambda x: str(x), [v.chrom,v.pos,v.id or '.',v.ref,','.join(v.alts),v.qual,FILTER,INFO])))



