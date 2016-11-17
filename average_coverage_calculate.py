

import tabix
import sys
import numpy
#import pandas





#sub to reduce coverage file in to many chunks. 
def chunkIt(seq, num):
    avg = seq / float(num)
    out = []
    last = 1.0

    while last < seq:
        out.append((int(last),int(last + avg)))
        last+=avg+1

    return out

  
#take a file with a list of paths to coverage files produced by the GATK. If you did it all in one go this will be a single file, if you split the coverage calculation in to multiple samples then it will have one GATK coverage output file path per line
  
coverage_file=(sys.argv[1])
coverage_list=open(coverage_file)

#coverage_list=open("/scratch/DGE/MOPOPGEN/dchubb/bedoptcheck/coverage/output/coverage2.list")
covfiles=[]
for i in coverage_list.read().splitlines():
    covfiles.append(tabix.open(i))

	
#buildsize=open("/data/cancgene-guard2/dchubb/chromosome_sizes.txt")

#use a capture file for the regions you want to calculate coverage over

capture_file=(sys.argv[2])
print capture_file
capture=open(capture_file)
outfile=open(coverage_file+".cov",'w')



positions={}

capture_regions=[]

for i in capture.read().splitlines():
    chr,posa,posb=i.split('\t')[0:3]
    capture_regions.append([chr,int(posa),int(posb)])

#print capture_regions
lenlist=(1,5,10,15,20,25,30,50,100)  
#outfile.write("#chrom\tpos\tmean\tmedian\t1\t5\t10\t15\t20\t25\t30\t50\t100\n") 
for i in capture_regions:
    cap={}
    for covfile in covfiles:
        record=covfile.query(i[0],int(i[1]),int(i[2]))
        for r in record:
            cap.setdefault(int(r[1]),[]).extend(r[2:])
    for k in sorted(cap.keys()):
        nr=numpy.array(cap[k]).astype(int)
        total=len(nr)
        mean=float("{0:.2f}".format(numpy.mean(nr)))
        median=float("{0:.2f}".format(numpy.median(nr)))
        depth_pcs=[]
        for lens in lenlist:
             number=sum(x>=lens for x in nr)
             pc=float("{0:.2f}".format(float(number)/float(total)))
             depth_pcs.append(pc)
        
        d="\t".join(str(x) for x in depth_pcs)
        outfile.write(i[0]+"\t"+str(k)+"\t"+str(mean)+"\t"+str(median)+"\t"+d+"\n")
        #print i[0]+"\t"+str(k)+"\t"+str(mean)+"\t"+str(median)+"\t"+d
        #positions.setdefault(chr,{})[k]=str(mean)+"\t"+str(median)+"\t"+d

    