'''
Created on 01/07/2016

@author: kibanez
'''

#!/usr/bin/python

import sys, getopt, os, subprocess,logging


def run(argv):

    try:
       
        opts, args = getopt.getopt(argv,"hi:o:",["idir=","odir="])

    except getopt.GetoptError:
        
        print 'changing_bam_chr_nomenclature.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
                    
    # Se leen las opciones aportadas por el usuario

    for opt, arg in opts:
        
        if opt == '-h':
            print 'changing_bam_chr_nomenclature.py -i <inputfile> -o <outputfile>'
            sys.exit()
           
        elif opt in ("-i", "--idir"):
            inputfile = arg
           
        elif opt in ("-o", "--odir"):
            outputfile = arg
           
    #Configure logger
    formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    console.setLevel(logging.INFO)
    logger = logging.getLogger("preprocess")
    logger.setLevel(logging.INFO)
    logger.addHandler(console)

    logger.info("Input file where the bam files (with chr nomenclature) are -->  %s" %(inputfile))
    
    logger.info("Output file where the bam files (with the new chr nomenclature) well be  -->  %s" %(outputfile))

    logger.info("The change in the chr nomenclature in the bam files starts now...")

    if not os.path.exists(inputfile):
        raise IOError("The path does not exist. %s" % (inputfile))
    
    if not os.path.exists(outputfile):
        os.mkdir(outputfile)

    # bam files (ending with *.bam)
    
    l_bam = []
    for file in os.listdir(inputfile):
        if file.endswith(".bam"):
            l_bam.append(os.path.join(inputfile,file))
            
    for b in l_bam:
        b_out = b + "chr_fixed.bam"
        #sh.Command("samtools view -h input.bam | awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0} $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -bS - > OUTPUT.bam")
                
        p1 = subprocess.Popen(["samtools view -h %s" %(b)], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["awk", 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0} $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0}'], stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(["sed 's/SN:/SN:chr/g'"], stdin=p2.stdout, stdout=subprocess.PIPE)
        p4 = subprocess.Popen(["sed 's/chrMT/chrM/g'"], stdin=p3.stdout,stdout=subprocess.PIPE)
        p5 = subprocess.Popen(["samtools view -bS - > %s" %(b_out)], stdin=p4.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        output,err = p5.communicate()
                    
    #samtools view -h input.bam | awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0} $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -bS - > OUTPUT.bam
    
    

    logger.info("The chr nomenclature in the bam files have been changed! ")    
        
    
############################################################################333

if __name__=='__main__':
    
    run(sys.argv[1:])



