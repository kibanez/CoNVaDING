'''
Created on 22/07/2016

@author: kibanez
'''

#!/usr/bin/python

import sys, shlex , os,subprocess, shutil,copy,re

from subprocess import Popen , PIPE



import ConfigParser

import optparse

import logging

from tempfile import mkstemp
from shutil import move
from os import remove, close



######################################################################

class OptionParser(optparse.OptionParser):

    def check_required (self, opt):

        option = self.get_option(opt)

        atrib = getattr(self.values, option.dest)
        
        if atrib is None:
#            self.error("%s option not supplied" % option)
            return False
        else:
            return True
            

######################################################################

def read_cfg_file(cfg_filename):
    
    fi = open(cfg_filename,'r')
    
    config = ConfigParser.ConfigParser()
    config.readfp(fi)
    
    hash_cfg = {}
        
    for field in config.options('INPUT'):
        hash_cfg[field] = config.get('INPUT',field)
   
    for field in config.options('OUTPUT'):
        hash_cfg[field] = config.get('OUTPUT',field)
     
    for field in config.options('SOFTWARE'):
        hash_cfg[field] = config.get('SOFTWARE',field)
        
    fi.close()
    
    return hash_cfg

#######################################################################

# functions that replaces in a given file, the pattern for the substitution
def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    close(fh)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

#######################################################################

def run(argv=None):
    
    if argv is None: argv = sys.argv    
   
    parser = OptionParser(add_help_option=True,description="The script performs CNV estimation within the regions of interest following CoNVaDING strategy\n")
    
    parser.add_option("--cfg",default=None,help="Config file with the complete information of the target regions and paths of the files needed for the calling",dest="f_cfg")

                    
    # Se leen las opciones aportadas por el usuario
    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)
    
    if not parser.check_required("--cfg"):
        raise IOError('The cfg file does not exist')
        
               
    try:
        
        if options.f_cfg <> None:
            
            cfg_file = options.f_cfg        
          
            if not os.path.exists(cfg_file):
                raise IOError('The file %s does not exist' % (cfg_file))
            
            hash_cfg = read_cfg_file(cfg_file)

            # INPUT           
            alignment_path  = hash_cfg.get('alignment_path','')
            l_samples  = hash_cfg.get("sample_names",'').split(',')            
            analysis_bed = hash_cfg.get('analysis_bed','')

            # OUTPUT
            results_path = hash_cfg.get('results','')

            # SOFTWARE (CoNIFER)
            convading_path = hash_cfg.get('convading_path','')

            
            if not os.path.exists(alignment_path):
                raise IOError('The alignment path does not exist. %s' % (alignment_path))

            if not os.path.exists(convading_path):
                raise IOError('The CoNVaDING main folder does not exist. %s' % (convading_path))
                
            if not os.path.isfile(analysis_bed):
                raise IOError('The file does not exist. %s' % (analysis_bed))
            
            if not os.path.exists(results_path):
                os.mkdir(results_path)
            
                
            #Configure logger
            formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
            console = logging.StreamHandler()
            console.setFormatter(formatter)
            console.setLevel(logging.INFO)
            logger = logging.getLogger("preprocess")
            logger.setLevel(logging.INFO)
            logger.addHandler(console)
            
            l_bams = []
            for bam_f in l_samples:
                abs_path = os.path.join(alignment_path,bam_f)
                if not os.path.exists(abs_path):
                    raise IOError("The bam file does not exist. Check if the introduced path is correct: %s" %(abs_path))
                else:
                    l_bams.append(abs_path)
                
            logger.info("CNV estimation will be done in the following files: \n %s \n" %(l_bams))
        
        
            convanding_script = os.path.join(convading_path,"CoNVaDING.pl")
            
            # 1 - Convading runs under bam files with the following nomenclature: 1,2,3,4...X,Y (instead of chr1,chr2,chr3,,..). And We must change also the bed file
            nomenclature_change_script = os.path.join(convading_path,"changing_bam_chr_nomenclature.sh")
            
            alignment_path2 = alignment_path + "_nomenclature"
            
            if not os.path.exists(alignment_path2):
                os.mkdir(alignment_path2)
            
            args = shlex.split("%s %s %s" %(nomenclature_change_script,alignment_path,alignment_path2))
            #subprocess.call(args)
            
#             nomenclature_output    = subprocess.Popen(args, stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,close_fds=True,bufsize=-1)
#             (trash,logdata) = nomenclature_output.communicate()   
#             nomenclature_output.close()
#             
#             if logdata <> "":
#                 if logdata.lower().find("error") <> -1:
#                     raise RuntimeError("c_CoverageControl.perform_coverage_control_bed: Error in bedtools genomecov: %s" % (logdata))
            
            
            #os.system("%s %s %s" %(nomenclature_change_script,alignment_path,alignment_path2))
            
            logger.info("%s %s %s" %(nomenclature_change_script,alignment_path,alignment_path2))
            
                              
            # 2 - StartWithBam => without remoing duplicates (already done)
            #  perl ./CoNVaDING.pl -mode StartWithBam -inputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/align_nomenclature/ 
            # -outputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/normalizedCoverageOutput/ 
            # -controlsDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/align_nomenclature/ 
            # -bed /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/DE_Roche_V4_NC_refseq_s20_update_june2015.bed --useSampleAsControl
    
            output1 = os.path.join(results_path,"normalizedCoverageOutput")
            if not os.path.exists(output1):
                os.mkdir(output1)
           
            
            args = shlex.split("perl %s -mode StartWithBam -inputDir %s -outputDir %s -controlsDir %s -bed %s --useSampleAsControl" %(convanding_script,alignment_path2,output1,alignment_path2,analysis_bed))
            subprocess.call(args)
            
            
    	    # 3 - replace the nan and -nan values with 0
            # sed -i -- 's/nan/0/g' *.txt
            # sed -i -- 's/-nan/0/g' *.txt
            # sed -i -- 's/-//g' *.txt
            
            l_norm_files = filter(lambda x: x.endswith(('.txt')), os.listdir(output1))
            
            for f in l_norm_files:
                f = os.path.join(output1,f)
                
                if not os.path.isfile(f):
                    raise IOError('The normalized file does not exist. %s' % (f))
                
                replace(f,"nan","0")
                replace(f,"-nan","0")
                replace(f,"-","")
            
            # 4 - Usage of the first five columns of the file as input for the next step of the pipeline. ==> Y ESTAS SERAN la entrada para ==> -mode StartWithAvgCount
            select_columns_script = os.path.join(convading_path,"select_columns.sh")
            
            output2 = output1 + "_columns"
            
            if not os.path.exists(output2):
                os.mkdir(output2)
                
            args = shlex.split("%s %s %s" %(select_columns_script,output1,output2))
            #####subprocess.call(args,shell=True)
            os.system("%s %s %s" %(select_columns_script,output1,output2))
            
            # 5 - StartWithAvgCount
            #perl ./CoNVaDING.pl -mode StartWithAvgCount -inputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/normalizedCoverageOutput/AverageCount/ 
            # -outputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/normalizedCoverageOutput_fromAverageCount 
            # -controlsDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/normalizedCoverageOutput/AverageCount/ 
            # -bed /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/DE_Roche_V4_NC_refseq_s20_update_june2015.bed --useSampleAsControl
            
            output3 = os.path.join(results_path,"normalizedCoverageOutput_fromAverageCount")
            if not os.path.exists(output3):
                os.mkdir(output3)
                
            args = shlex.split("perl %s -mode StartWithAvgCount -inputDir %s -outputDir %s -controlsDir %s -bed %s --useSampleAsControl" %(convanding_script,output2,output3,output2,analysis_bed))
            
            subprocess.call(args)
            
            
            # 6 -  StartWithMatchScore
            # perl ./CoNVaDING.pl -mode StartWithMatchScore -inputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/normalizedCoverageOutput_fromAverageCount/sample/ 
            # -controlsDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/normalizedCoverageOutput_fromAverageCount/controls/ 
            # -bed /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/DE_Roche_V4_NC_refseq_s20_update_june2015.bed 
            # -outputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/bestMatchScore 
            # -controlSamples 20
            
            
            # For each sample, we do the CNV calling
            l_bestMatch = []
            l_cov = os.listdir(output1)
            
            # l_go contains a list with the samples that can be continued through the CNV calling
            l_go = []
            for f_avg in l_cov:
                
                # list_others
                l_others = copy.copy(l_cov)
                l_others.remove(f_avg)
                
                  
                f_name = f_avg.split("_")[0] 
                #output4 = os.path.join(results_path,"bestMatchScore_"+f_name)
                output4 = os.path.join(results_path,f_name)
                
                if not os.path.exists(output4):
                    os.mkdir(output4)
                    
                output4 = os.path.join(output4,'bestMatchScore')
                if not os.path.exists(output4):
                    os.mkdir(output4)
                
                l_bestMatch.append(output4)
                
                # sample tmp folder
                tmp_sample = os.path.join(output4,"sample")
                if not os.path.exists(tmp_sample):
                    os.mkdir(tmp_sample)
                shutil.copy(os.path.join(output3,f_avg.split('.txt')[0] + ".normalized.coverage.txt"),tmp_sample) 
                

                # theOthers (controls) tmp folder
                tmp_controls = os.path.join(output4,"controls")
                if not os.path.exists(tmp_controls):
                    os.mkdir(tmp_controls)
                for c in l_others:
                    shutil.copy(os.path.join(output3,c.split('.txt')[0] + ".normalized.coverage.txt"),tmp_controls)

                num_controls = 20
                
                args = shlex.split("perl %s -mode StartWithMatchScore -inputDir %s -controlsDir %s -bed %s -outputDir %s -controlSamples %s" %(convanding_script,tmp_sample,tmp_controls,analysis_bed,output4,str(num_controls)))
                subprocess.call(args)
                
                tmp_file = "/tmp/kakafuti.log"
                f_log = open(tmp_file,'w')
                call_sal = Popen(args,stdin=PIPE,stdout=f_log,stderr=PIPE,close_fds=True,bufsize=1)
                call_sal.communicate()
                call_sal.wait()
                f_log.close()
                
                logdata = open(tmp_file).readlines()
                
                for item in logdata:
                    items = re.findall("Selecting best.*$",item,re.MULTILINE)
                    if len(items) > 0:
                        print items
                        aux = ''.join(items).split(' ')
                        logger.info("Control samples used in the %s analysis : %s" %(f_avg,aux[2]))
                        if aux[2] <> '0':
                            l_go.append(output4)
                            
            
                if os.path.isfile(tmp_file):
                    os.remove(tmp_file)
            
            # 7- StartWithBestScore
            #  perl ./CoNVaDING.pl -mode StartWithBestScore -inputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/bestMatchScore 
            # -outputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/bestScore

            # For each sample, we do the CNV calling
            
            # Only whether there is a sample with good or enough quality to continue in the CNV analysis
            if l_go <> [] :

                for bestMatch in l_bestMatch:

                    f_name = '/'.join(bestMatch.split('/')[:-1])  
                    
                    output_bestScore = os.path.join(f_name,'bestScore')
                              
                    if not os.path.exists(output_bestScore):
                        os.mkdir(output_bestScore)
                
                    args = shlex.split("perl %s -mode StartWithBestScore -inputDir %s -outputDir %s" %(convanding_script,bestMatch,output_bestScore)) 
                    subprocess.call(args)
                
                # 8 -  GenerateTargetQcList
                #  perl ./CoNVaDING.pl -mode GenerateTargetQcList -controlsDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/normalizedCoverageOutput_fromAverageCount/controls/ 
                # -outputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/targetQClist 
                # -inputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/normalizedCoverageOutput_fromAverageCount/controls/ 
                # -controlSamples 20
                
                output6 = os.path.join(f_name,"TargetQcList")
                
                args = shlex.split("perl %s -mode GenerateTargetQcList -controlsDir %s -outputDir %s -inputDir %s -controlSamples %s" %(convanding_script,'',output6,'',str(num_controls)))
                subprocess.call(args)
                
                # 9 -  CreateFinalList
                # perl ./CoNVaDING.pl -mode CreateFinalList -inputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/bestScore 
                # -outputDir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/FinalList 
                # -targetQcList /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNVaDING/DE/2014-06-14/results/targetQClist/targetQcList.txt
    
                output7 = os.path.join(f_name,"FinalList")
                
                target_list = os.path.join(output6,"targetQcList.txt")
                args = shlex.split("perl %s -mode CreateFinalList -inputDir %s -outputDir %s -targetQcList %s" %(convanding_script,f_name,output7,target_list))
                subprocess.call(args)
            
            else:
                logger.info("There quality control is not enough to do CoNVaDING CNV calling \n ")
            
            logger.info("CoNVaDING CNV estimation done! \n")    
        
    except:
        print >> sys.stderr , '\n%s\t%s' % (sys.exc_info()[0],sys.exc_info()[1])
        sys.exit(2)

############################################################################333

if __name__=='__main__':
    
    run()



