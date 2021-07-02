#!/usr/bin/env python

import os
import subprocess
# import argparse
import config
import multiprocessing as mp
from multiprocessing import Pool


rawdatadir = config.get_rawdatadir()
resultdir = config.get_resultdir()
inputfile = config.get_inputfile()
plink = config.get_plink()
rscript = config.get_rscript()
reference = config.get_reference()
output = config.get_output()
refdir = config.get_refdir()
impute2 = config.get_impute2()
shapeit = config.get_shapeit()

os.system("module load GCC/6.4.0-2.28  OpenMPI/2.1.1 R/3.4.3")
os.system("module load GCC/6.4.0-2.28 Python/3.6.3")
os.system("R --version")
os.system("python --version")

# download the reference dataset of 1000 Genome


def refdat():
        # Phase 1
        #os.system("wget https://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz")
        #os.system("tar xzvf ALL_1000G_phase1integrated_v3_impute.tgz --directory " + str(reference))
        #os.system("wget https://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_annotated_legends.tgz")
        #os.system("tar xzvf ALL_1000G_phase1integrated_v3_annotated_legends.tgz --directory " + str(reference))
        # Phase 3
        os.system("wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz")
        os.system("tar xzvf 1000GP_Phase3.tgz --directory " + str(reference))


def preimpute():
        # answer = input("Do you want to use 1000 Genome Phase3 or Phase1? Please answer Phase3 or Phase1: ")
        # if answer == "Phase1":
        #         # create a file with all the SNP names that are in the reference set
        #         SNP = "SNP"
        #         os.system("for i in `seq 1 22`; do " + "gunzip -c " + str(reference) + "ALL_1000G_phase1integrated_v3_annotated_legends/ALL_1000G_phase1integrated_v3_chr${i}_impute.legend.gz | " + "gawk -v chr=${i} '$5==\"" + str(SNP) + "\" " + "{print chr\"" + " \"$2}' >> " + str(inputfile) + "snpsref.txt; done")
        # else:
        #         # create a file with all the SNP names that are in the reference set
        #         SNP = "Biallelic_SNP"
        #         os.system("for i in `seq 1 22`; do " + "gunzip -c " + str(reference) + "1000GP_Phase3/1000GP_Phase3_chr${i}.legend.gz | " + "gawk -v chr=${i} '$5==\"" + str(SNP) + "\" " + "{print chr\"" + " \"$2}' >> " + str(inputfile) + "snpsref.txt; done")
        # get a list of positions of SNPs that are in the target set
        # os.system("gawk '{print $1\"" + " \"$4}' " + str(inputfile) + "totaldata_extractqc.bim > " + str(inputfile) + "snpsraw.txt")
        # get SNPs that are in both the target set and reference set, to make the format corresponding to the --extract range option in plink
        #os.system("Rscript " + str(rscript) + "snpref_raw.R " + str(inputfile))
        # since some SNPs in target set but not in reference set,SNPs that are in both the target set and reference set need to be extracted from the target set, according to the physical position, not the SNP name
        os.system(str(plink) +" --noweb --bfile " + str(inputfile) + "totaldata_extractqc --extract range " + str(inputfile) + "duplicatesnp --make-bed --out " + str(inputfile) + "cleantotaldata_extractqc")


def preshapeit():
        #split 22 chromosomes into vcf format(shapeit4)
        #os.system("for chr in `seq 9 9`; do " + str(plink) + " --noweb --bfile " + str(inputfile) + "cleantotaldata_extractqc --chr $chr --recode vcf --snps-only just-acgt --keep " + str(inputfile) + "extract_samples4.txt --threads 10 --out " + str(output) + "topmed_cleantotaldata_extractqc.chr${chr}_4; done")
        #os.system("for chr in `seq 9 9`; do /dat/Plink-Based-GWAS/software/htslib-1.11/bgzip < " + str(output) + "topmed_cleantotaldata_extractqc.chr${chr}_4.vcf > " + str(output) + "topmed_cleantotaldata_extractqc.chr${chr}_4.vcf.gz --threads 10; done")
        #os.system("for chr in `seq 9 9`; do /dat/Plink-Based-GWAS/software/htslib-1.11/tabix -p vcf " + str(output) + "cleantotaldata_extractqc.chr${chr}.vcf.gz; done")

        # os.system("sudo chmod u+x " + str(shapeit))
        # answer = input("Do you want to use 1000 Genome Phase3 or Phase1? Please answer Phase3 or Phase1: ")

        ## dir = str(inputfile)
        ## filename = "preshapeit.txt"
        ## filename = "%s%s" % (dir, filename)
        ## preshapeit = open(filename, "w")
        ## preshapeit.writelines(align)
        ## preshapeit.close()

        # if answer == "Phase1":
        #         # alignment of the SNPs between the target set and reference set
        #         align = "for chr in `seq 1 22`; do " + str(shapeit) + " -check -B " + str(output) + "cleantotaldata_extractqc.chr${chr} -M " + str(reference) + str(refdir) + "genetic_map_chr${chr}_combined_b37.txt --input-ref " + str(reference) + str(refdir) + "ALL_1000G_phase1integrated_v3_chr${chr}_impute.hap.gz " + str(reference) + str(refdir) + "ALL_1000G_phase1integrated_v3_chr${chr}_impute.legend.gz " + str(reference) + str(refdir) + "ALL_1000G_phase1integrated_v3.sample --output-log " + str(output) + "chr${chr}.alignment;done"
        #         os.system(align)
        # else:
        #         # alignment of the SNPs between the target set and reference set (should try to parallelize it)
        #         align = "for chr in `seq 9 9`; do echo " + str(shapeit) + " -check -B " + str(output) + "cleantotaldata_extractqc.chr${chr} -M " + str(reference) + "1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt --input-ref " + str(reference) + "1000GP_Phase3/1000GP_Phase3_chr${chr}.hap.gz " + str(reference) + "1000GP_Phase3/1000GP_Phase3_chr${chr}.legend.gz " + str(reference) + str(refdir) + "1000GP_Phase3.sample --thread 10 --output-log " + str(inputfile) + "chr${chr}.alignment >> " + str(inputfile) + "preshapeit_commands.txt; done"
        #         os.system(align)

        preshapeit_commands = "cat " + str(inputfile) + "preshapeit_commands.txt | xargs -P1 -n1 -L1 /dat/Plink-Based-GWAS/software/shapeit &" 
        os.system(preshapeit_commands)


def doshapeit():   
        
        #using shapeit4
        # align = "for chr in `seq 9 9`; do echo " + "--input %scleantotaldata_extractqc.chr${chr}.vcf.gz" % (output) + " --map" + " %schr${chr}.b37.gmap.gz" % (reference) + " --region ${chr}" + " --output" + " %sshapeit/cleantotaldata_extractqc.chr${chr}.vcf.gz" % (resultdir) + " --thread 30 --log %sshapeit/cleantotaldata_extractqc.chr${chr}.log >> " % (resultdir) + str(inputfile) + "doshapeit_commands.txt; done"
        # os.system(align)
        # doshapeit_commands = "cat " + str(inputfile) + "doshapeit_commands.txt | xargs -n1 -P1 -L1 /dat/Plink-Based-GWAS/software/shapeit4/bin/shapeit4 &" 
        # os.system(doshapeit_commands)

        #using shapeit3
        # align = "for chr in `seq 9 9`; do echo " + "-B %scleantotaldata_extractqc.chr${chr}" % (output) + " -M" + " %s%sgenetic_map_chr${chr}_combined_b37.txt" % (reference,refdir) + " -O" + " %sshapeit/cleantotaldata_extractqc.chr${chr}_shapeit3" % (resultdir) + " --exclude-snp %schr9.alignment.snp.strand.exclude" % (inputfile) + " --thread 30 --cluster-size 5000 –early-stopping --fast -L %sshapeit/cleantotaldata_extractqc.chr${chr}_shapeit3.log >> " % (resultdir) + str(inputfile) + "doshapeit3_commands.txt; done"
        # os.system(align)
        # os.system("chmod u+x " + "/dat/Plink-Based-GWAS/software/shapeit3.r884.1")
        # doshapeit_commands = "cat " + str(inputfile) + "doshapeit3_commands.txt | xargs -n1 -P1 -L1 /dat/Plink-Based-GWAS/software/shapeit3.r884.1 &" 
        # os.system(doshapeit_commands)

        #using shapeit2
        # shapeitcommand = "for chr in `seq 21 -1 1`; do echo " + "--force --input-bed %scleantotaldata_extractqc.chr${chr}.bed" % (output) + " %scleantotaldata_extractqc.chr${chr}.bim" % (output) + " %scleantotaldata_extractqc.chr${chr}.fam" % (output) + " --input-map" + " %s%sgenetic_map_chr${chr}_combined_b37.txt" % (reference, refdir) + " --input-ref" + " %s%s1000GP_Phase3_chr${chr}.hap.gz %s%s1000GP_Phase3_chr${chr}.legend.gz %s%s1000GP_Phase3.sample" % (reference, refdir,reference, refdir,reference, refdir) + " --exclude-snp" + " %schr${chr}.alignment.snp.strand.exclude" % (output) + " --output-max" + " %sshapeit/cleantotaldata_extractqc.chr${chr}.phased" % (resultdir) + " --no-mcmc --thread 20 --output-log %sshapeit/cleantotaldata_extractqc.chr${chr}.phased >> " % (resultdir) + str(inputfile) + "doshapeit2_commands.txt; done"
        # os.system(shapeitcommand)

        doshapeit_commands = "for chr in `seq 6 -1 4`; do cat " + str(inputfile) + "doshapeit2_commands_${chr}.txt; done" + " | xargs -n1 -P1 -L1 /dat/Plink-Based-GWAS/software/shapeit &" 
        os.system(doshapeit_commands)


# count = mp.cpu_count()
# pool = mp.Pool(processes=count)
# pool.map(doshapeit, range(9, 10))


def imputation():
        #================================using impute4 from shapeit4
        ##convert vcf.gz file from shapeit4 to haps/sample file
        #os.system("/dat/Plink-Based-GWAS/software/plink2 --vcf /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.vcf.gz --export haps --out /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9_shapeit4_phased")
        #using shapeit4 to do the haplotype estimation
        ## convert .gen file to vcf.gz file(convert gen to plink to vcf)
        ##no!!!os.system("/dat/Plink-Based-GWAS/software/htslib-1.11/bgzip < /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.chunk2_imputed4.gen > /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_imputed4.gen.gz --threads 10")
        ##no!!!os.system("/dat/Plink-Based-GWAS/software/bcftools/bcftools convert /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_imputed4.vcf.gz --vcf-ids --chrom --threads 10 --gensample2vcf /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_imputed4")
        #os.system("/dat/Plink-Based-GWAS/software/plink --data -gen /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.chunk2_imputed4.gen --sample /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9_shapeit4.sample --hard-call-threshold 0.2 --oxford-single-chr 9 --make-bed --out /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.chunk2_imputed4")
        #os.system("/dat/Plink-Based-GWAS/software/plink --bfile /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.chunk2_imputed4 --recode vcf bgz --out /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.chunk2_imputed4")
        #os.system("/dat/Plink-Based-GWAS/software/htslib-1.11/tabix -p vcf /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.chunk2_imputed4.vcf.gz")
        ##convert .gen to .vcf.gz file through bcftool
        #os.system("/dat/Plink-Based-GWAS/software/bcftools/bcftools convert --gensample2vcf /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.chunk2_imputed4.gen --vcf-ids --chrom --threads 30 -Oz -o /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_imputed4_gen")
        
        # os.system("/dat/Plink-Based-GWAS/software/plink2 --vcf /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.chunk2_imputed4.vcf.gz --export haps --out /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_imputed4")
        #os.system("/dat/Plink-Based-GWAS/software/shapeit4/bin/shapeit4 --input /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.vcf.gz.chunk2_imputed.vcf.gz --map /dat/Plink-Based-GWAS/reference/chr9.b37.gmap.gz --region 9 --output /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_imputed4_haplotype_haps.vcf.gz --thread 30 --log /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_imputed4_haplotype_haps.log")
        #os.system("/dat/Plink-Based-GWAS/software/plink2 --vcf /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_imputed4_haplotype_haps.vcf.gz --export haps --out /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_imputed4_haplotype_haps")

        #======================================using impute4 from shapeit3
        ##doing impute4 imputation from shapeit3
        #os.system("sudo chmod u+x /dat/Plink-Based-GWAS/software/impute4.1.2_r300.3")
        #os.system("for task in `seq 1 1`; do cat " + str(resultdir) + "impute5/chr9_task_${task}_impute4; done" + " | xargs -n1 -P1 -L1 /dat/Plink-Based-GWAS/software/impute4.1.2_r300.3 &")
        #doshapeit_commands = "cat " + str(inputfile) + "doshapeit3_imputed_commands.txt | xargs -n1 -P1 -L1 /dat/Plink-Based-GWAS/software/shapeit3.r884.1 &" 
        #os.system(doshapeit_commands)

        #======================================using impute4 from shapeit2
        ##doing impute4 imputation from shapeit3
        # for chr in range(22, 23):
        #         os.system("mkdir %simpute4/chr%s" % (resultdir, chr))
        # os.system("sudo chmod u+x /dat/Plink-Based-GWAS/software/impute4.1.2_r300.3")
        #os.system("Rscript " + str(rscript) + "splitlarge.R " + str(inputfile) + " " + str(resultdir) + " " + str(impute4) + " " + str(reference) + " " + str(refdir))
        # os.system("for task in `seq 2 3`; do cat " + str(resultdir) + "impute4/chr22_task_${task}; done" + " | xargs -n1 -P5 -L1 /dat/Plink-Based-GWAS/software/impute4.1.2_r300.3 &")
        # doshapeit_commands = "cat " + str(inputfile) + "doshapeit3_imputed_commands.txt | xargs -n1 -P1 -L1 /dat/Plink-Based-GWAS/software/shapeit3.r884.1 &" 
        # os.system(doshapeit_commands)

        #======================================using impute5 from shapeit2
        ##doing impute5 imputation from shapeit2
        # for chr in range(1, 23):
        #         os.system("mkdir %simpute4/chr%s" % (resultdir, chr))
        ##convert haps to vcf.gz
        # os.system("for chr in `seq 9 9`; do /dat/Plink-Based-GWAS/software/htslib-1.11/bgzip < /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.phased.haps > " + "/dat/Plink-Based-GWAS/result/shapeit/shapeitcleantotaldata_extractqc.chr9 --threads 8; done")
        # os.system("for chr in `seq 9 9`; do /dat/Plink-Based-GWAS/software/bcftools/bcftools convert --threads 8 --hapsample2vcf /dat/Plink-Based-GWAS/result/shapeit/shapeitcleantotaldata_extractqc.chr9 -Oz -o /dat/Plink-Based-GWAS/result/shapeit/shapeitcleantotaldata_extractqc.chr9.vcf.gz; done")
        # os.system("for chr in `seq 9 9`; do /dat/Plink-Based-GWAS/software/plink2 --haps /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.phased.haps --export vcf /dat/Plink-Based-GWAS/result/shapeit/shapeitcleantotaldata_extractqc.chr9.phased.vcf.gz")

        # os.system("Rscript " + str(rscript) + "splitlarge.R " + str(inputfile) + " " + str(resultdir) + " " + str(impute5) + " " + str(reference) + " " + str(refdir))      
        # os.system("chmod u+x " + str(impute5))
        #os.system("for chr in `seq 9 9`; do /dat/Plink-Based-GWAS/software/htslib-1.11/tabix -p vcf " + str(resultdir) + "shapeit/cleantotaldata_extractqc.chr${chr}.vcf.gz; done")
        #os.system("for task in `seq 1 1`; do cat " + str(resultdir) + "impute5/chr9_task_${task}; done" + " | xargs -n1 -P1 -L1 /dat/Plink-Based-GWAS/software/impute5_v1.1.3/impute5_v1.1.3_static &")
            #os.system("for chunk in `seq 2 2`; do /dat/Plink-Based-GWAS/software/htslib-1.11/tabix -p vcf /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.vcf.gz.chunk${chunk}_short_imputed.vcf.gz; done")

        ## ACCRE
        # os.system("for chr in `seq 9 9`; do %s --haps %sshapeit/cleantotaldata_extractqc.chr${chr}.phased.haps --export vcf bgz --out %sshapeit/cleantotaldata_extractqc.chr9.phased.haps.vcf.gz; done" % (plink,resultdir,resultdir))
        # os.system("for chr in `seq 9 9`; do /scratch/zhans23/GWAS_92k/software/htslib-1.11/tabix -p vcf %sshapeit/cleantotaldata_extractqc.chr${chr}.phased.haps.vcf.gz; done" % (resultdir))
        ##os.system("for chr in `seq 9 9`; do /scratch/zhans23/GWAS_92k/software/htslib-1.11/bgzip < " + str(resultdir) + "shapeit/cleantotaldata_extractqc.chr9.phased.haps > " + str(resultdir) + "shapeit/cleantotaldata_extractqc.chr9.phased.hap.gz --threads 70; done")
        ##os.system("for chr in `seq 9 9`; do /scratch/zhans23/GWAS_92k/software/bcftools/bcftools convert --threads 70 --hapsample2vcf " + str(resultdir) + "shapeit/cleantotaldata_extractqc.chr9.phased -Oz -o " + str(resultdir) + "shapeit/cleantotaldata_extractqc.chr9_phased_haps.vcf.gz; done")
        # os.system("Rscript " + str(rscript) + "splitlarge.R " + str(inputfile) + " " + str(resultdir) + " " + str(impute5) + " " + str(reference) + " " + str(refdir)) 
        ##convert vcf.gz file to bed
        # os.system("for chr in `seq 9 9`; do %s --vcf %simpute5/chr9/cleantotaldata_extractqc.chr${chr}.chunk2_imputed.vcf.gz --make-bed --out %simpute5/chr9/cleantotaldata_extractqc.chr${chr}.chunk2_imputed; done" % (plink, resultdir,resultdir))


        #======================================using impute2 from shapeit2
        ##os.system("/dat/Plink-Based-GWAS/software/plink2 --vcf /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.vcf.gz --export haps --out /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9_shapeit4_phased")
        #os.system("for task in `seq 1 1`; do cat " + str(resultdir) + "impute5/chr9_task_${task}_impute2; done" + " | xargs -n1 -P1 -L1 /dat/Plink-Based-GWAS/software/impute2 &")


        #     for chr in range(14, 15):
        #            os.system("mkdir %simpute5/chr%s" % (resultdir, chr))
                    
        #     for chr in range(9, 10):
        #              phasedfile="cleantotaldata_extractqc.chr%s.vcf.gz" % (chr)
        #              os.system("zcat " + str(resultdir) + "shapeit/" + str(phasedfile) + " | gawk 'NR>10 {print $2}' > " + str(inputfile) + "positions%s" % (chr))

        #     os.system("Rscript " + str(rscript) + "splitlarge.R " + str(inputfile) + " " + str(resultdir) + " " + str(impute5) + " " + str(reference) + " " + str(refdir))
            
        #     os.system("chmod u+x " + str(impute5))

        #     os.system("for chr in `seq 9 9`; do /dat/Plink-Based-GWAS/software/htslib-1.11/tabix -p vcf " + str(resultdir) + "shapeit/cleantotaldata_extractqc.chr${chr}.vcf.gz; done")

            #os.system("for task in `seq 1 1`; do cat " + str(resultdir) + "impute5/chr9_task_${task}; done" + " | xargs -n1 -P1 -L1 /dat/Plink-Based-GWAS/software/impute5_v1.1.3/impute5_v1.1.3_static &")
            #os.system("for chunk in `seq 2 2`; do /dat/Plink-Based-GWAS/software/htslib-1.11/tabix -p vcf /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.vcf.gz.chunk${chunk}_short_imputed.vcf.gz; done")
            #doing shapeit4 again to estimate haplotype
            #os.system("/dat/Plink-Based-GWAS/software/shapeit4/bin/shapeit4 --input /dat/Plink-Based-GWAS/result/impute5/chr9/cleantotaldata_extractqc.chr9.vcf.gz.chunk2_short_imputed.vcf.gz --map /dat/Plink-Based-GWAS/reference/chr9.b37.gmap.gz --region 9 --output /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_short_haplotype.vcf.gz --thread 30 --log /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_short_haplotype.log")
            #os.system("/dat/Plink-Based-GWAS/software/plink2 --vcf /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_short_haplotype.vcf.gz --export haps --out /dat/Plink-Based-GWAS/result/shapeit/cleantotaldata_extractqc.chr9.chunk2_short_haplotype_haps")

        #     for task in range(2,12):
        #             doimputation_commands = "cat " + str(resultdir) + "impute5/chr1_task_" + str(task) + " | xargs -n1 -P5 -L1 /dat/Plink-Based-GWAS/software/impute5_v1.1.3/impute5_v1.1.3_static &" 
        #             os.system(doimputation_commands)

        #     #os.system("gawk '{print $3}' " + str(resultdir) + "shapeit/" + str(phasedfile) + ".haps > " + str(inputfile) + "positions%s" % (chr))


        #==================================ACCRE
        # for chr in range(1, 23):
        #         os.system("mkdir %simpute2/chr%s" % (resultdir, chr))
        # os.system("chmod u+x " + str(impute4))
        # os.system("Rscript " + str(rscript) + "splitlarge.R " + str(inputfile) + " " + str(resultdir) + " " + str(impute2) + " " + str(reference) + " " + str(refdir))
        account = input("Please input the account: ")
        mail = input("Please input the mail: ")
        cpus = input("Please input cpus: ")
        time = input("Please input time: ")
        memory = input("Please input memory: ")

        for chr in range(8,9):
                for i in range(1,50):
                        dir = str(inputfile)
                        filename = "IMPUTE_TASK_%s_%s.slurm" % (chr, i)
                        filename = "%s%s" % (dir, filename)
                        impute = open(filename, "w")
                        impute.write("#!/bin/bash\n") 
                        imputelist1 = ["#SBATCH --account=%s\n" % (account), "#SBATCH --mail-user=%s\n" % mail, "#SBATCH --mail-type=ALL\n", "#SBATCH --ntasks=1\n", "#SBATCH --cpus-per-task=%s\n" % cpus, "#SBATCH --time=%s\n" % time, "#SBATCH --mem=%s\n" % memory, "#SBATCH --output=imputejob_%s_%i.out\n" % (chr,i)]
                        #imputelist2 = ["cat %simpute4/chr%s_task_%s | xargs -n1 -P3 -L1 %s &" % (resultdir, chr, i, impute4)]
                        # imputelist2 = ["cat %simpute2/chr%s_task_%s | parallel" % (resultdir, chr, i)]
                        imputelist2 = ["bash %simpute2/chr%s_task_%s" % (resultdir, chr, i)]
                        impute.writelines(imputelist1)
                        impute.writelines(imputelist2)
                        impute.close()
        # os.system("for task in `seq 2 3`; do cat " + str(resultdir) + "impute4/chr22_task_${task}; done" + " | xargs -n1 -P5 -L1 /dat/Plink-Based-GWAS/software/impute4.1.2_r300.3 &")
        # os.system("module load GCC/5.4.0-2.26  OpenMPI/1.10.3 R")
        # os.system("Rscript %ssplit.R %s %s %s %s %s" % (rscript, inputfile, resultdir, impute2, reference, refdir))


def imputechunk():
        os.system("chmod u+x " + str(impute2))
        # imputechunk = ("/dat/Plink-Based-GWAS/result/impute2/chr1_task_1")
        # os.system("cat %s" % imputechunk + " | parallel")

            #os.system("cd /dat/Plink-Based-GWAS/result/impute5/chr14/")

            #os.system("ls *.vcf.gz > vcfout.list")

            #index each chunk file
            #os.system("for chunk in `seq 1 22`; do /dat/Plink-Based-GWAS/software/htslib-1.11/tabix -p vcf /dat/Plink-Based-GWAS/result/impute5/chr14/cleantotaldata_extractqc.chr14.vcf.gz.chunk${chunk}_imputed.vcf.gz; done")
            #merge all chunks for each chromosome
        #    os.system("/dat/Plink-Based-GWAS/software/bcftools/bcftools merge -m none --file-list /dat/Plink-Based-GWAS/result/impute5/chr1/vcfout.list -Oz -o chr4_imputed.vcf.gz")
            #os.system("/dat/Plink-Based-GWAS/software/bcftools/bcftools concat --file-list /dat/Plink-Based-GWAS/result/impute5/chr14/vcfout.list --output-type z --output /dat/Plink-Based-GWAS/result/impute5/chr14/chr14_imputed.vcf.gz")

#refdat()
#preimpute()
#preshapeit()
#doshapeit()
#imputation()
# imputechunk()
for chr in range(7,8):
          for i in range(30,31):
                  for j in range(1,4):
                          os.system("sbatch " + str(inputfile) + "IMPUTE_TASK_%s_%s_%s.slurm" % (chr,i,j))

# os.system("for chr in `seq 9 9`; do " + str(plink) + " --bfile " + str(resultdir) + "impute4/cleanchr_qc --chr $chr --make-bed --out " + str(resultdir) + "cleanchr_qc_chr$chr; done")
        









