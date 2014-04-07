#-*- coding: utf-8 -*-

'''
GATK-based variant calling pipeline

It implements a workflow pipeline for next generation
sequencing variant detection using the Broad Institute's GATK
for variant calling.

The pipeline is configured by an options file in a python file
including the actual command which are run at each stage.

'''
import logging
from collections import defaultdict
from utils import parse_and_link
from commands import *
import subprocess
import sys
import os
import argparse
from yaml import load
import glob

fastq_metadata = defaultdict(dict)

def get_logger(logfile, level=logging.INFO):
    log_handler = open(logfile, 'w')
    logging.basicConfig(stream=log_handler, level=level)

    return log_handler

def get_fastq_files(fc_dir, work_dir):
    fastq_files = []
    fastqz_files = glob.glob('%s/*.fastq.gz' % fc_dir)

    if len(fastqz_files) < 1:
            exit('At least one sequence file must be specified')

    #now let's parse the metadata from each fastq.gz input files and
    #construct the symbolic links to them.
    for fastqz in fastqz_files:
        symb_link = parse_and_link(fastqz, work_dir, fastq_metadata)
        fastq_files.append(symb_link)

    return fastq_files

def extract_fastq_files(sample_folder):
    gz_files = glob.glob('%s/*.fastq.gz' % sample_folder)
    for gz_file in gz_files:
        subprocess.call(['gzip', '-f', '-d', gz_file])
    return glob.glob('%s/*.fastq' % sample_folder)

def next_sample(fastqz_files):
    all_sample_folders =  [fastq_metadata[os.path.splitext(os.path.basename(fastqz_file))[0]]['out_dir']
                                        for fastqz_file in fastqz_files]
    for sample in list(set(all_sample_folders)):
        yield sample

def run(global_config, fc_dir, work_dir, tools_dir, workflow_config, reference, dbsnp, is_picard_to_bam):
    #1. Get all fasta files and check if there is at least one to process.
    sequence_files = get_fastq_files(fc_dir, work_dir)

    #2. Create reference database
    make_reference_database(workflow_config['indexer']['command'],'bwtsw', reference)

    #3.Index the reference
    index_reference(reference)

    for sample_folder in next_sample(sequence_files):
        #4. create output subdirectories
        fastqc_dir = make_output_dir(os.path.join(sample_folder, 'fastqc'))
        sambam_dir =  make_output_dir(os.path.join(sample_folder, 'alignments'))
        variant_dir = make_output_dir(os.path.join(sample_folder, 'variant_calls'))
        coverage_dir = make_output_dir(os.path.join(sample_folder, 'coverage'))
        annovar_dir  = make_output_dir(os.path.join(sample_folder, 'annovar'))
        results_dir =  make_output_dir(os.path.join(sample_folder, 'results'))

        #5.Precheck -  extracting fastqfiles.
        sequence_files = sorted(extract_fastq_files(sample_folder))
        #5.Run fastqc on each fastq file.
        #fastqc(workflow_config['fastqc']['command'],  sequence_files, fastq_metadata, fastqc_dir)
        #6. Align sequence to the  reference database.
        if len(sequence_files) == 2:
            #two paired-end fastq files alignment.
            fastq_file, pair_file = sequence_files
            if  can_pipe(workflow_config['seqtk']['command'], fastq_file):
                sam_output = align_with_mem(workflow_config['bwamem']['command'], global_config['bwa']['threads'],
                                    reference, fastq_file, pair_file, fastq_metadata, sambam_dir)
            else:
                #two paired-end fastq files alignment not piped.
                sai_fastq_file = align(workflow_config['aligner']['command'], global_config['bwa']['threads'],
                                    reference, fastq_file, fastq_metadata, sambam_dir)
                sai_pair_file = align(workflow_config['aligner']['command'], global_config['bwa']['threads'],
                                    reference, pair_file, fastq_metadata, sambam_dir)

                sam_output = alignPE2sam(workflow_config['sampe']['command'], reference, fastq_file,
                                    pair_file, sai_fastq_file, sai_pair_file, fastq_metadata, sambam_dir)

        elif len(sequence_files) == 1:
            #one fastq file alignment.
            fastq_file = sequence_files[0]
            sai_fastq_file = align(workflow_config['aligner']['command'], global_config['bwa']['threads'],
                    reference, fastq_file, fastq_metadata, sambam_dir)
            sam_output = align2sam(workflow_config['samse']['command'], reference, fastq_file, sai_fastq_file,
                                    fastq_metadata, sambam_dir)
        else:
            print('Multiple or no fastq files, please check the input files at %s' % sample_folder)
            continue

        #7. Convert SAM to BAM
        if (is_picard_to_bam):
            #if it is picard the tool to convert to bam
            seq_bam = samP2bam(workflow_config['samP2bam']['command'], global_config['picard']['jvm_opts'],
                    os.path.join(tools_dir, 'picard-tools-1.109'), sam_output, sambam_dir)
        else:
            seq_bam = samS2bam(workflow_config['samS2bam']['command'], global_config['samtools']['memory'],
                    global_config['samtools']['threads'], sam_output, sambam_dir)
            #After using samtools to convert to SAM let's use the samtools index to index it.
            index_bam = indexbam(workflow_config['bamindexer']['command'], seq_bam, sambam_dir)

        #7. Mark PCR Duplicates
        marked_seq_bam = dedup(workflow_config['markduplicates']['command'], global_config['picard']['jvm_opts'],
                    os.path.join(tools_dir, 'picard-tools-1.109'), seq_bam, sambam_dir)

        '''
        #8. Find suspect intervals for realignment.
        bam_list = realign_intervals(workflow_config['realigner']['command'], global_config['gatk']['jvm_opts'],
                    tools_dir, reference, marked_seq_bam, work_dir)
        #9. Run local realignment around indels.
        realigned_bam = realign(workflow_config['indelrealigner']['command'], global_config['gatk']['jvm_opts'],
                    tools_dir, reference, marked_seq_bam, bam_list, work_dir)
        #10. Fix mate information
        realigned_bam = fix_mate(workflow_config['fixmates']['command'], global_config['picard']['jvm_opts'],
                    os.path.join(tools_dir, 'picard-tools-1.109'), realigned_bam, work_dir)
        #11. Count Covariates
        recal_file = base_qual_recal_count(workflow_config['countcovariates']['command'], global_config['gatk']['jvm_opts'],
                    tools_dir, reference, dbsnp, realigned_bam, work_dir)
        #12. Table Recalibration
        realigned_bam = base_qual_recal_tabulate(workflow_config['recaltabulate']['command'], global_config['gatk']['jvm_opts'],
                    tools_dir, reference, recal_file, realigned_bam, work_dir'

        #13. Call Snps
        #output_vcf = call_snps(workflow_config['callSNPs']['command'], global_config['gatk']['jvm_opts'], global_config['gatk']['threads'],
        #            tools_dir, reference, dbsnp, '10.0', '10.0', '5000', '3', realigned_bam, work_dir)
        output_vcf = '/Users/marcelcaraciolo/Projects/genomika/github/nextgen-pipeline/3806_S10_L001_R1_001.marked.realigned.fixed.recal.vcf'
        #14. Filter Snps
        filtered_vcf = filter_snps(workflow_config['filterSNPs']['command'], global_config['gatk']['jvm_opts'],
                    tools_dir, reference, output_vcf, '', work_dir)
        #15.Converting the vcf to .annovar format file
        annovar_file = convert2annovar(workflow_config['convertAnnovar']['command'], tools_dir, filtered_vcf, work_dir)
        #16. Annotate annnovar file using Annovar
        annotation_file = annotate(workflow_config['annotate']['command'], tools_dir, annovar_file, work_dir)
        #17. Summarize annovar file using Annovar
        summary_file = summarize(workflow_config['summarize']['command'], tools_dir, annovar_file,
                                        '1000g2012apr', '6500','137', 'refgene', 'hg19', work_dir)
        '''

def parse_cl_args():
    '''Parse input commandline arguments, handling multiple cases.

    Returns the main config file and set of kwargs
    '''

    parser = argparse.ArgumentParser(description =
            'Simple fully automated sequencing analysis')

    parser.add_argument('global_config', help='Global YAML configuration file specifying details '
                'about the system')
    parser.add_argument('fc_dir', help='A directory of fastq files to process (optional)',
                nargs = "?")
    parser.add_argument('workflow', help='YAML file with details about the pipeline workflow', nargs='?')
    parser.add_argument('--reference', help="Human genome Reference to use as base.",
                    default = 'hg19')
    parser.add_argument('--dbsnp', help="dbSnp Reference to use as base.")
    parser.add_argument('--workdir', help="Directory to process in. Defaults to current working directory",
                    default = os.getcwd())
    parser.add_argument('--tooldir', help="Directory where the tools are in. Defaults to current working directory",
                    default = os.getcwd())
    parser.add_argument('--bamconverter', help='Tool to convert sam 2 bam.', choices=['samtools', 'picard'], default = 'samtools')
    parser.add_argument('--logdir', help="Directory where the log files will be stored. Defaults to current working directory",
                    default = os.getcwd())
    return parser

def make_output_dir(dir):
    if not os.path.exists(dir):
        print('Creating folder %s' % dir)
        try:
            os.mkdir(dir, 0777)
        except IOError, e:
            raise IOError('%s\nFailed to make the directory %s' (e, dir))
    return dir

def main(args, sys_args , parser):
    #read the global config and extract info.
    if os.path.exists(args.global_config):
        with open(args.global_config) as f:
            contents = f.read()
            newConfig = load(contents)
    else:
        raise IOError('GLobal YAML config file not found')

    #Check the input directory
    if args.fc_dir:
        fc_dir = os.path.abspath(args.fc_dir)
    else:
        fc_dir = os.getcwd()

    #Create the output directory
    work_dir = os.path.abspath(args.workdir)
    make_output_dir(work_dir)

    is_picard_to_bam =  False if args.bamconverter == 'samtools' else True

    if args.workflow:
        with open(args.workflow) as f:
            contents = f.read()
            workflowConfig = load(contents)
    else:
        with open('workflow.yaml') as f:
            contents = f.read()
            workflowConfig = load(contents)

    if not args.dbsnp:
        raise IOError('dbsnp reference file not found')

    #check where tools dir is.
    tools_dir = os.path.abspath(args.tooldir)

    #start the logger
    make_output_dir(args.logdir)
    logger = get_logger(os.path.join(args.logdir, 'pipeline.log'))

    run(newConfig['resources'], fc_dir, work_dir, tools_dir, workflowConfig['stages']['algorithm'], args.reference,
                            os.path.abspath(args.dbsnp), is_picard_to_bam)

if __name__ == '__main__':
    parser = parse_cl_args()
    main(parser.parse_args(), sys.argv[1:], parser)
