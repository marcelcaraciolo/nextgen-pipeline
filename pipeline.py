#-*- coding: utf-8 -*-

'''
GATK-based variant calling pipeline

It implements a workflow pipeline for next generation
sequencing variant detection using the Broad Institute's GATK
for variant calling.

The pipeline is configured by an options file in a python file
including the actual command which are run at each stage.

'''
from commands import make_reference_database, index_reference, align
from commands import align2sam, sam2bam, dedup, realign_intervals
from commands import realign, fix_mate, base_qual_recal_count
from commands import base_qual_recal_tabulate, call_snps, filter_snps
from commands import convert2annovar, annotate, summarize
import sys
import os
import argparse
from yaml import load
import glob

def check_fasta_files(fc_dir):
    fasta_files = glob.glob('%s/*.fastq' % fc_dir)
    if len(fasta_files) <= 1:
        if len(fasta_files) == 0:
            exit('At least one sequence file must be specified')
        if 'hg19' in fasta_files[0]:
            exit('At least one sequence file must be specified')

    return fasta_files


def run(global_config, fc_dir, work_dir, tools_dir, workflow_config, reference, dbsnp):
    #1. Get all fasta files and check if there is at least one to process.
    sequence_files = check_fasta_files(fc_dir)

    #2. Create reference database
    make_reference_database(workflow_config['indexer']['command'],'bwtsw', reference)

    #3.Index the reference
    index_reference(reference)

    for seq in sequence_files:
        '''
        #4.Align sequence to the reference database.
        seq_align = align(workflow_config['aligner']['command'], global_config['bwa']['threads'],
                    reference, seq, work_dir)
        #5. Convert alignment to SAM format.
        seq_sam = align2sam(workflow_config['samse']['command'], reference, seq_align, seq, work_dir)
        #6. Convert SAM to BAM
        #@TODO: make picard-tools directory without version to normalize for any releases.
        seq_bam = sam2bam(workflow_config['sam2bam']['command'], global_config['picard']['jvm_opts'],
                    os.path.join(tools_dir, 'picard-tools-1.109'), seq_sam, work_dir)
        #7. Mark PCR Duplicates
        marked_seq_bam = dedup(workflow_config['markduplicates']['command'], global_config['picard']['jvm_opts'],
                    os.path.join(tools_dir, 'picard-tools-1.109'), seq_bam, work_dir)
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

        '''
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
    return parser

def make_output_dir(dir):
    if not os.path.exists(dir):
        print('Creating folder %s' % dir)
        try:
            os.mkdir(dir, 0777)
        except IOError, e:
            raise IOError('%s\nFailed to make the directory %s' (e, dir))

def main(args, sys_args, parser):
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

    run(newConfig['resources'], fc_dir, work_dir, tools_dir, workflowConfig['stages']['algorithm'], args.reference,
                        os.path.abspath(args.dbsnp))

if __name__ == '__main__':
    parser = parse_cl_args()
    main(parser.parse_args(), sys.argv[1:], parser)
