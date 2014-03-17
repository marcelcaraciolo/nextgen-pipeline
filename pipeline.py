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
import sys
import os
import argparse
from yaml import load
import glob

def check_fasta_files(fc_dir):
    fasta_files = glob.glob('%s/*.fasta' % fc_dir)
    if len(fasta_files) <= 1:
        if len(fasta_files) == 0:
            exit('At least one sequence file must be specified')
        if 'hg19' in fasta_files[0]:
            exit('At least one sequence file must be specified')

    return fasta_files

def run(global_config, fc_dir, work_dir, workflow_config, reference):
    #1. Get all fasta files and check if there is at least one to process.
    #sequence_files = check_fasta_files(fc_dir)
    sequence_files = glob.glob('%s/*.fasta' % fc_dir)

    #2. Create reference database
    make_reference_database(workflow_config['indexer']['command'],'bwtsw', reference)

    #3.Index the reference
    index_reference(reference)

    for seq in sequence_files:
        seq_align = align(workflow_config['aligner']['command'], global_config['bwa']['threads'],
                    reference, seq, work_dir)


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
    parser.add_argument('--workdir', help="Directory to process in. Defaults to current working directory",
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


    run(newConfig['resources'], fc_dir, work_dir, workflowConfig['stages']['algorithm'], args.reference)

if __name__ == '__main__':
    parser = parse_cl_args()
    main(parser.parse_args(), sys.argv[1:], parser)
