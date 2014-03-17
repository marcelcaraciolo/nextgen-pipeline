#-*- coding: utf-8 -*-

'''
GATK-based variant calling pipeline

It implements a workflow pipeline for next generation
sequencing variant detection using the Broad Institute's GATK
for variant calling.

The pipeline is configured by an options file in a python file
including the actual command which are run at each stage.

'''

import sys
import os
import argparse
from yaml import load
import glob


def run(global_config, fc_dir, work_dir, workflow_config):
    #1. Get all fasta files and check if there is at least one to process.
    fasta_files = glob.glob('%s/*.fasta' % fc_dir)
    if len(fasta_files) <= 1:
        if len(fasta_files) == 0:
            exit('At least one sequence file must be specified')
        if 'hg19' in fasta_files[0]:
            exit('At least one sequence file must be specified')

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


    run(newConfig, fc_dir, work_dir, workflowConfig)

if __name__ == '__main__':
    parser = parse_cl_args()
    main(parser.parse_args(), sys.argv[1:], parser)
