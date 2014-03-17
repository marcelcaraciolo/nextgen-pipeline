'''
All pipeline commands for NextGen Sequencing Pipeline (NGS)

'''
import sys
import os.path as path
from utils import runCommand, splitPath


def make_reference_database(command, algorithm, reference):
    '''
    Create reference database
    '''
    ref_database = reference + '.bwt'
    if path.exists(ref_database):
        print('Reference database already exists: using %s' % ref_database)
    else:
        command = command % {'algorithm': algorithm, 'prefix': reference, 'seq': reference + '.fasta'}
        runCommand('Creating Reference Database', command)

    return ref_database

def index_reference(reference):
    reference_index = reference + '.fai'
    if path.exists(reference_index):
        print('Reflist already exists: using %s' % reference_index)
    else:
       sys.exit('check if the bwa index ran successfully')


def align(command, threads, reference, sequence, output_dir):
    '''
    Align sequence reads to the reference genome. This is the bwa's first stage, bwa aln.
    Use -I for sequence files.
    '''
    (path, name, ext) = splitPath(sequence)
    if ext != '.fastq':
        sys.exit('align: sequence file %s does not have .fastq extension' % sequence)
    alignment_file = path.join(output_dir, name + '.sai')
    command = command % {'out': alignment_file, 'threads': int(threads), 'ref': reference,
                    'seq': sequence}
    runCommand('Running Alignment', command)

    return alignment_file

