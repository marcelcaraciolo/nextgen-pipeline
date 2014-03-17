'''
All pipeline commands for NextGen Sequencing Pipeline (NGS)

'''
import sys
import os.path as path
from utils import runCommand, splitPath


def make_metadata_string(metadata):
    return r'-r"@RG\tID:%s\tLB:%s\tSM:%s\tPL:%s"' % (metadata['ID'], metadata['LB'], metadata['SM'],
                    metadata['PL'])

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

def align2Sam(command, reference, alignment, sequence, output_dir):
    """
    Convert alignments to SAM format. Turn bwa sai alignments into a sam file
    """
    (path, name, ext) = splitPath(alignment)
    if ext != '.sai':
        sys.exit('align2Sam: alignment file %s does not have .sai extension' % alignment)
    sam_file = path.join(output_dir, name + '.sam')
    readgroup_metadata = {'PL': 'ILLUMINA', 'SM': name, 'LB': name, 'ID': name}
    metadata_str = make_metadata_string(readgroup_metadata)
    command =  command % {'out': sam_file, 'ref': reference, 'seq': sequence, 'meta': metadata_str}
    runCommand('Align to Sam', command)

    return sam_file
