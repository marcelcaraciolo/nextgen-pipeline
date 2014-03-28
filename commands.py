'''
All pipeline commands for NextGen Sequencing Pipeline (NGS)

'''
import sys
import os
from utils import runCommand, splitPath


def make_metadata_string(metadata):
    return r'-r "@RG\tID:%s\tLB:%s\tSM:%s\tPL:%s"' % (metadata['ID'], metadata['LB'], metadata['SM'],
                    metadata['PL'])

def make_reference_database(command, algorithm, reference):
    '''
    Create reference database
    '''
    ref_database = reference + '.bwt'
    if os.path.exists(ref_database):
        print('Reference database already exists: using %s' % ref_database)
    else:
        command = command % {'algorithm': algorithm, 'prefix': reference, 'seq': reference + '.fasta'}
        runCommand('Creating Reference Database', command)

    return ref_database

def index_reference(reference):
    reference_index = reference + '.fasta.fai'
    if os.path.exists(reference_index):
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
    alignment_file = os.path.join(output_dir, name + '.sai')
    command = command % {'out': alignment_file, 'threads': int(threads), 'ref': reference,
                    'seq': sequence, 'encodingflag': ''}
    runCommand('Running Alignment', command)

    return alignment_file

def align2sam(command, reference, alignment, sequence, output_dir):
    """
    Convert alignments to SAM format. Turn bwa sai alignments into a sam file
    """
    (path, name, ext) = splitPath(alignment)
    if ext != '.sai':
        sys.exit('align2Sam: alignment file %s does not have .sai extension' % alignment)
    sam_file = os.path.join(output_dir, name + '.sam')
    readgroup_metadata = {'PL': 'ILLUMINA', 'SM': name, 'LB': name, 'ID': name}
    metadata_str = make_metadata_string(readgroup_metadata)
    command =  command % {'out': sam_file, 'ref': reference, 'align': alignment,
                                      'seq': sequence, 'meta': metadata_str}
    runCommand('Align to Sam', command)

    return sam_file

def sam2bam(command, command_options, piccard_dir, alignment, output_dir):
    """
    Convert sam to bam and sort, using Picard.
    """
    (path, name, ext) = splitPath(alignment)
    command_options = command_options[1]
    if ext != '.sam':
        sys.exit('sam2bam: alignment file %s does not have .sam extension' % alignment)
    bam_file = os.path.join(output_dir, name + '.bam')
    command = command % {'out': bam_file, 'sam': alignment, 'jvmoptions': command_options,
            'picarddir': piccard_dir}
    runCommand('Sam to Sorted Bam', command)

    return bam_file

def dedup(command, command_options, piccard_dir, alignment, output_dir):
    """
    Remove apparent duplicates using Picard MarkDuplicates
    """
    (path, name, ext) = splitPath(alignment)
    command_options = command_options[1]
    if ext != '.bam':
        sys.exit('mark pcr duplicates: alignment file %s does not have .bam extension' % alignment)
    marked_bam_file = os.path.join(output_dir, name + '.marked.bam')
    command = command % {'out': marked_bam_file, 'bam': alignment, 'jvmoptions': command_options, \
        'picarddir': piccard_dir, 'log': 'metrics'}
    runCommand('Marking PCR duplicates', command)

    return marked_bam_file

def realign_intervals(command, command_options, gatk_dir, reference, alignment, output_dir):
    """
    Run GATK RealignTargetCreator to find suspect intervals for realignment.
    """
    (path, name, ext) = splitPath(alignment)
    command_options = command_options[1]
    if not alignment.endswith('marked.bam'):
        sys.exit('calculating realignment intervals: alignment file %s does not have .bam extension' % alignment)
    interval_file = os.path.join(output_dir, name + '.bam.list')
    command = command % {'out': interval_file, 'bam': alignment, 'jvmoptions': command_options, \
        'gatkdir': gatk_dir, 'ref': reference + '.fasta'}
    runCommand('Calculating realignment intervals', command)

    return interval_file

def realign(command, command_options, gatk_dir, reference, alignment, intervals, output_dir):
    '''
    Run GATK IndelRealigner for local realignment, using intervals found by realign_intervals
    '''
    (path, name, ext) =  splitPath(alignment)
    command_options = command_options[1]
    if not intervals.endswith('bam.list') or ext != '.bam':
        sys.exit('local realignment with intervals: intervals file %s does not have .list extension' % alignment)
    realigned_bam = os.path.join(output_dir, name + '.realigned.bam')
    command = command % {'jvmoptions': command_options, 'ref': reference + '.fasta', 'out': realigned_bam,
                            'bam': alignment, 'gatkdir': gatk_dir, 'intervals': intervals}
    runCommand('Running local realignment around indels', command)

    return realigned_bam

def fix_mate(command, command_options, piccard_dir, alignment, output_dir):
    '''
    Fix mate information in paired end data using picard
    '''
    (path, name, ext) =  splitPath(alignment)
    command_options = command_options[1]
    if ext != '.bam':
        sys.exit('mate information fix: alignment file %s does not have .bam extension' % alignment)
    fixed_bam = os.path.join(output_dir, name + '.fixed.bam')
    command = command % {'jvmoptions': command_options, 'out': fixed_bam,
                            'bam': alignment, 'picarddir': piccard_dir}
    runCommand('Fixing Mate information', command)

    return fixed_bam

def base_qual_recal_count(command, command_options, gatk_dir, reference, dbsnp, alignment, output_dir):
    '''
    GATK CountCovariates, first step of base quality score recalibration.
    '''
    (path, name, ext) =  splitPath(alignment)
    command_options = command_options[1]
    if ext != '.bam':
        sys.exit('count covariates: alignment file %s does not have .bam extension' % alignment)
    recal_file = os.path.join(output_dir, name + '.recal_data.csv')
    command = command % {'jvmoptions': command_options, 'out': recal_file, 'dbsnp': dbsnp,
                            'bam': alignment, 'gatkdir': gatk_dir, 'ref': reference + '.fasta'}
    runCommand('count covariates for base quality score', command)

    return recal_file

def base_qual_recal_tabulate(command, command_options, gatk_dir, reference, recal_file, alignment, output_dir):
    '''
    GATK TableRecalibration: recalibrate base quality scores using the output of CountCovariates.
    '''
    (path, name, ext) =  splitPath(alignment)
    command_options = command_options[1]
    if ext != '.bam':
        sys.exit('table recalibration: alignment file %s does not have .bam extension' % alignment)
    recal_bam = os.path.join(output_dir, name + '.recal.bam')
    command = command % {'jvmoptions': command_options, 'out': recal_bam, 'recalfile': recal_file,
                            'bam': alignment, 'gatkdir': gatk_dir, 'ref': reference + '.fasta'}
    runCommand('recalibrate base quality scores', command)

    return recal_bam

def call_snps(command, command_options, threads, gatk_dir, reference, dbsnp, standard_emit_conf,
                standard_call_conf, dcov, alleles, alignment, output_dir):
    """
    Use GATK HaplotypeGenotyper to call SNPs from recalibrated bams.
    """
    (path, name, ext) =  splitPath(alignment)
    command_options = command_options[1]
    if ext != '.bam':
        sys.exit('call snp : alignment file %s does not have .bam extension' % alignment)
    out_vcf = os.path.join(output_dir, name + '.vcf')
    command = command % {'jvmoptions': command_options, 'out': out_vcf, 'dbsnp': dbsnp, 'alleles': alleles,
                            'threads': threads, 'scf': standard_call_conf, 'sec':standard_emit_conf,
                            'dcov': dcov, 'bam': alignment, 'gatkdir': gatk_dir, 'ref': reference + '.fasta'}
    runCommand('Calling snps', command)

    return out_vcf

def filter_snps(command, command_options, gatk_dir, reference, vcf, filter_expression, output_dir):
    '''
    Use GATK VariantFiltration to filter raw SNP calls.
    '''
    (path, name, ext) =  splitPath(vcf)
    command_options = command_options[1]
    if ext != '.vcf':
        sys.exit('filtering SNPs: vcf file %s does not have .vcf extension' % vcf)
    out_vcf = os.path.join(output_dir, name + '.filtered.vcf')
    command = command % {'jvmoptions': command_options, 'out': out_vcf,
                    'vcf': vcf, 'gatkdir': gatk_dir, 'ref': reference + '.fasta', 'expression': filter_expression}
    runCommand('Calling snps', command)

    return out_vcf

