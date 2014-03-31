'''
All pipeline commands for NextGen Sequencing Pipeline (NGS)

'''
import sys
import os
from utils import runCommand, splitPath
import subprocess

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


def fastqc(command, sequences, fastq_metadata, output_dir):
    '''
    Run FastQC on each fastq file.
    '''
    for fastq_file in sequences:
        command = command % {'outdir': output_dir, 'seq': fastq_file}
        runCommand('Checking fastq quality', command)

def can_pipe(command, fastq_file):
    '''
    bwa-mem handles longer (> 70bp) reads with improved piping.
    Randomly samples 5000 reads from the first two million.
    Default to no piping if more than 75% of the sampled reads are small.
    '''
    min_size = 70
    thresh = 0.75
    head_count = 8000000
    tocheck = 5000
    cat_cmd = 'cat {fastq_file}'
    cmd = (cat_cmd + " | head -n {head_count} | "
           "{command} sample -s42 - {tocheck} | "
           "awk '{{if(NR%4==2) print length($1)}}' | sort | uniq -c")
    count_out = subprocess.check_output(cmd.format(**locals()), shell=True,
                                        executable="/bin/bash", stderr=open("/dev/null", "w"))
    if not count_out.strip():
        raise IOError("Failed to check fastq file sizes with: %s" % cmd.format(**locals()))

    shorter = 0
    for count, size in (l.strip().split() for l in count_out.strip().split("\n")):
        if int(size) < min_size:
            shorter += int(count)
    return (float(shorter) / float(tocheck)) <= thresh

def align_with_mem(command, threads, reference, fastq_file, pair_file, fastq_metadata, output_dir):
    '''
    Perform alignment on two paired-end fastq files to a reference genome to produce a sam file.
    '''
    (path, name, ext) = splitPath(fastq_file)
    (pathP, nameP, extP) = splitPath(pair_file)

    if ext != '.fastq' or extP != '.fastq':
        sys.exit('align: one of the fastq file %s or %s does not have .fastq extension' % (fastq_file, pair_file))

    sam_file = os.path.join(output_dir, os.path.splitext(os.path.basename(fastq_file))[0]) +  '.sam'
    sample =  fastq_metadata[os.path.basename(fastq_file)]['sample']
    run_id =  fastq_metadata[os.path.basename(fastq_file)]['run_id']
    lane =   fastq_metadata[os.path.basename(fastq_file)]['lane']
    identifier =  fastq_metadata[os.path.basename(fastq_file)]['identifier']
    readgroup_metadata = {'PL': 'ILLUMINA', 'SM': sample,
                            'LB': '%s_%s_%s_Lane%s' % (identifier, sample, run_id, lane),
                            'ID':  '%s_%s_%s_Lane%s' % (identifier, sample, run_id, lane) }
    metadata_str = make_metadata_string(readgroup_metadata)

    command = command % {'threads': threads, 'meta': metadata_str, 'ref': reference,
                                'seq': fastq_file , 'pair': pair_file, 'out': sam_file}
    runCommand('bwa mem alignment from fastq: %s' % sample, command)

    return sam_file

def align(command, threads, reference, sequence, fastq_metadata, output_dir):
    '''
    Align sequence reads to the reference genome. This is the bwa's first stage, bwa aln.
    '''
    (path, name, ext) = splitPath(sequence)
    if ext != '.fastq':
        sys.exit('align: sequence file %s does not have .fastq extension' % sequence)
    alignment_file = os.path.join(output_dir, name + '.sai')
    command = command % {'out': alignment_file, 'threads': int(threads), 'ref': reference,
                    'seq': sequence, 'encodingflag': ''}
    runCommand('Running Alignment', command)

    return alignment_file

def alignPE2sam(command, reference, fastq_file, pair_file, sai_fastq_file, sai_pair_file,
                        fastq_metadata, output_dir):
    '''
    Convert alignments to SAM format. Turn bwa sai alignments into a sam file.
    It uses bwa sampe commandline. (Pair End only)
    '''
    (path, name, ext) = splitPath(sai_fastq_file)
    (pathP, nameP, extP) = splitPath(sai_pair_file)
    if ext != '.sai' or extP != '.sai':
        sys.exit('alignPE2sam: one .sai file %s or %s does not have .sai extension' % (sai_fastq_file, sai_pair_file))

    sam_file = os.path.join(output_dir, os.path.splitext(os.path.basename(fastq_file))[0]) +  '.sam'
    sample =  fastq_metadata[os.path.basename(fastq_file)]['sample']
    run_id =  fastq_metadata[os.path.basename(fastq_file)]['run_id']
    lane =   fastq_metadata[os.path.basename(fastq_file)]['lane']
    identifier =  fastq_metadata[os.path.basename(fastq_file)]['identifier']
    readgroup_metadata = {'PL': 'ILLUMINA', 'SM': sample,
                            'LB': '%s_%s_%s_Lane%s' % (identifier, sample, run_id, lane),
                            'ID':  '%s_%s_%s_Lane%s' % (identifier, sample, run_id, lane) }
    metadata_str = make_metadata_string(readgroup_metadata)

    command = command % {'meta': metadata_str, 'ref': reference, 'align': sai_fastq_file, 'alignP': sai_pair_file,
                                'seq': fastq_file , 'pair': pair_file, 'out': sam_file}
    runCommand('bwa sampe alignment from fastq: %s' % sample, command)

    return sam_file

def align2sam(command, reference, fastq_file, sai_fastq_file, fastq_metadata, output_dir):
    """
    Convert alignments to SAM format. Turn bwa sai alignments into a sam file.
    It uses bwa samse commandline.
    """
    (path, name, ext) = splitPath(sai_fastq_file)
    if ext != '.sai':
        sys.exit('align2Sam: alignment file %s does not have .sai extension' % sai_fastq_file)

    sam_file = os.path.join(output_dir, os.path.splitext(os.path.basename(fastq_file))[0]) +  '.sam'
    sample =  fastq_metadata[os.path.basename(fastq_file)]['sample']
    run_id =  fastq_metadata[os.path.basename(fastq_file)]['run_id']
    lane =   fastq_metadata[os.path.basename(fastq_file)]['lane']
    identifier =  fastq_metadata[os.path.basename(fastq_file)]['identifier']
    readgroup_metadata = {'PL': 'ILLUMINA', 'SM': sample,
                            'LB': '%s_%s_%s_Lane%s' % (identifier, sample, run_id, lane),
                            'ID':  '%s_%s_%s_Lane%s' % (identifier, sample, run_id, lane) }
    metadata_str = make_metadata_string(readgroup_metadata)

    command =  command % {'out': sam_file, 'ref': reference, 'align': sai_fastq_file,
                                      'seq': fastq_file, 'meta': metadata_str}

    runCommand('bwa samse alignment from fastq: %s' % sample, command)

    return sam_file

def samP2bam(command, command_options, piccard_dir, alignment, output_dir):
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

def samS2bam(command, command_options, piccard_dir, alignment, output_dir):
    """
    Convert sam to bam and sort, using Samtools.
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

def convert2annovar(command, annovar_dir, vcf, output_dir):
    '''
    Convert vcf file to Annovar variant caller .annovar
    '''
    (path, name, ext) =  splitPath(vcf)
    if ext != '.vcf':
        sys.exit('Converting to .annovar: vcf file %s does not have .vcf extension' % vcf)
    out_prefix = os.path.join(output_dir, name.split('.')[0])
    command = command % {'out': out_prefix, 'vcf': vcf, 'annovardir': annovar_dir}
    runCommand('Coverting to .annovar format', command)

    return '.'.join([out_prefix, name.split('.')[0],  'avinput'])

def annotate(command, annovar_dir, annovar_file, output_dir):
    '''
    Annotate vcf using Annovar variant caller.
    '''
    (path, name, ext) =  splitPath(annovar_file)
    if ext != '.avinput':
        sys.exit('Annotating vcf: vcf file %s does not have .avinput extension' % annovar_file)
    out_prefix = os.path.join(output_dir, name)
    command = command % {'out': out_prefix, 'annovarfile': annovar_file, 'annovardir': annovar_dir}
    runCommand('Annotating with Annovar', command)

    return out_prefix

def summarize(command, annovar_dir, annovar_file, ver1000g, veresp, verdbsnp, genetype, buildver, output_dir):
    '''
    Summarize information with Annovar.
    '''
    (path, name, ext) =  splitPath(annovar_file)
    if ext != '.avinput':
        sys.exit('Summarizing annotations: vcf file %s does not have .avinput extension' % annovar_file)
    out_prefix = os.path.join(output_dir, name)
    command = command % {'out': out_prefix, 'annovarfile': annovar_file, 'annovardir': annovar_dir,
                'ver1000g': ver1000g, 'veresp': veresp, 'verdbsnp': verdbsnp, 'genetype': genetype, 'buildver': buildver}
    runCommand('Summarizing with Annovar', command)

    return out_prefix






