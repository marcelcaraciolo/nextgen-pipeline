import os
import subprocess
import re
import sys

def runCommand(message, command):
    subprocess.check_call(command, shell=True)

def splitPath(path):
    (prefix, base) = os.path.split(path)
    (name, ext) = os.path.splitext(base)
    return (prefix, name, ext)

def make_sample_dir(dir):
    if not os.path.exists(dir):
        print('Creating folder %s' % dir)
        try:
            os.mkdir(dir, 0777)
        except IOError, e:
            raise IOError('%s\nFailed to make the directory %s' (e, dir))

def parse_and_link(filename, out_dir, metadata):
    '''
    Parse metadata out of input filenames and construct the symbolic links.
    It takes the fastq filename, destination directory and metadata info (dict).
    Parse the filename to extract the sample name, run, read, etc.
    '''

    fastq_new_pattern = re.match(r".*?/([a-zA-Z0-9-.]+)_S([0-9]+)_L([0-9]+)_R(1|2)_([0-9]+).fastq.gz",filename)
    if fastq_new_pattern:
        run_id = fastq_new_pattern.group(5)
        sample = fastq_new_pattern.group(2)
        lane = (fastq_new_pattern.group(3))
        pair = fastq_new_pattern.group(4)
        identifier = fastq_new_pattern.group(1)
    else:
        sys.exit('Unable to parse the name of the fastq.gz file %s' % filename)
    new_fastq = '%s_%s_L%s_%s.fastq.gz' % (identifier +'_'+ sample, run_id, lane, pair)
    out_dir = out_dir + '/' + '_'.join([identifier, sample])
    make_sample_dir(out_dir)
    subprocess.call(['cp','-f', filename, os.path.join(out_dir, new_fastq)])
    new_fastq = os.path.join(out_dir, new_fastq)
    basename =  os.path.splitext(os.path.basename(new_fastq))[0]
    metadata[basename]['sample'] = sample
    metadata[basename]['run_id'] = run_id
    metadata[basename]['lane'] = lane
    metadata[basename]['pair'] = pair
    metadata[basename]['identifier'] = identifier
    metadata[basename]['out_dir'] = out_dir

    return new_fastq



