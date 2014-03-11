#-*- coding: utf-8 -*-


'''
Automatically install required tools and data to run human variant pipelines.

This automates the steps required for installation and setup to make it
easier to get started with next-gen pipeline.

Requires: git, Python 2.7 or argparse for earlier versions.

@author marcel@genomika.com.br

'''

import os
import sys
import subprocess
import contextlib
import shutil
import platform
import datetime

remotes = {"anaconda": "http://repo.continuum.io/miniconda/Miniconda-3.0.0-%s-x86_64.sh",
            'system_config':  'nextgen_system.yaml',
            'requirements': 'https://raw.github.com/genomika/nextgen-pipeline/master/requirements.txt'

}


def main(args, sys_argv):
    check_dependencies()
    with nextgen_tmpdir():
        setup_data_dir(args)
        print("Installing isolated base python installation")
        anaconda = install_anaconda_python(args, remotes)
        print("Installing nextgen-pipeline")
        install_conda_packages(anaconda)
        nextgen = bootstrap_nextgen(anaconda, args, remotes)

    print("Installing data and third party dependencies")
    install_data_tools(sys_argv, args, nextgen)
    system_config = write_system_config(remotes['system_config'], args.datadir, args.tooldir)

    print("Finished: nextgen-pipeline, tools and data installed")
    print(" Genome data installed in: \n %s" % args.datadir)
    if args.tooldir:
        print(" Tools installed in:\n %s" %  args.tooldir)
    print("  Ready to use system configuration at:\n %s" % system_config)
    print("Edit configuration file as need to match  your machine")


def install_data_tools(sys_argv, args, nextgen):
    '''
    Handle installation and updates of nextgen, third party software and data.
    Automated installation tool and in-place updates to install additional data
    and software.
    '''
    base = [ x for x in sys_argv if x.startswith('-') or not
                args.datadir == os.path.abspath(os.path.expanduser(x))]

    print('Installing nextgen tools')
    #install homebrew, BWA, Picard, SamTools, GATK


    print('Installing nextgen data files')
    install_nextgen_data(args, remotes)


def install_nextgen_data(args, REMOTES):
    pass

def write_system_config(base_url, datadir, tooldir):
    """
    Write a nextgen_system.yaml configuration file with tool information.
    """
    out_file = os.path.join(datadir, "galaxy", os.path.basename(base_url))
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    if os.path.exists(out_file):
        #if not tool directory and exists, do not overwrite
        if tooldir is None:
            return out_file
        else:
            bak_file = out_file + '.bak%s' % (datetime.datetime.now().strftime('%Y%M%d_%H%M'))
            shutil.copy(out_file, bak_file)

    if tooldir:
        java_basedir = os.path.join(tooldir, "share", "java")

    rewrite_ignore = ("log",)
    #WHY NOT USE @YAML?! REFACTOR!
    with contextlib.closing(open(base_url)) as in_handler:
        with open(out_file, 'w') as out_handler:
            in_resources = False
            in_prog = None
            for line in in_handler:
                if  line[0] != ' ':
                    in_resources = line.startswith('resources')
                    in_prog = None
                elif(in_resources and line[:2] == '  ' and line[2]  != " " \
                            and not line.strip().startswith(rewrite_ignore)):
                    in_prog = line.split(':')[0].strip()
                elif line.strip().startswith('dir:') and in_prog:
                    final_dir = os.path.basename(line.split()[-1])
                    if tooldir:
                        line = '%s: %s\n' % (line.split(':')[0], \
                                os.path.join(java_basedir, final_dir))
                    in_prog = None
                out_handler.write(line)

    return out_file



def bootstrap_nextgen(anaconda, args, remotes):
    """Install nextgen to bootstrap rest of installation process.
    """
    subprocess.check_call([anaconda["pip"], "install", "fabric"])
    subprocess.check_call([anaconda['pip'], 'install', '-r', remotes['requirements']])

    out = {}
    for script in ["pipeline.py"]:
        ve_script = os.path.join(anaconda["dir"], "bin", script)
        if args.tooldir:
            final_script = os.path.join(args.tooldir, "bin", script)
            sudo_cmd = ["sudo"]
            subprocess.check_call(sudo_cmd + ['mkdir', '-p', os.path.dirname(final_script)])
            if os.path.lexists(final_script):
                cmd = ['rm', '-f', final_script]
                subprocess.check_call(sudo_cmd + cmd)
            cmd = ["ln", '-s', ve_script, final_script]
            subprocess.check_call(sudo_cmd + cmd)
        out[script] = ve_script

    return out




def install_conda_packages(anaconda):
    packages = ["pip" ]
    subprocess.check_call([anaconda["conda"], 'install', '--yes'] +  packages)

def install_anaconda_python(args, remotes):
    """Provide insolated installation of Anaconda python for running nextgen.
    http://docs.continuum.io/anaconda/index.html
    """
    anaconda_dir = os.path.join(args.datadir, "anaconda")
    bin_dir = os.path.join(anaconda_dir, "bin")
    conda = os.path.join(bin_dir, "conda")

    if not os.path.exists(anaconda_dir) or not os.path.exists(conda):
        if os.path.exists(anaconda_dir):
            shutil.rmtree(anaconda_dir)

        dist = args.distribution if args.distribution else _guess_distribution()
        url = remotes['anaconda'] % ("MacOSX" if dist.lower() == 'macosx' else 'Linux')
        if not os.path.exists(os.path.basename(url)):
            subprocess.check_call(['wget', url])
        subprocess.check_call("bash %s -b -p %s" % (os.path.basename(url), anaconda_dir),
                        shell=True)

    return {"conda": conda,  "pip": os.path.join(bin_dir, "pip"), "dir": anaconda_dir}

def _guess_distribution():
    """
    Approach to identify if the OS active is a  MacOSX or Linux System for Anaconda.
    """
    if platform.mac_ver()[0]:
        return "macosx"
    else:
        return "linux"

@contextlib.contextmanager
def nextgen_tmpdir():
    orig_dir = os.getcwd()
    work_dir  = os.path.join(os.getcwd(), 'tmpnextgen-install')
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    os.chdir(work_dir)
    yield work_dir
    os.chdir(orig_dir)
    shutil.rmtree(work_dir)

def setup_data_dir(args):
    if not os.path.exists(args.datadir):
        cmd = ["mkdir", "-p", args.datadir]
        cmd.insert(0, "sudo")
        subprocess.check_call(cmd)
        subprocess.check_call(["sudo", "chown", "-R", os.environ["USER"], args.datadir])

def check_dependencies():
    """Ensure required tools for installation are present"""
    print("Checking required dependencies")
    try:
        subprocess.check_call(['git', "--version"])
    except OSError:
        raise OSError("nextgen installer requires Git. please install it from"
                        "http://git-scm.com/")

if __name__ == '__main__':
    try:
        import argparse
    except ImportError:
        raise ImportError('nextgen-pipeline requires "argparse", included '
                            'in Python 2.7\n'
                            'Install for earlier versions with pip install argparse'
                            'or easy_install argparse')

    parser = argparse.ArgumentParser(
            description="Automatic installation for next-gen pipelines")

    parser.add_argument("datadir", help='Directory to install genome data',
            type=lambda y: (os.path.abspath(os.path.expanduser(y))))

    parser.add_argument("--tooldir", help="Directory to install 3rd pary software tools. \
                            Leave unspecified for no tools",
                            type=lambda y: (os.path.abspath(os.path.expanduser(y))), default=None)

    parser.add_argument("--genomes", help='Genomes to download',
                            action='append', default=['hg19'], choices=['hg19', 'GRCh37'])

    parser.add_argument("--aligners", help="Aligner indexes to download",
                            action="append", default=["bwa"], choices=["bwa"])

    parser.add_argument("--distribution", help="Operating system distribution",
                            default="", choices=["ubuntu", "macosx"])

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args(), sys.argv[1:])

