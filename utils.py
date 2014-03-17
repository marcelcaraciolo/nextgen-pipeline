import os
import subprocess

def runCommand(message, command):
    subprocess.check_call(command, shell=True)

def splitPath(path):
    (prefix, base) = os.path.split(path)
    (name, ext) = os.path.splitext(base)
    return (prefix, name, ext)

