#!/usr/bin/python

import os, glob
from fabricate import *

flags = ['-std=c++11', '-I../src', '-lGEL', '-c', '-fPIC']
build_dir = 'build'
target = 'libPyGEL.so'
dirs = ['../src/PyGEL']
sources = []

for dir in dirs:
    for file in glob.glob(dir + '/*.cpp'):
        base_file_name, ext = os.path.splitext(file)
        sources.append(base_file_name)

def build():
    if not os.path.exists(build_dir):
        os.mkdir(build_dir)
    compile()
    link()

def oname(build_dir, filename):
    return os.path.join(build_dir, os.path.basename(filename))

def compile():
    for source in sources:
        run('gcc', flags, source+'.cpp', '-o', oname(build_dir, source+'.o'))

def link():
    objects = [oname(build_dir, s+'.o') for s in sources]
    after()
    run('gcc', '-shared', '-o', oname(build_dir, target), objects)

def clean():
    autoclean()

main(parallel_ok=True, jobs=10)
