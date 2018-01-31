#!/usr/bin/python

import os
import glob
from fabricate import *

flags = ['-std=c++11', '-I../src', '-lGEL', '-fPIC']
build_dir = 'build/PyGEL'
target = 'libPyGEL.so'
dirs = ['../src/PyGEL']
dependencies = ['-lGEL', '-lGL', '-lGLEW', '-lm', '-lglfw', '-lstdc++']
sources = []

for dir in dirs:
    for file in glob.glob(dir + '/*.cpp'):
        base_file_name, ext = os.path.splitext(file)
        sources.append(base_file_name)


def build():
    if not os.path.exists(build_dir):
        os.makedirs(build_dir)
    compile()
    link()


def oname(build_dir, filename):
    return os.path.join(build_dir, os.path.basename(filename))


def compile():
    for source in sources:
        run('gcc', flags, '-c', source+'.cpp', '-o', oname(build_dir, source+'.o'))


def link():
    objects = [oname(build_dir, s+'.o') for s in sources]
    after()
    run('gcc', flags, '-shared', '-Wl,-unresolved-symbols=report-all', objects, '-o', oname(build_dir, target), dependencies)


def clean():
    autoclean()


main(parallel_ok=True, jobs=10)
