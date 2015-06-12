#!/usr/bin/python

import os, glob
from fabricate import *

flags = ['-c', '-std=c++11', '-w', '-I../src/GEL']
build_dir = 'build'
target = 'libGEL.a'
dirs = ['../src/GEL/CGLA', '../src/GEL/GLGraphics','../src/GEL/Geometry','../src/GEL/HMesh','../src/GEL/Util']
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
        run('g++', flags, source+'.cpp', '-o', oname(build_dir, source+'.o'))

def link():
    objects = [oname(build_dir, s+'.o') for s in sources]
    after()
    run('ar', '-cr', oname(build_dir, target), objects)

def clean():
    autoclean()

main(parallel_ok=True, jobs=10)
