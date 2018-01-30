#!/usr/bin/python

import os, glob
from fabricate import *

flags = ['-std=c++11', '-I../src/GEL', '-DNOT_HAVE_SA_LEN', '-c', '-fPIC']
build_dir = 'build'
target = 'libGEL.so'
dirs = ['../src/GEL/CGLA', '../src/GEL/GLGraphics','../src/GEL/Geometry','../src/GEL/HMesh','../src/GEL/Util']
sources = []
for dir in dirs:
    for file in glob.glob(dir + '/*.c*'):
        sources.append(file)

def build():
    if not os.path.exists(build_dir):
        os.mkdir(build_dir)
    compile()
    link()

def oname(build_dir, filename):
    return os.path.join(build_dir, os.path.basename(filename))

def compile():
    for source in sources:
        base, ext = os.path.splitext(source)
        run('gcc', flags, source, '-o', oname(build_dir, base+'.o'))

def link():
    objects = []
    for source in sources:
        base, ext = os.path.splitext(source)
        objects.append(oname(build_dir, base+'.o'))
    after()
    run('gcc', '-shared', '-o', oname(build_dir, target), objects)

def clean():
    autoclean()

main(parallel_ok=True, jobs=10)
