#!/usr/bin/env python
"""
Created on Mon 9 Sep 2013

Script used for making project 2 in computational physics.

author: Benedicte Emilie Braekken
"""
import os,sys,time

if len(sys.argv) == 1:
    print 'Give name for output exe.'
    sys.exit(1)

xName = sys.argv[1]
buildPath = 'build/'
binaryPath = 'src/'
compileFlags = []
linkFlags = []
libLocations = []
libraries = ['armadillo']
fileList = ['main.cpp',
            'JacobiRotationProblem.cpp',
            'potentials.cpp']

start = time.time()

# Create link and compile flagstrings
cmpFlagString = ''; lnkFlagString = ''
for flag in compileFlags:
    cmpFlagString += '-%s' % flag
    if not flag == compileFlags[-1]:
        cmpFlagString += ' '
for flag in linkFlags:
    lnkFlagString += '-%s' % flag
    if not flag == linkFlags[-1]:
        lnkFlagString += ' '

# Check if executable already present
for existing_file in os.listdir(buildPath):
    if existing_file == '%s.x' % xName:
        # If it is, delete it
        os.remove('%s%s.x' % (buildPath,xName))
        print 'Deleted old executable.\n'

print 'Compiling files:'

# Compile
for cfile in fileList:
    name = cfile.split('.')[0]

    if len(compileFlags) != 0:
        runstring = 'c++ %s -c %s%s -o %s%s.o' % (cmpFlagString,binaryPath,
                                                cfile,buildPath,name)
    else:
        runstring = 'c++ -c %s%s -o %s%s.o' % (binaryPath,cfile,buildPath,name)

    print runstring, '\n'
    os.system(runstring)

# Link
if len(linkFlags) != 0:
    linkstring = 'c++ %s -o %s%s.x ' % (lnkFlagString,buildPath,xName)
else:
    linkstring = 'c++ -o %s%s.x ' % (buildPath,xName)

if len(libLocations) != 0:
    for libPath in libLocations:
        linkstring += '-L%s ' % libPath

if len(libraries) != 0:
    for lib in libraries:
        linkstring += '-l%s ' % lib

linkstring += '%s*.o' % buildPath
print 'Now linking:\n',linkstring,'\n'
os.system(linkstring)
print ''

# Delete *.o files
for cfile in fileList:
    name = cfile.split('.')[0]
    rmstring = '%s%s.o' % (buildPath,name)
    print 'Removing', rmstring
    os.remove(rmstring)

end = time.time()

print '\nSucessfully compiled in %g s.' % (end - start)
