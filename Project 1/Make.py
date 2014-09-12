#!/usr/bin/env python
"""
Created on Mon 9 Sep 2013

Script makes the project 1 c++ project.

author: Benedicte Emilie Braekken
"""
import os,sys,time

if len(sys.argv) == 1:
    print 'Give name for output exe.'
    sys.exit(1)

xName = sys.argv[1]
buildPath = '../build/'
binaryPath = 'FYS3150-Prj01/'
libraries = ['armadillo']
fileList = ['GaussianElimination.cpp',
            'ThomasAlgorithm.cpp',
            'relativeError.cpp',
            'main.cpp']

start = time.time()

# Check if executable already present
for existing_file in os.listdir(buildPath):
    if existing_file == '%s.x' % xName:
        # If it is, delete it
        os.remove('%s%s.x' % (buildPath,xName))
        print 'Deleted old executable.'

# Compile
for cfile in fileList:
    name = cfile.split('.')[0]
    os.system('c++ -O3 -c %s%s -o %s%s.o' % (binaryPath,cfile,buildPath,name))

# Link
linkstring = 'c++ -o %s%s.x ' % (buildPath,xName)
for lib in libraries:
    linkstring += '-l%s ' % lib
linkstring += '%s*.o' % buildPath
os.system(linkstring)

# Delete *.o files
for cfile in fileList:
    name = cfile.split('.')[0]
    os.remove('%s%s.o' %(buildPath,name))

end = time.time()

print 'Sucessfully compiled in %g s.' % (end - start)
