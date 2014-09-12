#!/usr/bin/env python
"""
Created on Mon 9 Sep 2013

Script used for making project 2 in computational physics.

author: Benedicte Emilie Braekken
"""
import os,sys,time

class Project(object):
    """
    Contains the metainformation from the project.
    """
    def __init__(self, metafile):
        """
        @param metafile The file containing metainformation for the project.
        """
        self.data = {}
        self._readData(metafile)

    def __call__(self):
        """
        Can get hold of dictionairy by calling. This results in shorter code
        usage.
        """
        return self.data

    def _readData(self, metafile):
        """
        Stores the data in the instance. Private method.

        @param metafile File containing the data.
        """
        reader = open(metafile,'r')
        reader.readline() # Skip header
        for line in reader:
            if len(line.strip().split(':')) == 2:
                if line.strip().split(':')[0] in ['files','compileFlags','linkFlags',
                        'libLocations','libs','headerSearchDirs']:
                    # Some params are lists
                    identity = line.strip().split(':')[0]
                    value = line.strip().split(':')[1]
                    value = [i.strip() for i in value.strip('[').strip(']').split(',')]
                    # If only contains empty string, make it empty
                    if value[0] == '' and value[-1] == '':
                        value = []
                else:
                    identity,value = line.strip().split(':')
                self.data[identity] = value
            else:
                print 'Found line longer or shorter than one separation.'
                print 'Please correct metafile.'
                reader.close()
                sys.exit(1)

        reader.close()

# Clear screen first
os.system('clear')

# Check that metafile is present
for contents in os.listdir('.'):
    if contents == 'meta.dat':
        break;
    if contents == os.listdir('.')[-1] and contents != 'meta.dat':
        print 'No meta.dat file found in root of project (or here).'
        sys.exit(1)
# Load metafile
p = Project('meta.dat')

start = time.time()

# Create link and compile flagstrings
cmpFlagString = ''; lnkFlagString = ''
for flag in p()['compileFlags']:
    cmpFlagString += '-%s' % flag
    if not flag == p()['compileFlags'][-1]:
        cmpFlagString += ' '
for path in p()['headerSearchDirs']:
    cmpFlagString += '-I%s' % path
    if not path == p()['headerSearchDirs'][-1]:
        cmpFlagString += ' '
for flag in p()['linkFlags']:
    lnkFlagString += '-%s' % flag
    if not flag == p()['linkFlags'][-1]:
        lnkFlagString += ' '

# Check if executable already present
for existing_file in os.listdir(p()['buildDir']):
    if existing_file == '%s.x' % p()['exe']:
        # If it is, delete it
        os.remove('%s%s.x' % (p()['buildDir'],p()['exe']))
        print 'Deleted old executable.\n'

print 'Compiling files:'

# Compile
for cfile in p()['files']:
    name = cfile.split('.')[0]

    if cmpFlagString != '':
        runstring = 'c++ %s -c %s%s -o %s%s.o' % (cmpFlagString,p()['sourceDir'],
                                                cfile,p()['buildDir'],name)
    else:
        runstring = 'c++ -c %s%s -o %s%s.o' % (p()['sourceDir'],cfile,p()['buildDir'],name)

    print runstring, '\n'
    os.system(runstring)

# Link
if len(p()['linkFlags']) != 0:
    linkstring = 'c++ %s -o %s%s ' % (lnkFlagString,p()['buildDir'],p()['exe'])
else:
    linkstring = 'c++ -o %s%s ' % (p()['buildDir'],p()['exe'])

linkstring += '%s*.o ' % p()['buildDir']

if len(p()['libLocations']) != 0:
    for libPath in p()['libLocations']:
        linkstring += '-L%s ' % libPath

if len(p()['libs']) != 0:
    for lib in p()['libs']:
        if lib == p()['libs'][-1]:
            linkstring += '-l%s' % lib
            continue
        linkstring += '-l%s ' % lib

print 'Now linking:\n',linkstring,'\n'
os.system(linkstring)
print ''

# Delete *.o files
for cfile in p()['files']:
    name = cfile.split('.')[0]
    rmstring = '%s%s.o' % (p()['buildDir'],name)
    print 'Removing', rmstring
    os.remove(rmstring)

end = time.time()

print '\nSucessfully compiled in %g s.' % (end - start)
