##################################################
# Copyright (C) 2012 by 
# Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
# ... and other members of the Reconstruction Team of David Haussler's 
# lab (BME Dept. UCSC).
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE. 
##################################################
import os
import platform
import random
import shutil
import string
import subprocess
import xml.etree.ElementTree as ET
import xml.parsers.expat

def makeTempDir(name=None):
    """
    make the typical directory where all temporary test files will be stored.
    """
    if not os.path.exists(os.path.join(os.curdir, 'tempTestDir')):
        os.mkdir(os.path.join(os.curdir, 'tempTestDir'))
    charSet = string.ascii_lowercase + '123456789'
    if name is None:
        while True:
            name = '%s_%s' % (''.join(random.choice(charSet) for x in xrange(4)), 
                              ''.join(random.choice(charSet) for x in xrange(4)))
            if not os.path.exists(os.path.join(os.curdir, 'tempTestDir', name)):
                break
    if not os.path.exists(os.path.join(os.curdir, 'tempTestDir', name)):
        os.mkdir(os.path.join(os.curdir, 'tempTestDir', name))
    return os.path.join(os.curdir, 'tempTestDir', name)
def removeDir(dirpath):
    """
    destroy a directory
    """
    if os.path.exists(dirpath):
        shutil.rmtree(dirpath)
def runCommandsS(cmds, localTempDir, inPipes=[], outPipes=[]):
    """ 
    runCommandsS uses the subprocess module
    to issue serial processes from the cmds list.
    """
    if not len(inPipes):
        inPipes = [None] * len(cmds)
    if not len(outPipes):
        outPipes = [None] * len(cmds)
    for i, c in enumerate(cmds, 0):
        if inPipes[i] is None:
            sin = None
        else:
            sin = subprocess.PIPE
        if outPipes[i] is None:
            sout = None
        else:
            sout = subprocess.PIPE
        p = subprocess.Popen(c, cwd=localTempDir, stdin=sin, stdout=sout)
            
        if inPipes[i] is None:
            sin = None
        else:
            if not os.path.exists(inPipes[i]):
                raise IOError('Unable to locate inPipe file: %s for command %s' % (inPipes[i], ' '.join(c)))
            sin = open(inPipes[i], 'r').read()
        if outPipes[i] is None:
            pout, perr = p.communicate(sin)
            handleReturnCode(p.returncode, cmds[i])
        else:
            f = open(outPipes[i], 'w')
            f.write(p.communicate(sin)[0])
            f.close()
            handleReturnCode(p.returncode, cmds[i])
def handleReturnCode(retcode, cmd):
    """
    handle the return codes from runCommandsS
    """
    if not isinstance(retcode, int):
        raise TypeError('handleReturnCode takes an integer for '
                        'retcode, not a %s.' % retcode.__class__)
    if not isinstance(cmd, list):
        raise TypeError('handleReturnCode takes a list for '
                        'cmd, not a %s.' % cmd.__class__)
    if retcode:
        if retcode < 0:
            raise RuntimeError('Experienced an error while trying to execute: '
                               '%s SIGNAL:%d' %(' '.join(cmd), -retcode))
        else:
            raise RuntimeError('Experienced an error while trying to execute: '
                               '%s retcode:%d' %(' '.join(cmd), retcode))
def which(program):
    """which() acts like the unix utility which, but is portable between os.
    If the program does not exist in the PATH then 'None' is returned. 
    """
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath != '':
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None
def noMemoryErrors(xml):
    """
    parse the valgrind output xml file looking for any signs of memory errors.
    returns False if memory errors exist, True if no memory errors.
    """
    if not os.path.exists(xml):
        raise RuntimeError('Input xml, %s does not exist.' % xml)
    try:
        tree = ET.parse(xml)
    except xml.parsers.expat.ExpatError:
        raise RuntimeError('Input xml, %s is not a well formed xml document.' % xml)
    root = tree.getroot()
    errors = root.findall('error')
    if len(errors):
        return False
    errorcounts = root.find('errorcounts')
    if len(errorcounts):
        return False
    return True
def mafIsEmpty(maf, headers):
    """ take a path to a maf file and a list of possible header strings.
    an "empty" maf in this context will contain only a header line.
    """
    f = open(maf)
    s = f.read()
    for h in headers:
        if s == h:
            return True
    return False
def fileIsEmpty(filename):
    """
    given a filename, if the file is empty returns True, returns False otherwise
    """
    f = open(filename)
    s = f.read()
    if s == '':
        return True
    return False
def processHeader(f):
    """
    takes in a file handle, f, and reads off lines until it comes to the first
    alignment block line, or the first blank line. If blank it returns nothing,
    if an alignment line it returns that line.
    """
    for line in f:
        line = line.strip()
        if line == '':
            return None
        if line.startswith('a'):
            return line
def extractBlockStr(f, lastLine=None):
    """
    read one block from a maf file handle and turn it into a single string
    """
    if lastLine is None:
        block = ''
    else:
        block = lastLine + '\n'
    for line in f:
        line = line.strip()
        if line == '':
            return block + '\n'
        block += '%s\n' % line
    return block + '\n'
def testFile(mafFile, s, headers):
    """
    given a path to a maffile, a string containing the desired contents of 
    the file and a list of different possible headers, pick one header at
    random and write the string to the file.
    """
    f = open(mafFile, 'w')
    header = random.choice(headers)
    f.write(header)
    f.write(s)
    f.close()
    return mafFile, header
def genericValgrind(tmpDir):
    """ 
    returns a list (in the subprocess command style) containing
    a generic call to valgrind.
    """
    valgrind = which('valgrind')
    if platform.mac_ver() == ('', ('', '', ''), ''):
        return [valgrind, '--leak-check=full', '--show-reachable=yes', '--track-origins=yes', 
                '--xml=yes', '--xml-file=' + os.path.join(tmpDir, 'valgrind.xml')]
    else:
        # --dsymutil is for mac os x builds
        return [valgrind, '--leak-check=full', '--show-reachable=yes', '--track-origins=yes', 
                '--dsymutil=yes', '--xml=yes', '--xml-file=' + os.path.join(tmpDir, 'valgrind.xml')]
