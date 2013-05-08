#!/usr/bin/env python2.7
import os
import subprocess
import sys
import time
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '../../include/')))
import mafToolsTest as mtt

g_BoilerPlate = '''/* 
 * Copyright (C) 2009-2013 by 
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
 * Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
 * Mark Diekhans (markd@soe.ucsc.edu)
 * ... and other members of the Reconstruction Team of David Haussler's 
 * lab (BME Dept. UCSC).
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE. 
 */
'''
g_git = mtt.which('git')
def runCommand(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    pout, perr = p.communicate()
    mtt.handleReturnCode(p.returncode, cmd)
    return pout
def getBranch():
    branchList = runCommand([g_git, 'branch']).split('\n')
    for b in branchList:
        if b.startswith('* '):
            return b[2:]
def getSha():
    return runCommand([g_git, 'rev-parse', 'HEAD']).strip()
def writeHeader(location):
    f = open(os.path.join(location, 'buildVersion.h'), 'w')
    f.write(g_BoilerPlate)
    f.write('#ifndef _BUILD_VERSION_H_\n')
    f.write('#define _BUILD_VERSION_H_\n')
    f.write('extern const char g_build_date[];\n')
    f.write('extern const char g_build_git_branch[];\n')
    f.write('extern const char g_build_git_sha[];\n')
    f.write('#endif // _BUILD_VERSION_H_\n')
    f.close()
def writeSource(location, buildDate, buildBranch, buildSha):
    f = open(os.path.join(location, 'buildVersion.c'), 'w')
    f.write(g_BoilerPlate)
    f.write('#include "buildVersion.h"\n\n')
    f.write('const char g_build_date[] = "%s";\n' % buildDate)
    f.write('const char g_build_git_branch[] = "%s";\n' % buildBranch)
    f.write('const char g_build_git_sha[] = "%s";\n' % buildSha)
    f.close()
def main():
    location = os.path.join(os.curdir, 'src')
    buildDate = time.strftime('%Y-%m-%dT%H:%M%Z', time.localtime()) # gmtime()
    buildBranch = getBranch()
    buildSha = getSha()
    writeHeader(location)
    writeSource(location, buildDate, buildBranch, buildSha)

if __name__ == '__main__':
    main()
