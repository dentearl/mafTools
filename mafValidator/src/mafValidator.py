#!/usr/bin/env python2.6
"""
mafValidator
10 Oct 2011
dent earl, dearl (a) soe ucsc edu

Script to validate Multpile Alignment Format (maf) 
files.

"""
##############################
# Copyright (C) 2009-2012 by
# Dent Earl (dearl@soe.ucsc.edu, dent.earl@gmail.com)
#
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
##############################
from optparse import OptionParser
import os
import re
import sys

class ValidatorError(Exception): pass
class SourceLengthError(ValidatorError): pass
class SpeciesFieldError(ValidatorError): pass
class AlignmentLengthError(ValidatorError): pass
class FieldNumberError(ValidatorError): pass
class FooterError(ValidatorError): pass
class StrandCharacterError(ValidatorError): pass
class StartFieldError(ValidatorError): pass
class SourceSizeFieldError(ValidatorError): pass
class OutOfRangeError(ValidatorError): pass
class HeaderError(ValidatorError): pass

def initOptions(parser):
   parser.add_option('--maf', dest = 'filename', 
                     help = 'path to maf file to validate.')
   parser.add_option('--testChromNames', dest = 'testChromNames', 
                     action = 'store_true',
                     default = False,
                     help = ('Test that all species fields contain chrom name, i.e.: '
                             's hg19.chr1 ... default = %default'))
   
def checkOptions(options, args, parser):
   if options.filename is None:
      parser.error('specify --maf')
   if not os.path.exists(options.filename):
      parser.error('--maf %s does not exist.' % options.filename)

def validateMaf(filename, testChromNames = False):
   """ returns true on valid maf file
   """ 
   nameRegex = r'(.+?)\.(chr.+)'  
   namePat = re.compile(nameRegex)
   f = open(filename, 'r')
   header = f.next()
   validateHeader(header, filename)
   sources = {}
   for lineno, line in enumerate(f, 2):
      line = line.strip()
      if line.startswith('#'):
         continue
      if line.startswith('s'):
         species, chrom, length = validateSeqLine(namePat, testChromNames, lineno, line, filename)
         if (species, chrom)  in sources:
            if sources[(species, chrom)] != length:
               raise SourceLengthError('maf %s has different source lengths for '
                                       'species %s lineno %d: %s'
                                       % (filename, species, lineno, line))
         else:
            sources[(species, chrom)] = length
         
   if line != '':
      raise FooterError('maf %s has a bad footer, should end with two new lines.' 
                        % filename)
   return True

def validateSeqLine(namePat, testChromNames, lineno, line, filename):
   data = line.split()
   if len(data) != 7:
      raise FieldNumberError('maf %s has incorrect number of fields on lineno %d: %s' 
                             % (filename, lineno, line))
   if data[4] not in ('-', '+'):
      raise StrandCharacterError('maf %s has unexpected character in strand field "%s" lineno %d: %s' 
                                 % (filename, data[4], lineno, line))
   if int(data[3]) != len(data[6].replace('-', '')):
      raise AlignmentLengthError('maf %s has incorrect seq len or alignment field lineno %d: %s'
                                 % (filename, lineno, line))
   if int(data[2]) < 0:
      raise StartFieldError('maf %s has bad start field lineno %d: %s'
                           % (filename, lineno, line))
   if int(data[5]) < 0:
      raise SourceSizeFieldError('maf %s has bad srcSize field lineno %d: %s'
                                 % (filename, lineno, line))
   if int(data[2]) + int(data[3]) > int(data[5]):
      raise OutOfRangeError('maf %s out of range sequence lineno %d: %s'
                            % (filename, lineno, line))
   if testChromNames:
      m = re.match(namePat, data[1])
      if m is None:
         raise SpeciesFieldError('maf %s has name (src) field without ".chr" suffix: "%s" lineno %d: %s' 
                                 % (filename, data[1], lineno, line))
      return m.group(1), m.group(2), data[5]
   return data[1], None, data[5]

def validateHeader(header, filename):
   """ tests the first line of the maf file to make sure it is valid
   """
   if not header.startswith('##'):
      raise HeaderError('maf %s has bad header, fails to start with #: %s' 
                           % (filename, header))
   data = header.split()
   version = False
   for d in data[1:]:
      if d.startswith('=') or d == '=' or d.endswith('='):
         raise HeaderError('maf %s has bad header, there may be '
                              'no whitespace surrounding "=": %s' 
                              % (filename, header))
      try:
         k,v = d.split('=')
      except ValueError:
         raise HeaderError('maf %s has bad header, there may be '
                              'no whitespace surrounding "=": %s' 
                              % (filename, header))
      if k == 'version':
         version = True
   if not version:
      raise HeaderError('maf %s has bad header, no version information: %s' 
                           % (filename, header))
         
def main():
   parser = OptionParser()
   initOptions(parser)
   options, args = parser.parse_args()
   checkOptions(options, args, parser)
   
   validateMaf(options.filename, options.testChromNames)

if __name__ == '__main__':
   main()
