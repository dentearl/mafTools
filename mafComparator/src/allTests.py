#!/usr/bin/env python
"""

Copyright (C) 2009-2011 by 
Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
Mark Diekhans (markd@soe.ucsc.edu)
... and other members of the Reconstruction Team of David Haussler's 
lab (BME Dept. UCSC).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

"""
import setDiffDroppedMissing as setD
import unittest

class Intersections( unittest.TestCase ):
    def test_intersect_1( self ):
        """intersect() should return intersecting pairs. Testing starting positions"""
        missing = []
        m = setD.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 100
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = {}
        d = setD.Drop()
        d.file = 'banana.maf'
        d.name = 'banana.chr1'
        d.start = 100
        d.end = 110
        dropped[ d.name ] = [ d ]

        result = setD.intersect( dropped, missing )
        self.assertTrue( m in result )
        
    def test_intersect_2( self ):
        """intersect() should return empty list when things do not intersect. checking apple.chr1 seq position."""
        missing = []
        m = setD.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 100
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = {}
        d = setD.Drop()
        d.file = 'apple.maf'
        d.name = 'apple.chr1'
        d.start = 100
        d.end = 200
        dropped[ d.name ] = [ d ]

        result = setD.intersect( dropped, missing ) 
        self.assertTrue( len(result) == 0 )

    def test_intersect_3( self ):
        """intersect() should return empty list when the Missing sample pair is self-self."""
        missing = []
        m = setD.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 100
        m.seq2 = 'banana.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = {}
        d = setD.Drop()
        d.file = 'apple.maf'
        d.name = 'apple.chr1'
        d.start = 100
        d.end = 200
        dropped[ d.name ] = [ d ]

        result = setD.intersect( dropped, missing ) 
        self.assertTrue( len(result) == 0 )

    def test_setDiff_1( self ):
        """setDiff() should return things that do not intersect."""
        missing = []
        m = setD.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 100
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = {}
        d = setD.Drop()
        d.file = 'apple.maf'
        d.name = 'apple.chr1'
        d.start = 100
        d.end = 200
        dropped[ d.name ] = [ d ]

        iSect = setD.intersect( dropped, missing )
        result = setD.setDiff( missing, iSect ) 
        self.assertEqual( m.file, result[0].file )
        self.assertEqual( m.seq1, result[0].seq1 )
        self.assertEqual( m.pos1, result[0].pos1 )
        self.assertEqual( m.seq2, result[0].seq2 )
        self.assertEqual( m.pos2, result[0].pos2 )
    def test_intersect_4( self ):
        """intersect() should return intersecting pairs. Testing ending edges."""
        missing = []
        m = setD.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 109
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = {}
        d = setD.Drop()
        d.file = 'banana.maf'
        d.name = 'banana.chr1'
        d.start = 100
        d.end = 110
        dropped[ d.name ] = [ d ]

        result = setD.intersect( dropped, missing )
        self.assertTrue( m in result )

    def test_intersect_5( self ):
        """intersect() should return intersecting pairs. Testing ending edges."""
        missing = []
        m = setD.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 110
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = {}
        d = setD.Drop()
        d.file = 'banana.maf'
        d.name = 'banana.chr1'
        d.start = 100
        d.end = 110
        dropped[ d.name ] = [ d ]

        result = setD.intersect( dropped, missing )
        self.assertTrue( len(result) == 0 )
        
    def test_intersect_6( self ):
        """intersect() should return intersecting pairs. Testing starting positions."""
        missing = []
        m = setD.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 99
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = {}
        d = setD.Drop()
        d.file = 'banana.maf'
        d.name = 'banana.chr1'
        d.start = 100
        d.end = 110
        dropped[ d.name ] = [ d ]

        result = setD.intersect( dropped, missing )
        self.assertTrue( len(result) == 0 )

if __name__ == '__main__':
    unittest.main()
