#!/usr/bin/env python

import eval_intersectDroppedMissing as idm
import unittest

class Intersections( unittest.TestCase ):
    def test_intersect_1( self ):
        """intersect() should return intersecting pairs. Testing starting positions"""
        missing = []
        m = idm.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 100
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = []
        d = idm.Drop()
        d.file = 'banana.maf'
        d.name = 'banana.chr1'
        d.start = 100
        d.end = 110
        dropped.append( d )

        result = idm.intersect( dropped, missing ) 
        self.assertEqual( m.file, result[0].file )
        self.assertEqual( m.seq1, result[0].seq1 )
        self.assertEqual( m.pos1, result[0].pos1 )
        self.assertEqual( m.seq2, result[0].seq2 )
        self.assertEqual( m.pos2, result[0].pos2 )
        
    def test_intersect_2( self ):
        """intersect() should return empty list when things do not intersect. checking apple.chr1 seq position."""
        missing = []
        m = idm.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 100
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = []
        d = idm.Drop()
        d.file = 'apple.maf'
        d.name = 'apple.chr1'
        d.start = 100
        d.end = 200
        dropped.append( d )

        result = idm.intersect( dropped, missing ) 
        self.assertTrue( len(result) == 0 )

    def test_intersect_3( self ):
        """intersect() should return empty list when the Missing sample pair is self-self."""
        missing = []
        m = idm.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 100
        m.seq2 = 'banana.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = []
        d = idm.Drop()
        d.file = 'apple.maf'
        d.name = 'apple.chr1'
        d.start = 100
        d.end = 200
        dropped.append( d )

        result = idm.intersect( dropped, missing ) 
        self.assertTrue( len(result) == 0 )

    def test_intersectComplement_1( self ):
        """intersectComplement() should return things that do not intersect."""
        missing = []
        m = idm.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 100
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = []
        d = idm.Drop()
        d.file = 'apple.maf'
        d.name = 'apple.chr1'
        d.start = 100
        d.end = 200
        dropped.append( d )

        iSect = idm.intersect( dropped, missing )
        result = idm.intersectComplement( dropped, missing, iSect ) 
        self.assertEqual( m.file, result[0].file )
        self.assertEqual( m.seq1, result[0].seq1 )
        self.assertEqual( m.pos1, result[0].pos1 )
        self.assertEqual( m.seq2, result[0].seq2 )
        self.assertEqual( m.pos2, result[0].pos2 )
    def test_intersect_4( self ):
        """intersect() should return intersecting pairs. Testing ending edges."""
        missing = []
        m = idm.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 109
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = []
        d = idm.Drop()
        d.file = 'banana.maf'
        d.name = 'banana.chr1'
        d.start = 100
        d.end = 110
        dropped.append( d )

        result = idm.intersect( dropped, missing ) 
        self.assertEqual( m.file, result[0].file )
        self.assertEqual( m.seq1, result[0].seq1 )
        self.assertEqual( m.pos1, result[0].pos1 )
        self.assertEqual( m.seq2, result[0].seq2 )
        self.assertEqual( m.pos2, result[0].pos2 )
    def test_intersect_5( self ):
        """intersect() should return intersecting pairs. Testing ending edges."""
        missing = []
        m = idm.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 110
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = []
        d = idm.Drop()
        d.file = 'banana.maf'
        d.name = 'banana.chr1'
        d.start = 100
        d.end = 110
        dropped.append( d )

        result = idm.intersect( dropped, missing )
        self.assertTrue( len(result) == 0 )
        
    def test_intersect_6( self ):
        """intersect() should return intersecting pairs. Testing starting positions."""
        missing = []
        m = idm.Miss()
        m.file = 'banana.maf'
        m.seq1 = 'banana.chr1'
        m.pos1 = 99
        m.seq2 = 'apple.chr1'
        m.pos2 = 200
        missing.append( m )

        dropped = []
        d = idm.Drop()
        d.file = 'banana.maf'
        d.name = 'banana.chr1'
        d.start = 100
        d.end = 110
        dropped.append( d )

        result = idm.intersect( dropped, missing )
        self.assertTrue( len(result) == 0 )

if __name__ == '__main__':
    unittest.main()
