/*
 * Copyright (C) 2011-2013 by
 * Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
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
#include <assert.h>
#include <errno.h> // file existence via ENOENT, errno
#include <getopt.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "sharedMaf.h"
#include "mafCoverageAPI.h"
#include "bioioC.h" // benLine()


bool is_wild(const char *s) {
  // return true if char array ends in *, false otherwise
  if (s[strlen(s) - 1] == '*')
    return true;
  return false;
}
bool searchMatched(mafLine_t *ml, const char *seq) {
  // report false if search did not match, true if it did
  if (maf_mafLine_getType(ml) != 's')
    return false;
  return searchMatched_(maf_mafLine_getSpecies(ml), seq);
}
bool searchMatched_(const char *target, const char *seq) {
  if (is_wild(seq)) {
    // only compare up to the wildcard character
    if (!(strncmp(target, seq, strlen(seq) - 1) == 0)) {
      return false;
    }
  } else {
    if (!(strcmp(target, seq) == 0)) {
      return false;
    }
  }
  return true;
}
