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
#include <getopt.h>
#include <math.h>  // ceil()
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "sharedMaf.h"
#include "mafCoverage.h"
#include "mafCoverageAPI.h"
#include "buildVersion.h"

const char *g_version = "version 0.1 May 2013";
uint64_t getRegionSize(char *seq1, stHash *intervalsHash);

void version(void) {
  fprintf(stderr, "mafCoverage, %s\nbuild: %s, %s, %s\n\n", g_version,
          g_build_date, g_build_git_branch, g_build_git_sha);
}

void usage(void) {
  version();
  fprintf(stderr, "Usage: mafCoverage --maf [maf file] \n\n");
  fprintf(stderr, "Options: \n");
  usageMessage('h', "help", "show this help message and exit.");
  usageMessage('m', "maf", "path to maf file.");
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
  return EXIT_SUCCESS;
}
