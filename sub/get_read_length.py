### * Description

# Script to get read length from a fastq file (sequences must be one-liners).
#
# Usage:
# python script.py inputFile

### * Setup

### ** Import

import sys

### ** Parameters

INPUT_FILE = sys.argv[1]

### * Run

with open(INPUT_FILE, "r") as fi :
    with open(INPUT_FILE + ".lengths", "w") as fo :
        for l in fi :
            r = fi.next()
            fo.write(str(len(r.strip())) + "\n")
            fi.next()
            fi.next()
