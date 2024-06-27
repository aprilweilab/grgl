#!/bin/bash
# Convert an input file (.vcf, .vcf.gz, .bgen, .igd) to an .igd file, and then compress
# the result into an .igd.lz4 file.

set -ev
echo "Converting file"
gconverter $1 $2
echo "Compressing result"
lz4 --rm -1 $2 $2.lz4
