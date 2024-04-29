#!/usr/bin/env python3

################################################################################
# Roll-seq measurements of RNA polymerase and reverse transcriptase fidelity
# Copyright (C) 2024 New England Biolabs, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import argparse
import pysam
import sys

parser = argparse.ArgumentParser()
parser.add_argument('bam')
args = parser.parse_args()

bamfile = pysam.AlignmentFile(args.bam, "rb")

print('qname', 'qstart', 'qend', 'qlen', 'rlen', 'strand', 'flag', sep=',')

for read in bamfile.fetch():
    if read.flag & 0x4:
        print('UNMAPPED', read.query_name, file=sys.stderr)
        continue

    if read.flag & 0x100:
        print('SECONDARY', read.query_name, file=sys.stderr)
        continue

    strand = 0
    if read.flag & 0x16:
        strand = 1

    print(
        read.query_name,
        read.query_alignment_start,
        read.query_alignment_end,
        read.query_alignment_end - read.query_alignment_start,
        read.reference_length,
        strand,
        read.flag,
        sep=',')
