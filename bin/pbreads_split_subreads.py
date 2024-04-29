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
import copy

parser = argparse.ArgumentParser()
parser.add_argument('bam')
parser.add_argument('coord')
parser.add_argument('--out', default='split_reads.bam')

args = parser.parse_args()

data = {}

# load concatemer coordinates for each PacBio subread
f = open(args.coord)

for line in f:
    movie, qstart, qend, qlen, rlen, strand, flag = line.strip().split(',')

    if movie not in data:
        data[movie] = []
    
    data[movie].append([int(qstart),int(qend),int(qlen),int(strand),int(flag)])

f.close()


bamfile_in = pysam.AlignmentFile(args.bam, "rb", check_sq=False)
bamfile_out = pysam.AlignmentFile(args.out, "wb", template=bamfile_in)

#          1         2         3         4         5         6
#0123456789012345678901234567890123456789012345678901234567890
#CTTTATACTCGTTTTTTTGTCTTTTTTTTTTTTCTTTGGCCACAGATCCAGGGGGCCGGTG

# go through each PacBio subread
for bamread in bamfile_in:
    if bamread.query_name not in data:
        continue

    # go through each concatemer in PacBio subread
    for concatemer in sorted(data[bamread.query_name], key=lambda x: x):
        read = copy.deepcopy(bamread)

        # concatemer coordinates and length
        qstart = concatemer[0]
        qend   = concatemer[1]
        qlen   = concatemer[2]
        strand = concatemer[3]

        if strand == 1:
            qstart1 = len(read.query_sequence) - qend
            qend1 = len(read.query_sequence) - qstart
        
            qstart = qstart1
            qend = qend1

        # m54025_211102_222528/4194815/6605_12458
        movie, zmw, region = read.query_name.split('/')

        # rename concatemer based on coordinates in PacBio subread (likely optional)
        read.query_name = '%s/%s/%i_%i' % (movie, zmw, qstart, qstart + qlen)

        # extract concatemer sequence and base quality values
        query_sequence = read.query_sequence[qstart:(qend+1)]
        query_qualities = read.query_qualities[qstart:(qend+1)]

        # update read
        read.query_sequence = query_sequence
        read.query_qualities = query_qualities

        # save concatemer to PacBio BAM file
        bamfile_out.write(read)

bamfile_in.close()
bamfile_out.close()
