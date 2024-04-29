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
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument('bam')
parser.add_argument('reference')
parser.add_argument('--mapq', type=int, default=60, help='minimum mapping quality')
parser.add_argument('--np', type=int, default=15, help='minimum number of passes')
parser.add_argument('--rq', type=float, default=1.0, help='minimum read quality')
parser.add_argument('--qual', type=int, default=93, help='minimum variant call quality')
parser.add_argument('--positional-output', type=str, default='positional.tsv')
parser.add_argument('--spectrum-output', type=str, default='spectrum.tsv')
parser.add_argument('--output-file', type=str, default='variants.tsv')
parser.add_argument('--exclude-file', type=str, default='')
parser.add_argument('--region', type=str, nargs='+')
parser.add_argument('--zmw', type=int, default=-1)
parser.add_argument('--verbose', action='store_true', default=False)

args = parser.parse_args()

def reformat_insertions(insertions):
    data = {}

    for mutation in insertions:
        mutype, pos, mubase = mutation[0], int(mutation[1:-1]), mutation[-1]

        if pos not in data:
            data[pos] = ''
        
        data[pos] += mubase

    data_sorted = ['-%i%s' % (x, data[x]) for x in sorted(data.keys())]

    return(data_sorted)

def reformat_deletions(deletions):
    data = {}
    count = 0

    for mutation in deletions:
        mubase, pos, mutype = mutation[0], int(mutation[1:-1]), mutation[-1]

        if pos - count not in data:
            count = 0
            data[pos - count] = ''

        data[pos - count] += mubase
        count += 1

    data_sorted = ['%s%i-' % (data[x], x) for x in sorted(data.keys())]

    return(data_sorted)

reference_sequence = None

with pysam.FastxFile(args.reference) as fh:
    for entry in fh:
        reference_sequence = entry.sequence.upper()
        break

df_exclude = pd.DataFrame()

if args.exclude_file != '':
    df_exclude = pd.read_csv(args.exclude_file, header=None)

allowed_positions = []

for region in args.region:
    start, end = region.split('-')
    start = int(start)
    end = int(end)

    # reference positions are zero-based, while region defintions are 1-based
    for pos in range(start-1,end):
        allowed_positions.append(pos)

fvariants = open(args.output_file, 'wt')

print(
    'query_name',
    'query_length',

    'np',
    'rq',

    'reference_start',
    'reference_end',
    'query_start',
    'query_end',
    'cigar',

    'num_matches',
    'num_substitutions',
    'num_deletions',
    'num_insertions',

    'substitutions',
    'deletions',
    'insertions',

    'reference_seq',
    'query_seq',
    'qual_str',
    'alignment_str',

    sep='\t',
    file=fvariants
)

### predefined substitution types
pairs = []

### mutational spectrum (substitutions)
spectrum = {}

for rb in 'ACGT':
    for qb in 'ACGT':
        spectrum[rb + qb] = 0
        pairs.append(rb + qb)

### positional frequencies for matches, subs, dels, ins
positional = {}

for pos, base in enumerate(reference_sequence):
    ### use 1-based numbering to tabular output
    positional[pos+1] = {
        'base'          : base,
        'matches'       : 0,
        'substitutions' : 0,
        'deletions'     : 0,
        'insertions'    : 0,
    }

bamfile = pysam.AlignmentFile(args.bam, "rb", check_sq=False)

while True:
    try:
        read = next(bamfile)
    except StopIteration:
        break
    
    movie, zmw, *other = read.query_name.split('/')
    zmw = int(zmw)

    if (not df_exclude.empty) and (zmw in df_exclude.iloc[:, 0].values):
        continue

    if read.is_unmapped:
        if args.verbose: print('WARN: unmapped')
        continue

    if read.is_secondary:
        if args.verbose: print('WARN: secondary')
        continue
    
    if read.is_supplementary:
        if args.verbose: print('WARN: supplementary')
        continue

    if read.mapping_quality < args.mapq:
        if args.verbose: print('WARN: mapping_quality')
        continue
    
    if read.get_tag('np') < args.np:
        if args.verbose: print('WARN: np')
        continue
    
    if read.get_tag('rq') < args.rq:
        if args.verbose: print('WARN: rq')
        continue

    # make sure that the read covers the entire reference sequence
    if read.reference_start > 5 or read.reference_end < len(reference_sequence) - 5:
        if args.verbose: print('WARN: reference span')
        continue

    # make sure that almost the entire read maps to reference
    if read.query_alignment_start > 5 or read.query_alignment_end < len(read.query_sequence) - 5:
        if args.verbose: print('WARN: alignment span')
        continue
    
    if args.zmw != -1 and args.zmw != zmw:
        if args.verbose: print('WARN: zmw')
        continue

    aligned_pairs = read.get_aligned_pairs()

    qseq = read.query_sequence
    rseq = reference_sequence
    qual = read.query_qualities

    query_aligned     = ''
    reference_aligned = ''
    quality_aligned   = ''
    alignment         = ''

    matches       = 0
    substitutions = []
    deletions     = []
    insertions    = []

    ### BAM file has 0-based numbering
    last_rp = 0
    last_qp = 0

    for (qp,rp) in aligned_pairs:
        # print(rp, qp, last_rp)

        if qp != None and rp != None:
            # ---------- Substitutions ---------- #
            # if rp >= (args.start-1) and rp <= (args.end-1):
            if rp in allowed_positions:
                ### alignment strings
                query_aligned += qseq[qp]
                reference_aligned += rseq[rp]
                quality_aligned += chr(qual[qp] + 33)

                pair = rseq[rp] + qseq[qp]

                if qseq[qp] != rseq[rp]:
                    alignment += 'X'

                    ### we use 1-based numbering for the tabular output
                    substitutions.append('%s%i%s' % (rseq[rp], (rp+1), qseq[qp]))

                    spectrum[pair] += 1
                    positional[(rp+1)]['substitutions'] += 1
                else:
                    alignment += '.'
                    matches += 1
                    spectrum[pair] += 1
                    positional[(rp+1)]['matches'] += 1

            last_rp = rp
            last_qp = qp

            pass
        elif qp != None and rp == None:
            # ---------- Insertions ---------- #
            # if last_rp >= (args.start-1) and last_rp <= (args.end-1):
            if last_rp in allowed_positions:
                ### alignment strings
                query_aligned += qseq[qp]
                reference_aligned += '-'
                quality_aligned += chr(qual[qp] + 33)
                alignment += '.'

                ### we use 1-based numbering for the tabular output
                insertions.append('%s%i%s' % ('-', (last_rp+1), qseq[qp]))

                positional[(last_rp+1)]['insertions'] += 1

            last_qp = qp

            pass
        elif qp == None and rp != None:
            # ---------- Deletions ---------- #
            # if last_rp >= (args.start-1) and last_rp <= (args.end-1):
            if rp in allowed_positions:
                ### alignment strings
                query_aligned += '-'
                reference_aligned += rseq[rp]
                quality_aligned += ' '
                alignment += '.'

                ### we use 1-based numbering for the tabular output
                deletions.append('%s%i%s' % (rseq[rp], (rp+1), '-'))

                positional[(rp+1)]['deletions'] += 1

            last_rp = rp

            pass
        elif qp == None and rp == None:

            pass

    # insertions = reformat_insertions(insertions)
    # deletions = reformat_deletions(deletions)

    print(
        read.query_name,
        len(read.query_sequence),

        read.get_tag('np'),
        read.get_tag('rq'),

        read.reference_start,
        read.reference_end,
        read.query_alignment_start,
        read.query_alignment_end,
        read.cigarstring,

        matches,
        len(substitutions),
        len(deletions),
        len(insertions),

        ','.join(substitutions),
        ','.join(deletions),
        ','.join(insertions),

        reference_aligned,
        query_aligned,
        quality_aligned,
        alignment,

        sep='\t',
        file=fvariants
        )

bamfile.close()

fvariants.close()

### write to file spectrum information
f = open(args.spectrum_output, 'wt')

print(','.join(pairs), file=f)
print(','.join([ str(spectrum[pair]) for pair in pairs]), file=f)

f.close()

### write to file positional information
f = open(args.positional_output, 'wt')

for pos in sorted(positional.keys()):
    print(
        pos,
        positional[pos]['base'],
        positional[pos]['matches'],
        positional[pos]['substitutions'],
        positional[pos]['deletions'],
        positional[pos]['insertions'],
        sep='\t', file=f)

f.close()
