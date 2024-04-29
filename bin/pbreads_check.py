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
import pandas as pd
import sys
import re
import time
import pysam

parser = argparse.ArgumentParser()
parser.add_argument('reference')
parser.add_argument('tsv', nargs='+')
parser.add_argument('--output-file-counts', default='concatemers_counts.csv')
parser.add_argument('--output-file-positional', default='concatemers_positional.tsv')
parser.add_argument('--output-file-zmw', default='concatemers_zmw.tsv')
args = parser.parse_args()

def reformat_deletions(str_deletions):
    if pd.isna(str_deletions):
        return(None)
    
    deletions = str_deletions.split(',')

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

    return(','.join(data_sorted))

def reformat_insertions(str_insertions):
    if pd.isna(str_insertions):
        return(None)

    insertions = str_insertions.split(',')

    data = {}

    for mutation in insertions:
        mutype, pos, mubase = mutation[0], int(mutation[1:-1]), mutation[-1]

        if pos not in data:
            data[pos] = ''
        
        data[pos] += mubase

    data_sorted = ['-%i%s' % (x, data[x]) for x in sorted(data.keys())]

    return(','.join(data_sorted))

def find_num_matches(concats, num_concats):
    num_matches = 0

    for index, row in concats.iterrows():
        if num_matches < row['num_matches']:
            num_matches = row['num_matches']
    
    return(num_matches)

def combine_mutations(concats, num_concats, mut_type):
    # store all mutations (in all concatemers)
    combined = {}

    # store mutations per concatemer
    data = []

    for index, row in concats.iterrows():
        if not pd.isna(row[mut_type]):
            concat = []

            # combine mutations
            for m in row[mut_type].split(','):
                match = re.match(r'^([ACTG-]+)(\d+)([ACTG-]+)$', m)
                mm = (match[1],int(match[2]),match[3])
                combined[mm] = 0
                concat.append(mm)
            
            # per concatemer mutations
            data.append(concat)
        else:
            # store empty array if concatemer contains no mutations
            data.append([])
        
    for m1 in combined:
        for concat in data:
            for m2 in concat:
                if m1 == m2:
                    combined[m1] += 1
    
    return(combined, data)

def check_positions(combined):
    # check: no different mutations can be observed in the SAME position.
    # - count number of mutations per position
    # - if count > 1, flag as abberant

    counts = {}

    for m in combined:
        pos = m[1]

        if pos not in counts:
            counts[pos] = 0
        
        counts[pos] += 1

        if counts[pos] > 1:
            return(False)
    
    return(True)

def check_mutations(combined, num_concats):
    # check: mutation must be present in ONE concatemer or ALL concatemers.

    for m in combined:
        if combined[m] != 1 and combined[m] != num_concats:
            return(False)
    
    return(True)

def count_errors(concats, num_concats, combined_sub, combined_del, combined_ins, reference, positional, rp, rt):
    count_errors_SDI(combined_sub, combined_del, combined_ins, num_concats, positional, rp, rt)
    count_errors_matches(concats, reference, positional, rp, rt)
    return()

def count_errors_SDI(combined_sub, combined_del, combined_ins, num_concats, positional, rp, rt):
    for m in combined_sub:
        if combined_sub[m] == num_concats:
            rp['substitutions'] += 1
            rp[m[0]+m[2]] += 1
            positional[m[1]]['rp_substitutions'] += 1
            if m[2] not in positional[m[1]]['rp_spectrum_sub']:
                positional[m[1]]['rp_spectrum_sub'][m[2]] = 0
            positional[m[1]]['rp_spectrum_sub'][m[2]] += 1

        if combined_sub[m] == 1:
            rt['substitutions'] += 1
            rt[m[0]+m[2]] += 1
            positional[m[1]]['rt_substitutions'] += 1
            if m[2] not in positional[m[1]]['rt_spectrum_sub']:
                positional[m[1]]['rt_spectrum_sub'][m[2]] = 0
            positional[m[1]]['rt_spectrum_sub'][m[2]] += 1

    for m in combined_del:
        if combined_del[m] == num_concats:
            rp['deletions'] += 1
            positional[m[1]]['rp_deletions'] += 1
            if m[0] not in positional[m[1]]['rp_spectrum_del']:
                positional[m[1]]['rp_spectrum_del'][m[0]] = 0
            positional[m[1]]['rp_spectrum_del'][m[0]] += 1

        if combined_del[m] == 1:
            rt['deletions'] += 1
            positional[m[1]]['rt_deletions'] += 1
            if m[0] not in positional[m[1]]['rt_spectrum_del']:
                positional[m[1]]['rt_spectrum_del'][m[0]] = 0
            positional[m[1]]['rt_spectrum_del'][m[0]] += 1

    for m in combined_ins:
        if combined_ins[m] == num_concats:
            rp['insertions'] += 1
            positional[m[1]]['rp_insertions'] += 1
            if m[2] not in positional[m[1]]['rp_spectrum_ins']:
                positional[m[1]]['rp_spectrum_ins'][m[2]] = 0
            positional[m[1]]['rp_spectrum_ins'][m[2]] += 1

        if combined_ins[m] == 1:
            rt['insertions'] += 1
            positional[m[1]]['rt_insertions'] += 1
            if m[2] not in positional[m[1]]['rt_spectrum_ins']:
                positional[m[1]]['rt_spectrum_ins'][m[2]] = 0
            positional[m[1]]['rt_spectrum_ins'][m[2]] += 1
    
    return()

def count_errors_matches(concats, reference, positional, rp, rt):
    ref_seq = None
    query_seq = None

    rp_matches = 0
    rt_matches = 0

    concat_no = 0

    for index, row in concats.iterrows():
        ref_seq = row['reference_seq']
        query_seq = row['query_seq']

        # count number of times each position was analyzed
        # !! assumption that reference sequence always starts with position 1
        pos = reference.find(ref_seq.replace('-','')[:20])

        if pos == -1:
            print(['ERROR reference mismatch'])
            print('reference:', reference)
            print('bam refseq:', ref_seq)
            exit()

        for b1,b2 in zip(ref_seq,query_seq):
            if b1 != '-':
                pos += 1

            if b1 == b2 and b1 != '-':
                rt_matches += 1
                rt[b1+b2] += 1
                positional[pos]['rt_matches'] += 1

                if concat_no == 0:
                    rp_matches += 1
                    rp[b1+b2] += 1
                    positional[pos]['rp_matches'] += 1

        concat_no += 1

    rp['matches'] += rp_matches
    rt['matches'] += rt_matches

    return()

def save_split_errors(movie, zmw, combined_sub, combined_del, combined_ins, num_concats, fout=sys.stdout):
    rp_sub = []
    rp_del = []
    rp_ins = []

    rt_sub = []
    rt_del = []
    rt_ins = []

    for m in combined_sub:
        if combined_sub[m] == num_concats:
            rp_sub.append('%s%i%s' % (m[0], m[1], m[2]))
        if combined_sub[m] == 1:
            rt_sub.append('%s%i%s' % (m[0], m[1], m[2]))

    for m in combined_del:
        if combined_del[m] == num_concats:
            rp_del.append('%s%i%s' % (m[0], m[1], m[2]))
        if combined_del[m] == 1:
            rt_del.append('%s%i%s' % (m[0], m[1], m[2]))

    for m in combined_ins:
        if combined_ins[m] == num_concats:
            rp_ins.append('%s%i%s' % (m[0], m[1], m[2]))
        if combined_ins[m] == 1:
            rt_ins.append('%s%i%s' % (m[0], m[1], m[2]))
    
    print(
        movie,
        zmw,

        len(rp_sub),
        len(rp_del),
        len(rp_ins),
        ','.join(rp_sub),
        ','.join(rp_del),
        ','.join(rp_ins),

        len(rt_sub),
        len(rt_del),
        len(rt_ins),
        ','.join(rt_sub),
        ','.join(rt_del),
        ','.join(rt_ins),

        sep='\t',
        file=fout)
    
    return()

def save_positional_data(output_file, positional):
    fpos = open(output_file, 'wt')

    print(
        'position',
        'base',
        'rp_matches',
        'rp_substitutions',
        'rp_deletions',
        'rp_insertions',
        'rp_spectrum_sub',
        'rp_spectrum_del',
        'rp_spectrum_ins',
        'rt_matches',
        'rt_substitutions',
        'rt_deletions',
        'rt_insertions',
        'rt_spectrum_sub',
        'rt_spectrum_del',
        'rt_spectrum_ins',
        sep='\t',
        file=fpos
        )

    for pos in sorted(positional.keys()):
        print(
            pos,
            positional[pos]['base'],
            positional[pos]['rp_matches'],
            positional[pos]['rp_substitutions'],
            positional[pos]['rp_deletions'],
            positional[pos]['rp_insertions'],
            ','.join(['%s=%i' % (k,v) for (k,v) in positional[pos]['rp_spectrum_sub'].items()]),
            ','.join(['%s=%i' % (k,v) for (k,v) in positional[pos]['rp_spectrum_del'].items()]),
            ','.join(['%s=%i' % (k,v) for (k,v) in positional[pos]['rp_spectrum_ins'].items()]),
            positional[pos]['rt_matches'],
            positional[pos]['rt_substitutions'],
            positional[pos]['rt_deletions'],
            positional[pos]['rt_insertions'],
            ','.join(['%s=%i' % (k,v) for (k,v) in positional[pos]['rt_spectrum_sub'].items()]),
            ','.join(['%s=%i' % (k,v) for (k,v) in positional[pos]['rt_spectrum_del'].items()]),
            ','.join(['%s=%i' % (k,v) for (k,v) in positional[pos]['rt_spectrum_ins'].items()]),
            sep='\t',
            file=fpos
            )

    fpos.close()

    return()

def save_counts(rp, rt, spectrum, fout=sys.stdout):
    values = []
    cols = ['matches','substitutions','deletions','insertions'] + spectrum

    for col in cols:
        values.append(rp[col])

    for col in cols:
        values.append(rt[col])

    print(','.join(cols + cols), file=fout)
    print(','.join([str(x) for x in values]), file=fout)

    return()

positional = {}
reference = ''

# initialize reference and positional counts
with pysam.FastxFile(args.reference) as fh:
    for entry in fh:
        reference = entry.sequence.upper()

        for i, b in enumerate(reference):
            positional[i+1] = {
                'base':b,

                'rp_matches':0,
                'rp_substitutions':0,
                'rp_deletions':0,
                'rp_insertions':0,
                'rp_spectrum_sub':{},
                'rp_spectrum_del':{},
                'rp_spectrum_ins':{},

                'rt_matches':0,
                'rt_substitutions':0,
                'rt_deletions':0,
                'rt_insertions':0,
                'rt_spectrum_sub':{},
                'rt_spectrum_del':{},
                'rt_spectrum_ins':{},
                }
        break

# ZMW stats
num_correct = 0
num_incorrect = 0
num_insufficient_num_concat = 0

# RP and RT error counts
rp = {'matches': 0, 'substitutions': 0, 'deletions': 0, 'insertions': 0}
rt = {'matches': 0, 'substitutions': 0, 'deletions': 0, 'insertions': 0}

# init mutational spectrum
spectrum = []

for b1 in 'ACTG':
    for b2 in 'ACTG':
        rp[b1+b2] = 0
        rt[b1+b2] = 0
        spectrum.append(b1+b2)

fout_zmw = open(args.output_file_zmw, 'wt')

print(
    'movie',
    'zmw',

    'rp_substitutions',
    'rp_deletions',
    'rp_insertions',
    'rp_spectrum_sub',
    'rp_spectrum_del',
    'rp_spectrum_ins',

    'rt_substitutions',
    'rt_deletions',
    'rt_insertions',
    'rt_spectrum_sub',
    'rt_spectrum_del',
    'rt_spectrum_ins',

    sep='\t',
    file=fout_zmw
)
# process multiple files one-by-one
for tsv in args.tsv:
    # load and reformat data
    print('Loading RT concatamers from "%s" ...' % tsv, file=sys.stderr)
    df = pd.read_csv(tsv, sep='\t')

    print('Reformatting deletions...', file=sys.stderr)
    df['deletions2'] = df.apply(lambda row: reformat_deletions(row['deletions']), axis=1)

    print('Reformatting insertions...', file=sys.stderr)
    df['insertions2'] = df.apply(lambda row: reformat_insertions(row['insertions']), axis=1)

    # extract ZMW ID from query name
    df['zmw'] = df['query_name'].str.split('/', expand=True)[1].astype(int)

    # process all ZMWs
    correct_zmw = []

    for zmw in pd.unique(df['zmw']):
        # select concatemers for a given zmw
        concats = df[df['zmw'] == zmw]
        num_concats = concats.shape[0]

        # unique PacBio movie identifier
        movie = concats['query_name'].head(1).values[0].split('/')[0]

        # determine the number of matches
        num_matches = find_num_matches(concats, num_concats)

        # we require at least 3 concatemers
        if num_concats >= 3:
            # test substitutions for consistency
            combined_sub, data_sub = combine_mutations(concats, num_concats, 'substitutions')
            ret1_sub = check_positions(combined_sub)
            ret2_sub = check_mutations(combined_sub, num_concats)

            # test deletions for consistentcy
            combined_del, data_del = combine_mutations(concats, num_concats, 'deletions2')
            ret1_del = check_positions(combined_del)
            ret2_del = check_mutations(combined_del, num_concats)

            # test insertions for consistency
            combined_ins, data_ins = combine_mutations(concats, num_concats, 'insertions2')
            ret1_ins = check_positions(combined_ins)
            ret2_ins = check_mutations(combined_ins, num_concats)

            # test if all results are successful
            correct = ret1_sub and ret2_sub and ret1_del and ret2_del and ret1_ins and ret2_ins

            if not correct:
                # log info on inconsistent ZMWs
                print('', file=sys.stderr)
                print(movie, zmw, num_concats, combined_sub, data_sub, ret1_sub, ret2_sub, sep=',', file=sys.stderr)
                print(movie, zmw, num_concats, combined_del, data_del, ret1_del, ret2_del, sep=',', file=sys.stderr)
                print(movie, zmw, num_concats, combined_ins, data_ins, ret1_ins, ret2_ins, sep=',', file=sys.stderr)

                # count inconsistent ZMWs
                num_incorrect += 1
            else:
                correct_zmw.append(zmw)

                num_correct += 1

            # count_errors_rp(combined_sub, num_concats, rp)
            count_errors(concats, num_concats, combined_sub, combined_del, combined_ins, reference, positional, rp, rt)

            # save diagnostic info (per zmw)
            save_split_errors(movie, zmw, combined_sub, combined_del, combined_ins, num_concats, fout_zmw)
        else:
            num_insufficient_num_concat += 1

    print('', file=sys.stderr)
    print('Number of correct ZMWs   :', num_correct, file=sys.stderr)
    print('Number of incorrect ZMWs :', num_incorrect, file=sys.stderr)
    print('Insufficient concatemers :', num_insufficient_num_concat, file=sys.stderr)

fout_zmw.close()

fout = open(args.output_file_counts, 'wt')
save_counts(rp, rt, spectrum, fout)
fout.close()

save_positional_data(args.output_file_positional, positional)
