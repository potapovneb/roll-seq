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

parser = argparse.ArgumentParser()
parser.add_argument('summary')
parser.add_argument('type', type=str, default='RP')
parser.add_argument('samples', type=str, default='13,14,15,16')
args = parser.parse_args()

df = pd.read_csv(args.summary)

samples = [int(x) for x in args.samples.split(',')]

df = df[df['sample'].isin(samples)]
df['total']  = df['substitutions'] + df['deletions'] + df['insertions']

df['rp_sub'] = df['substitutions'] / df['matches']
df['rp_del'] = df['deletions'] / df['matches']
df['rp_ins'] = df['insertions'] / df['matches']
df['rp_err'] = df['total'] / df['matches']

df['total.1'] = df['substitutions.1'] + df['deletions.1'] + df['insertions.1']

df['rt_sub'] = df['substitutions.1'] / df['matches.1']
df['rt_del'] = df['deletions.1'] / df['matches.1']
df['rt_ins'] = df['insertions.1'] / df['matches.1']
df['rt_err'] = df['total.1'] / df['matches.1']

total = df.sum()
std = df.std()

if args.type == 'RP':
    sub_total = round(total['substitutions'] / total['matches'] * 1e6)
    del_total = round(total['deletions'] / total['matches'] * 1e6)
    ins_total = round(total['insertions'] / total['matches'] * 1e6)
    err_total = round((total['total']) / total['matches'] * 1e6)

    sub_std = round(df['rp_sub'].std() * 1e6)
    del_std = round(df['rp_del'].std() * 1e6)
    ins_std = round(df['rp_ins'].std() * 1e6)
    err_std = round(df['rp_err'].std() * 1e6)

elif args.type == 'RT':
    sub_total = round(total['substitutions.1'] / total['matches.1'] * 1e6)
    del_total = round(total['deletions.1'] / total['matches.1'] * 1e6)
    ins_total = round(total['insertions.1'] / total['matches.1'] * 1e6)
    err_total = round(total['total.1'] / total['matches.1'] * 1e6)

    sub_std = round(df['rt_sub'].std() * 1e6)
    del_std = round(df['rt_del'].std() * 1e6)
    ins_std = round(df['rt_ins'].std() * 1e6)
    err_std = round(df['rt_err'].std() * 1e6)

print(
    args.samples,
    '%i ± %i' % (sub_total, sub_std),
    '%i ± %i' % (del_total, del_std),
    '%i ± %i' % (ins_total, ins_std),
    '%i ± %i' % (err_total, err_std),
    sep='\t')
