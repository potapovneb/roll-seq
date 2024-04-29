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
parser.add_argument('coord')

args = parser.parse_args()

df = pd.read_csv(args.coord, sep=',')

subreads = {}

# check the mapping direction for each concatemer in PacBio subread.
# keep the number of mapping directions per subread.
for index, row in df.iterrows():
    flag = 0

    if row['flag'] & 0x10:
        flag = 1

    qname = row['qname']

    if qname not in subreads:
        subreads[qname] = set()

    subreads[qname].add(flag)

# go again through the list of concatemers.
# keep only concatemers from "good" subreads.
for index, row in df.iterrows():
    qname = row['qname']

    # we expect all concatemers map in the same direction
    if len(subreads[qname]) == 1:
        print(','.join(str(x) for x in row.values))
