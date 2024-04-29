#!/bin/bash -l

conda env create -p $PWD/conda --file conda_env.yml

conda activate $PWD/conda

bindir="$PWD/bin"

# generate ccs subreads
ccs --num-threads 16 $subreads subreads_ccs.bam

# align ccs subreads
TMPDIR=$PWD pbmm2 align $reference subreads_ccs.bam aligned_subreads_ccs.bam \
  --sort \
  --min-concordance-perc 70.0 \
  --min-length 500 \
  --sample "" \
  --num-threads 16 \
  --log-level INFO >aligned_subreads_ccs.log 2>&1

# extract RT-concatemer coordinates
$bindir/pbreads_extract.py aligned_subreads_ccs.bam > coordinates_all.txt
$bindir/pbreads_filter_subreads.py coordinates_all.txt > coordinates_filtered.txt

# split reads to RT-concatemers
$bindir/pbreads_split_subreads.py subreads_ccs.bam coordinates_filtered.txt --out concatemers.bam

# index RT-concatemers
pbindex concatemers.bam

# align RT-concatemers
TMPDIR=$PWD pbmm2 align "${reference}" concatemers.bam aligned_concatemers.bam \
  --sort \
  --min-concordance-perc 70.0 \
  --min-length 500 \
  --sample "" \
  --num-threads 16 \
  --log-level INFO >aligned_concatemers.log 2>&1

# call all RT variants
options=""

case $rname in
    DNA1_spliced)
        options="--region 118-1208"
        ;;
    IRES-Cas)
        options="--region 23-2030"
        ;;
    *)
        # catch-all category
        ;;
esac

# call RT variants
$bindir/pbreads_variant.py $options \
    --positional-output RT_positional.tsv \
    --spectrum-output RT_spectrum.csv \
    --output-file RT_variants.tsv \
    aligned_concatemers.bam $reference

# Concatmer-based calls
time $bindir/pbreads_check.py \
    --output-file-counts output.concatemers_counts.csv \
    --output-file-positional output.concatemers_positional.tsv \
    --output-file-zmw output.concatemers_zmw.tsv \
    $reference RT_variants.tsv 1>output.concatemers.log 2>&1
