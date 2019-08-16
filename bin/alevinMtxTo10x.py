#!/usr/bin/env python

# Read the results of Alevin and write outputs to a .mtx file readable by tools
# expecting 10X outputs. Adapted from
# https://github.com/k3yavi/vpolo/blob/master/vpolo/alevin/parser.py 

from __future__ import print_function
from collections import defaultdict
from struct import Struct
import pandas as pd
import gzip
import sys
import os
from scipy.io import mmread,mmwrite
from scipy.sparse import *
from shutil import copyfile
import pathlib
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Convert Alevin outputs to 10X .mtx.')
parser.add_argument('alevin_out', help = 'Alevin output directory')
parser.add_argument('mtx_out', help = 'Output directory for converted results')
parser.add_argument('--cell_prefix', dest='cell_prefix', default='', help = 'Prefix to apply to cell barcodes')
args = parser.parse_args() 

alevin_out=args.alevin_out
mtx_out=args.mtx_out
cell_prefix=args.cell_prefix

# Run some checks in the Alevin output

if not os.path.isdir(alevin_out):
    print("{} is not a directory".format( alevin_out ))
    sys.exit(1)

alevin_out = os.path.join(alevin_out, "alevin")

if not os.path.exists(alevin_out):
    print("{} directory doesn't exist".format( alevin_out ))
    sys.exit(1)

quant_file = os.path.join(alevin_out, "quants_mat.mtx.gz")
if not os.path.exists(quant_file):
    print("quant file {} doesn't exist".format( quant_file ))
    sys.exit(1)

cb_file = os.path.join(alevin_out, "quants_mat_rows.txt")
if not os.path.exists(cb_file):
    print("quant file's index: {} doesn't exist".format( cb_file ))
    sys.exit(1)

gene_file = os.path.join(alevin_out, "quants_mat_cols.txt")
if not os.path.exists(gene_file):
    print("quant file's header: {} doesn't exist".format( gene_file))
    sys.exit(1)

# Read gene and cell labels, apply cell prefix

cb_names = [cell_prefix + s for s in pd.read_csv(cb_file, header=None)[0].values]
gene_names = pd.read_csv(gene_file, header=None)[0].values
umi_counts = mmread( quant_file )
    
# Write outputs to a .mtx file readable by tools expecting 10X outputs.
# Barcodes file works as-is, genes need to be two-column, duplicating the
# identifiers. Matrix itself needs to have genes by row, so we transpose. 

pathlib.Path(mtx_out).mkdir(parents=True, exist_ok=True)
mmwrite('%s/matrix.mtx' % mtx_out, umi_counts.transpose()) 

genes_frame = pd.DataFrame([ gene_names, gene_names]).transpose()
genes_frame.to_csv(path_or_buf='%s/genes.tsv' % mtx_out, index=False, sep="\t", header = False)

with open('%s/barcodes.tsv' % mtx_out, 'w') as f:
    f.write("\n".join(cb_names))    
    f.write("\n")    
