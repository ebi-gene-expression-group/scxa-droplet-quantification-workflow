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
args = parser.parse_args() 

alevin_out=args.alevin_out
mtx_out=args.mtx_out

# Run some checks in the Alevin output

if not os.path.isdir(alevin_out):
    print("{} is not a directory".format( alevin_out ))
    sys.exit(1)

alevin_out = os.path.join(alevin_out, "alevin")
print(alevin_out)
if not os.path.exists(alevin_out):
    print("{} directory doesn't exist".format( alevin_out ))
    sys.exit(1)

quant_file = os.path.join(alevin_out, "quants_mat.gz")
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

# Read gene and cell labels

cb_names = pd.read_csv(cb_file, header=None)[0].values
gene_names = pd.read_csv(gene_file, header=None)[0].values

header_struct = Struct( "d" * len(gene_names) )

with gzip.open( quant_file ) as f:
    count = 0
    tot_read_count = 0

    # Initialise a sparse matrix with dimensions as per genes and barcodes

    umi_counts = dok_matrix((len(cb_names),len(gene_names)), dtype=np.float32)

    cell_no = 0
    gene_no = 0

    while True:

        if cell_no >0 and cell_no%100 == 0:
            print ("\r Done reading " + str(cell_no) + " cells", end= "")
            sys.stdout.flush()
        
        try:
            cellCounts = header_struct.unpack_from( f.read(header_struct.size) ) 

        except:
            print ("\nRead total " + str(cell_no) + " cells")
            print ("Found total " + str(tot_read_count) + " reads")
            break
        
        for i in range(len(cellCounts)):
            countVal = cellCounts[i]

            # Set the non-zero values in the sparse matrix

            if countVal > 0:    
                umi_counts[cell_no, i] = countVal
                tot_read_count += float(countVal)
        
        cell_no += 1

    # Write outputs to a .mtx file readable by tools expecting 10X outputs.
    # Barcodes file works as-is, genes need to be two-column, duplicating the
    # identifiers. Matrix itself needs to have genes by row, so we transpose. 

    pathlib.Path(mtx_out).mkdir(parents=True, exist_ok=True)
    mmwrite('%s/matrix.mtx' % mtx_out, umi_counts.transpose()) 
    copyfile(cb_file, '%s/barcodes.tsv' % mtx_out)
    
    genes_frame = pd.DataFrame([ gene_names, gene_names]).transpose()
    genes_frame.to_csv(path_or_buf='%s/genes.tsv' % mtx_out, index=False, sep="\t", header = False)
