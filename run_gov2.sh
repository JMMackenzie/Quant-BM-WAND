#!/bin/bash

# ---------------------- README ----------------------------

# If an 'ordered' index is desired, please rebuild atire with:
# USE_PARALLEL_INDEXING to 0 in GNUmakefile.defns, line 37
# in the external/atire/ directory, otherwise the index will
# be built using the pipeline and the order will be 
# psuedo-random

# ----------------------------------------------------------

set -f

# Atire binary
ATIRE_BIN=external/atire/bin/index

# Set collection location: This is the root of the GOV2 collection
COLLECTION_DIR=/path/to/your/gov2/collection

# Atire flag for recursive trecweb files -- This is how GOV2 is stored.
COLLECTION_TYPE="-rrtrec"

# 'find' the collection files
COLL_FILES=$(find $COLLECTION_DIR -mindepth 1 -maxdepth 1 -type d -name 'GX*' -printf '%p/*.gz ')

# Build GOV2, notify us for every 1million docs indexed, s-stem the terms,
# used quantized BM25.
$ATIRE_BIN -N1000000 -sa $COLLECTION_TYPE -iscrub:an -ts -kt -QBM25 -q ${COLL_FILES[@]}

# Build WAND index
./bin/mk_impact_idx -findex index.aspt wand-gov2 WAND

# Build BMW index
./bin/mk_impact_idx -findex index.aspt bmw-gov2 BMW

# Run WAND queries for top-1000
./bin/wand_search -q ir-repo/gov2.qry -k 1000 -z 1.0 -c wand-gov2 -o wand-k_1000-t_1.0-gov2

# Run BMW queries for top-1000
./bin/wand_search -q ir-repo/gov2.qry -k 1000 -z 1.0 -c bmw-gov2 -o bmw-k_1000-t_1.0-gov2

# Example of aggressive pruning
./bin/wand_search -q ir-repo/gov2.qry -k 1000 -z 1.2 -c wand-gov2 -o wand-k_1000-t_1.2-gov2

./bin/wand_search -q ir-repo/gov2.qry -k 1000 -z 1.2 -c bmw-gov2 -o bmw-k_1000-t_1.2-gov2
