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
mv index.aspt gov2_quantized.aspt

# IMPORTANT: If you want a frequency index instead, use the following command instead of the above one
$ATIRE_BIN -N1000000 -sa $COLLECTION_TYPE -iscrub:an -ts -kt ${COLL_FILES[@]}
mv index.aspt gov2_frequency.aspt

# Build quant WAND index
./bin/build_index -findex gov2_quantized.aspt wand-gov2-quant WAND

# Build quant BMW index
./bin/build_index -findex gov2_quantized.aspt bmw-gov2-quant BMW

# Build freq WAND index
./bin/build_index -findex gov2_frequency.aspt wand-gov2-freq WAND

# Build freq BMW index
./bin/build_index -findex gov2_frequency.aspt bmw-gov2-freq BMW

# Run WAND queries for top-1000 on the quant index for disjunctive queries
./bin/search_index -q ir-repo/gov2.qry -k 1000 -z 1.0 -c wand-gov2-quant -t OR -o test-disjunctive

# Run BMW queries for top-1000
./bin/search_index -q ir-repo/gov2.qry -k 1000 -z 1.0 -c bmw-gov2-quant -t OR -o test-disjunctive

# Example of aggressive pruning (also works for BMW)
./bin/search_index -q ir-repo/gov2.qry -k 1000 -z 1.2 -c wand-gov2-quant -t OR -o test-aggro-disjunctive

# Example of ranked conjunctive processing on the frequency index (also works for BMW)
./bin/search_index -q ir-repo/gov2.qry -k 1000 -z 1.0 -c wand-gov2-freq -t AND -o test-conjunctive

