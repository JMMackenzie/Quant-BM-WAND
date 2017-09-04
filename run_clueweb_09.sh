#!/bin/bash

# ---------------------- README ----------------------------

# If an 'ordered' index is desired, please rebuild atire with:
# USE_PARALLEL_INDEXING to 0 in GNUmakefile.defns, line 37
# in the external/atire/ directory, otherwise the index will
# be built using the pipeline and the order will be 
# psuedo-random

# Please note that this binary will also suffice for Clueweb12.
# Just change where you point, and it should just_work_(tm).

# ----------------------------------------------------------

set -f

# Atire binary
ATIRE_BIN=external/atire/bin/index

# Set collection location and format for atire. 
# For Clueweb09B, this is the root of the English 1 portion.
COLLECTION_DIR=/research/remote/collections/TREC/ClueWeb09/disk1/ClueWeb09_English_1/

# Atire flag for recursive warcgz files -- This is how the Clueweb documents are stored
COLLECTION_TYPE="-rrwarcgz"

# 'find' the required files
COLL_FILES=$(find $COLLECTION_DIR -mindepth 1 -maxdepth 2 -type d -name 'en*' -printf '%p/*warc.gz ')

# Build CW*, notify us for every 1million docs indexed, s-stem the terms,
# used quantized BM25. Note that -iscrub:un is used instead of :an (as used in Gov2).
$ATIRE_BIN -N1000000 -sa $COLLECTION_TYPE -iscrub:un -ts -kt -QBM25 -q ${COLL_FILES[@]}

# Important: To build the frequency version, use the following line instead of the above line
# $ATIRE_BIN -N1000000 -sa $COLLECTION_TYPE -iscrub:un -ts -kt ${COLL_FILES[@]}

# Build WAND index
./bin/build_index -findex index.aspt wand-cw09b-quant WAND

# Build BMW index
./bin/build_index -findex index.aspt bmw-cw09b-quant BMW

# Run WAND queries for top-1000
./bin/search_index -q ir-repo/cw09b.qry -k 1000 -z 1.0 -c wand-cw09b-quant -t OR -o wand

# Run BMW queries for top-1000
./bin/search_index -q ir-repo/cw09b.qry -k 1000 -z 1.0 -c bmw-cw09b-quant -t OR -o bmw

