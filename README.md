ATIRE-derived Quantized WAND/BMW
======================

Reference
---------
If you use this code for your own experiments, please cite the following work:


Matt Crane, J. Shane Culpepper, Jimmy Lin, Joel Mackenzie and Andrew Trotman,
**A Comparison of Document-at-a-Time and Score-at-a-Time Query Evaluation.**
In Proceedings of the Tenth ACM International Conference on Web Search and Data Mining (WSDM 2017)

Acknowledgements
----------------
This codebase was derived from many seperate codebases, including:

- [ATIRE](https://github.com/snapbug/atire)
- [WandBL](https://github.com/jsc/WANDbl)
- [JASS](https://github.com/lintool/JASS)
- [FastPFor](https://github.com/lemire/FastPFor)
- [FastDifferentialCoding](https://github.com/lemire/FastDifferentialCoding)
- [SDSL](https://github.com/simongog/sdsl-lite)

Please see each of these repo's for further information.

Licence
-------
The licence is provided, see the LICENCE file.

Usage
=====
Since all of the systems were implemented to use the same underlying index,
which is an index from the ATIRE search engine, this ATIRE index must be 
built before the Wand/BMW indexes can be built.

We have provided an end-to-end script which will build the appropriate ATIRE
index, construct the Wand/BMW indexes from this, and then run the subsequent
queries. Please see: `run_gov2.sh` and `run_clueweb_09.sh`. 
Note that the ATIRE syntax differs between GOV2 and ClueWeb collections,
and that the `run_clueweb_09.sh` script will work for ClueWeb12 too.

WAND/BMW
========
Index building
--------------
Building an index is simple:

`./bin/build_index -findex <atire_index> <output_index_directory> <type>`

Type can either be `BMW` or `WAND`.

Important: Now that frequency indexes are supported, the type of index that will be output by
the `build_index` binary will depend on the type of ATIRE index in which you pass in. For example,
if you build a frequency index in ATIRE, `build_index` will automatically deduce this and thus
build a tf index for `WAND` or `BMW`. Take a look at `run_gov2.sh` for specific examples on
how to build a frequency/quantized ATIRE index.

Running queries
---------------
Running queries should also be simple. We have provided some query files in
the `ir-repo/` directory.

`./bin/search_index -c <collection> -k <no. results> -q <query_file> -o <output_prefix> -t <traversal strategy: AND|OR>`
- collection corresponds to the index built above (output_index_directory)
- `<output prefix>`: Two files will be output: `*-time.log` and `*-trec.run`.
The `*-trec.run` file is directly usable with `trec_eval`.
- `-z` specifies the aggression parameter: A float between 1.0 and infinity.
- `-t` specifies whether you want conjunctive or disjunctive processing. If you have a block-max index and use -t AND, this will run block-max AND (and so on).

JASS
====
The instructions and code for the JASS engine can be found on [Github](https://github.com/lintool/JASS)

Updates
=======
* 24/5/2017: Added support for conjunctive Wand/Block-Max AND querying. Please
note that BM-AND is currently untested.
* 26/5/2017: Fix added for Block-Max AND. BMA tested and working as expected.
* 04/9/2017: Support added for frequency indexes, so pre-quantization is not required (that is, we calculate our ranking function at query processing time).

Wish List/To Do
===============
* In the original [WandBL](https://github.com/jsc/WANDbl) repo, the index is built using Indri as the initial indexer. This repo has the tools and code
necessary to use Indri to build consumable indexes for this particular repo. The aim is to integrate this code here on a rainy day, but would also require
the quantizer to be built in here too. 

