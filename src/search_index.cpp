#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <ctime>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "query.hpp"
#include "invidx.hpp"
#include "generic_rank.hpp"
#include "impact.hpp"
#include "bm25.hpp"
#include "util.hpp"

typedef struct cmdargs {
    std::string collection_dir;
    std::string query_file;
    std::string postings_file;
    std::string doclen_file;
    std::string global_file;
    std::string output_prefix;
    std::string index_type_file;
    uint64_t k;
    double F_boost;
    query_traversal traversal;
    std::string traversal_string;
} cmdargs_t;

void print_usage(std::string program) {
  std::cerr << program << " -c <collection>"
                       << " -q <query_file>"
                       << " -k <no. items to retrieve>"
                       << " -z <F: aggression parameter. 1.0 is rank-safe>"
                       << " -o <output file handle>"
                       << " -t <traversal type: AND|OR>"
                       << std::endl;
  exit(EXIT_FAILURE);
}

cmdargs_t
parse_args(int argc, char* const argv[])
{
  cmdargs_t args;
  int op;
  args.collection_dir = "";
  args.output_prefix = "";
  args.traversal = UNKNOWN;
  args.traversal_string = "";
  args.k = 10;
  args.F_boost = 1.0;
  while ((op=getopt(argc,argv,"c:q:k:z:o:t:")) != -1) {
    switch (op) {
      case 'c':
        args.collection_dir = optarg;
        args.postings_file = args.collection_dir + "/WANDbl_postings.idx";
        args.doclen_file = args.collection_dir +"/doc_lens.txt";
        args.global_file = args.collection_dir +"/global.txt";
        args.index_type_file = args.collection_dir + "/index_info.txt";
        break;
      case 'o':
        args.output_prefix = optarg;
        break;
      case 'q':
        args.query_file = optarg;
        break;
      case 'k':
        args.k = std::strtoul(optarg,NULL,10);
        break;
      case 'z':
        args.F_boost = atof(optarg);
        break;
      case 't':
        args.traversal_string = optarg;
        if (args.traversal_string == "OR")
          args.traversal = OR;
        else if (args.traversal_string == "AND")
          args.traversal = AND;
        else 
          print_usage(argv[0]);
        break;
      case '?':
      default:
        print_usage(argv[0]);
    }
  }
  if (args.collection_dir=="" || args.query_file=="" || args.F_boost < 1 ||
      args.traversal == UNKNOWN) {
    std::cerr << "Missing/Incorrect command line parameters.\n";
    print_usage(argv[0]);
  }
  return args;
}

int 
main (int argc,char* const argv[])
{
  /* define types */
  using plist_type = block_postings_list<128>;
  using my_index_t = idx_invfile<plist_type, generic_rank>;
  using clock = std::chrono::high_resolution_clock;

  /* parse command line */
  cmdargs_t args = parse_args(argc,argv);

  std::cerr << "NOTE: Global F boost = " << args.F_boost << std::endl;

  // Read the index and traversal type
  std::ifstream read_type(args.index_type_file);
  std::string t_traversal, t_postings;
  read_type >> t_traversal;
  read_type >> t_postings;
  
  // Wand or BMW index? 
  index_form t_index_type;
  if (t_traversal == STRING_WAND)
    t_index_type = WAND;
  else if (t_traversal == STRING_BMW)
    t_index_type = BMW;
  else {
    std::cerr << "Index is corrupted. Please rebuild." << std::endl;
    exit(EXIT_FAILURE);
  }

  // TF or a quant index?
  postings_form t_postings_type;
  if (t_postings == STRING_FREQ) {
    t_postings_type = FREQUENCY;
  }
  else if (t_postings == STRING_QUANT) {
    t_postings_type = QUANTIZED;
  }
  else {
    std::cerr << "Index is corrupted. Please rebuild." << std::endl;
    exit(EXIT_FAILURE);
  }

 
  /* parse queries */
  std::cout << "Parsing query file '" << args.query_file << "'" << std::endl;
  auto queries = query_parser::parse_queries(args.collection_dir,args.query_file);
  std::cout << "Found " << queries.size() << " queries." << std::endl;

  std::string index_name(basename(strdup(args.collection_dir.c_str())));

  /* load the index */
  my_index_t index;
 
  auto load_start = clock::now();
  // Construct index instance.
  construct(index, args.postings_file, args.F_boost);

  // Prepare Ranker
  uint64_t temp;
  std::vector<uint64_t>doc_lens;
  ifstream doclen_file(args.doclen_file);
  if(!doclen_file.is_open()){
    std::cerr << "Couldn't open: " << args.doclen_file << std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout << "Reading document lengths." << std::endl;
  /*Read the lengths of each document from asc file into vector*/
  while(doclen_file >> temp){
    doc_lens.push_back(temp);
  }
  ifstream global_file(args.global_file);
  if(!global_file.is_open()) {
    std::cerr << "Couldn't open: " << args.global_file << std::endl;
    exit(EXIT_FAILURE);
  } 
  // Load the ranker
  uint64_t total_docs, total_terms;
  global_file >> total_docs >> total_terms;
  index.load(doc_lens, total_terms, total_docs, t_postings_type);

  auto load_stop = clock::now();
  auto load_time_sec = std::chrono::duration_cast<std::chrono::seconds>(load_stop-load_start);
  std::cout << "Index loaded in " << load_time_sec.count() << " seconds." << std::endl;

  /* process the queries */
  std::map<uint64_t,std::chrono::microseconds> query_times;
  std::map<uint64_t,result> query_results;
  std::map<uint64_t,uint64_t> query_lengths;

  size_t num_runs = 3;
  std::cerr << "Times are the average across " << num_runs << " runs." << std::endl;
  for(size_t i = 0; i < num_runs; i++) {
    // For each query
    for(const auto& query: queries) {
      auto id = std::get<0>(query);
      auto qry_tokens = std::get<1>(query);
      std::cout << "[" << id << "] |Q|=" << qry_tokens.size(); 
      std::cout.flush();

      // run the query
      auto qry_start = clock::now();
      auto results = index.search(qry_tokens,args.k, t_index_type, args.traversal);
      auto qry_stop = clock::now();

      auto query_time = std::chrono::duration_cast<std::chrono::microseconds>(qry_stop-qry_start);
      std::cout << " TIME = " << std::setprecision(5)
                << query_time.count() / 1000.0 
                << " ms" << std::endl;

      auto itr = query_times.find(id);
      if(itr != query_times.end()) {
        itr->second += query_time;
      } else {
        query_times[id] = query_time;
      }

      if(i==0) {
        query_results[id] = results;
        query_lengths[id] = qry_tokens.size();
      }
    }
  }


  // generate output string
  args.output_prefix = args.output_prefix + "-" // user specified
                       + t_postings + "-"  // quantized or frequency
                       + t_traversal + "-" // wand or bmw 
                       + args.traversal_string + "-" // OR, AND, etc
                       + std::to_string(args.k) + "-" // no. results
                       + std::to_string(args.F_boost);

  // Average the times
  for(auto& timing : query_times) {
    timing.second = timing.second / num_runs;
  }

  std::string time_file = args.output_prefix + "-time.log";

  /* output */
  std::cout << "Writing timing results to '" << time_file << "'" << std::endl;     
  std::ofstream resfs(time_file);
  if(resfs.is_open()) {
    resfs << "query;num_results;postings_eval;docs_fully_eval;docs_added_to_heap;threshold;num_terms;time_ms;traversal_type" << std::endl;
    for(const auto& timing: query_times) {
      auto qry_id = timing.first;
      auto qry_time = timing.second;
      auto results = query_results[qry_id];
      resfs << qry_id << ";" << results.list.size() << ";" 
            << results.postings_evaluated << ";"
            << results.docs_fully_evaluated << ";" 
            << results.docs_added_to_heap << ";" 
            << results.final_threshold << ";" 
            << query_lengths[qry_id] << ";" 
            << qry_time.count() / 1000.0 << ";"
            << args.traversal_string << std::endl;
    }
  } else {
    perror ("Could not output results to file.");
  }

  // Write TREC output file.

  /* load the docnames map */
  std::unordered_map<uint64_t,std::string> id_mapping;
  std::string doc_names_file = args.collection_dir + "/doc_names.txt";
  std::ifstream dfs(doc_names_file);
  size_t j=0;
  std::string name_mapping;
  while( std::getline(dfs,name_mapping) ) {
    id_mapping[j] = name_mapping;
    j++;
  }
 

  std::string trec_file = args.output_prefix + "-trec.run";
  std::cout << "Writing trec output to " << trec_file << std::endl;
  std::ofstream trec_out(trec_file);
  if(trec_out.is_open()) {
    for(const auto& result: query_results) {
      auto qry_id = result.first;
      auto qry_res = result.second.list;
      for(size_t i=1;i<=qry_res.size();i++) {
        trec_out << qry_id << "\t"
                 << "Q0" << "\t"
                 << id_mapping[qry_res[i-1].doc_id] << "\t"
                 << i << "\t"
                 << qry_res[i-1].score << "\t"  
                 << "WANDbl" << std::endl;
      }
    }
  } else {
    perror ("Could not output results to file.");
  }

  return EXIT_SUCCESS;
}
