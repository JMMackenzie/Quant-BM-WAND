#include <iostream>
#include "ant_param_block.h"
#include "search_engine.h"
#include "search_engine_btree_leaf.h"
#include "btree_iterator.h"
#include "memory.h"
#include "include/impact.hpp"
#include "sdsl/int_vector_buffer.hpp"
#include "include/block_postings_list.hpp"
#include "include/util.hpp"

const static size_t INIT_SZ = 4096; 
const static size_t INDRI_OFFSET = 2; // Indri offsets terms 0 and 1 as special


int main(int argc, char **argv)
{
	ANT_ANT_param_block params(argc, argv);
	long last_param = params.parse();

	if (last_param == argc)
	{
		std::cout << "USAGE: " << argv[0];
		std::cout << " [ATIRE options] <collection folder> <index_type>\n" 
              << " index type can be `BMW` or `WAND`" << std::endl;
		return EXIT_FAILURE;
	}
	using clock = std::chrono::high_resolution_clock;

	std::string collection_folder = argv[argc-2];
  std::string s_index_type = argv[argc-1];
	create_directory(collection_folder);
	std::string dict_file = collection_folder + "/dict.txt";
	std::string doc_names_file = collection_folder + "/doc_names.txt";
	std::string postings_file = collection_folder + "/WANDbl_postings.idx";
	std::string global_info_file = collection_folder + "/global.txt";
	std::string doclen_tfile = collection_folder + "/doc_lens.txt";
  std::string index_type_file = collection_folder + "/index_info.txt";

	std::ofstream doclen_out(doclen_tfile);

	auto build_start = clock::now();

  // Select the index format - BMW or WAND?
  index_form index_format;
  if (s_index_type == STRING_BMW) {
    index_format = BMW;
  }
  else if (s_index_type == STRING_WAND) {
    index_format = WAND;
  }
  else {
    std::cerr << "Incorrect index type specified. Exiting." << std::endl;
    return EXIT_FAILURE;
  }

  // For reference (later, for a user), write out which index type this is
  std::ofstream index_file_output(index_type_file);
  index_file_output << s_index_type;

  // load stuff
  ANT_memory memory;
  ANT_search_engine search_engine(&memory);
  search_engine.open(params.index_filename);

  if (!search_engine.quantized()) {
    std::cerr << "Must use a quantized index for conversion" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Keep track of term ordering
  unordered_map<string, uint64_t> map;

  std::cout << "Writing global info to " << global_info_file << "."
            << std::endl;
  std::vector<std::string> document_names;

  // dump global info; num documents in collection, num of all terms
  std::ofstream of_globalinfo(global_info_file);
  of_globalinfo << search_engine.document_count() << " "
                << search_engine.term_count() << std::endl;

  // write the lengths and names
  {
    std::cout << "Writing document lengths to " << doclen_tfile << "."
      << std::endl;
    std::cout << "Writing document names to " << doc_names_file << "." 
      << std::endl;
    std::ofstream of_doc_names(doc_names_file);
    
    long long start = search_engine.get_variable(ATIRE_DOCUMENT_FILE_START);
    long long end = search_engine.get_variable(ATIRE_DOCUMENT_FILE_END);
    unsigned long bsize = end - start;
    char *buffer = (char *)malloc(bsize);
    auto filenames = search_engine.get_document_filenames(buffer, &bsize);

    uint64_t uniq_terms = search_engine.get_unique_term_count();
    // Shift all IDs from ATIRE by 2 so \0 and \1 are free.
    uniq_terms += 2; 

    double mean_length;
    auto lengths = search_engine.get_document_lengths(&mean_length);
    {
      for (long long i = 0; i < search_engine.document_count(); i++)
      {
        doclen_out << lengths[i] << std::endl;
        of_doc_names << filenames[i] << std::endl;;
      }
    }

    free(buffer);
  }
  // write dictionary
  {
    std::cout << "Writing dictionary to " << dict_file << "." << std::endl;
    std::ofstream of_dict(dict_file);

    ANT_search_engine_btree_leaf leaf;
    ANT_btree_iterator iter(&search_engine);

    size_t j = 2;
    for (char *term = iter.first(NULL); term != NULL; term = iter.next()) {
      iter.get_postings_details(&leaf);
      of_dict << term << " " << j << " "
        << leaf.local_document_frequency << " "
        << leaf.local_collection_frequency << " "
        << "\n";
      map.emplace(strdup(term), j);
      j++;
    }
  }

  // write inverted files
  {
    using plist_type = block_postings_list<128>;
    vector<plist_type> m_postings_lists;
    vector<vector<pair<uint64_t, uint64_t>>> temp_postings_lists;
    uint64_t a = 0, b = 0;
    uint64_t n_terms = search_engine.get_unique_term_count() + INDRI_OFFSET; // + 2 to skip 0 and 1
 
    vector<pair<uint64_t, uint64_t>> post; 
    post.reserve(INIT_SZ);

    // Open the files
    filebuf post_file;
    post_file.open(postings_file, std::ios::out);
    ostream ofs(&post_file);

    std::cerr << "Generating postings lists ..." << std::endl;

    m_postings_lists.resize(n_terms);

    ANT_search_engine_btree_leaf leaf;
    ANT_btree_iterator iter(&search_engine);
    ANT_impact_header impact_header;
    ANT_compression_factory factory;

    ANT_compressable_integer *raw;
    long long impact_header_size = ANT_impact_header::NUM_OF_QUANTUMS * sizeof(ANT_compressable_integer) * 3;
    ANT_compressable_integer *impact_header_buffer = (ANT_compressable_integer *)malloc(impact_header_size);
    auto postings_list_size = search_engine.get_postings_buffer_length();
    auto raw_list_size = sizeof(*raw) * (search_engine.document_count() + ANT_COMPRESSION_FACTORY_END_PADDING);
    unsigned char *postings_list = (unsigned char *)malloc((size_t)postings_list_size);
    raw = (ANT_compressable_integer *)malloc((size_t)raw_list_size);
    uint64_t term_count = 0;

    size_t num_lists = n_terms;
    cout << "Writing " << num_lists << " postings lists." << endl;
    sdsl::serialize(num_lists, ofs);

    // take the 0 and 1 terms with dummies
    sdsl::serialize(block_postings_list<128>(), ofs);
    sdsl::serialize(block_postings_list<128>(), ofs);

     for (char *term = iter.first(NULL); term != NULL; term_count++, term = iter.next())
    {
	// don't capture ~ terms, they are specific to ATIRE
      if (*term == '~')
        break;

      iter.get_postings_details(&leaf);
      postings_list = search_engine.get_postings(&leaf, postings_list);

      auto the_quantum_count = ANT_impact_header::get_quantum_count(postings_list);
      auto beginning_of_the_postings = ANT_impact_header::get_beginning_of_the_postings(postings_list);
      factory.decompress(impact_header_buffer, postings_list + ANT_impact_header::INFO_SIZE, the_quantum_count * 3);

      if (term_count % 100000 == 0) {
      /* if (true) { */
        std::cout << term << " @ " << leaf.postings_position_on_disk << " (cf:" << leaf.local_collection_frequency << ", df:" << leaf.local_document_frequency << ", q:" << the_quantum_count << ")" << std::endl;
		fflush(stdout);
      }

      long long docid, max_docid, sum;
      ANT_compressable_integer *impact_header = (ANT_compressable_integer *)impact_header_buffer;
      ANT_compressable_integer *current, *end;

      max_docid = sum = 0;
      ANT_compressable_integer *impact_value_ptr = impact_header;
      ANT_compressable_integer *doc_count_ptr = impact_header + the_quantum_count;
      ANT_compressable_integer *impact_offset_start = impact_header + the_quantum_count * 2;
      ANT_compressable_integer *impact_offset_ptr = impact_offset_start;

      post.clear();
      post.reserve(leaf.local_document_frequency);


      while (doc_count_ptr < impact_offset_start) {
        factory.decompress(raw, postings_list + beginning_of_the_postings + *impact_offset_ptr, *doc_count_ptr);
        docid = -1;
        current = raw;
        end = raw + *doc_count_ptr;
        while (current < end) {
          docid += *current++;
          post.emplace_back(docid, *impact_value_ptr);
        }
        impact_value_ptr++;
        impact_offset_ptr++;
        doc_count_ptr++;
      }

      // The above will result in sorted by impact first, so re-sort by docid
      std::sort(std::begin(post), std::end(post));

      plist_type pl(post, index_format);
      sdsl::serialize(pl, ofs);

    }
    //close output files
    post_file.close();
  }

	auto build_stop = clock::now();
	auto build_time_sec = std::chrono::duration_cast<std::chrono::seconds>(build_stop-build_start);
	std::cout << "Index built in " << build_time_sec.count() << " seconds." << std::endl;

	return EXIT_SUCCESS;
}
