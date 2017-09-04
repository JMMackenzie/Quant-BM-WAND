#ifndef QUERY_HPP
#define QUERY_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>

#include "util.hpp"


struct doc_score {
	uint64_t doc_id;
	double score;
  bool operator>(const doc_score& rhs) const {
  	if(score == rhs.score)
    	return doc_id > rhs.doc_id;
      return score > rhs.score;
    }
  doc_score() {};
  doc_score(uint64_t did, double s) : doc_id(did) , score(s) {};
};

struct result {
  std::vector<doc_score> list;
  uint64_t qry_id = 0;
  uint64_t wt_search_space = 0;
  uint64_t wt_nodes = 0;
  uint64_t postings_evaluated = 0; // Sum calls to "evaluate_pivot"
  uint64_t postings_total = 0; // Sum of lengths of postings lists from query
  uint64_t docs_fully_evaluated = 0;
  uint64_t docs_added_to_heap = 0; 
  double final_threshold = 0; // Final top-k heap threshold
};

struct query_token{
    uint64_t token_id;
    std::string token_str;
    uint64_t f_qt;
	query_token(const uint64_t id,
              const std::string str,
              uint64_t f) : token_id(id), token_str(str), 
              f_qt(f) 
    {
    }
};

using query_t = std::tuple<uint64_t,std::vector<query_token>>;


struct query_parser {
    query_parser() = delete;
    using mapping_t = std::pair<std::unordered_map<std::string,uint64_t>,
                     std::unordered_map<uint64_t,std::string>
                     >;

    static mapping_t
         load_dictionary(const std::string& collection_dir)
    {
        std::unordered_map<std::string,uint64_t> id_mapping;
        std::unordered_map<uint64_t,std::string> reverse_id_mapping;
        {
            auto dict_file = collection_dir + "/" + DICT_FILENAME;
            std::ifstream dfs(dict_file);
            if(!dfs.is_open()) {
                std::cerr << "cannot load dictionary file.";
                exit(EXIT_FAILURE);
            }
            std::string term_mapping;
            while( std::getline(dfs,term_mapping) ) {
                auto sep_pos = term_mapping.find(' ');
                auto term = term_mapping.substr(0,sep_pos);
                auto idstr = term_mapping.substr(sep_pos+1);
                uint64_t id = std::stoull(idstr);
                id_mapping[term] = id;
                reverse_id_mapping[id] = term;
            }
        }
        return {id_mapping,reverse_id_mapping};
    }

    static std::tuple<bool,uint64_t,std::vector<uint64_t>> 
        map_to_ids(const std::unordered_map<std::string,uint64_t>& id_mapping,
                   std::string query_str,bool only_complete,bool integers)
    {
        auto id_sep_pos = query_str.find(';');
        auto qryid_str = query_str.substr(0,id_sep_pos);
        auto qry_id = std::stoull(qryid_str);
        auto qry_content = query_str.substr(id_sep_pos+1);

        std::vector<uint64_t> ids;
        std::istringstream qry_content_stream(qry_content);
        for(std::string qry_token; std::getline(qry_content_stream,qry_token,' ');) {
            if(integers) {
                uint64_t id = std::stoull(qry_token);
                ids.push_back(id);
            } else {
                auto id_itr = id_mapping.find(qry_token);
                if(id_itr != id_mapping.end()) {
                    ids.push_back(id_itr->second);
                } else {
                    std::cerr << "ERROR: could not find '" 
                              << qry_token << "' in the dictionary." 
                              << std::endl;
                    if(only_complete) {
                        return std::make_tuple(false,qry_id,ids);
                    }
                }
            }
        }
        return std::make_tuple(true,qry_id,ids);
    }

    static std::pair<bool,query_t> parse_query(const mapping_t& mapping,
                const std::string& query_str,
                bool only_complete = false,bool integers = false)
    {

        const auto& id_mapping = mapping.first;
        const auto& reverse_mapping = mapping.second;

        auto mapped_qry = map_to_ids(id_mapping,query_str,only_complete,integers);

        bool parse_ok = std::get<0>(mapped_qry);
        auto qry_id = std::get<1>(mapped_qry);
        
        if(parse_ok) {
            std::unordered_map<uint64_t,uint64_t> qry_set;
            const auto& tids = std::get<2>(mapped_qry);
            for(const auto& tid : tids) {
                qry_set[tid] += 1;
            }
            std::vector<query_token> query_tokens;
            size_t index = 0;
            for(const auto& qry_tok : qry_set) {
                uint64_t term = qry_tok.first;
                auto rmitr = reverse_mapping.find(term);
                std::string term_str;
                if(rmitr != reverse_mapping.end()) {
                    term_str = rmitr->second;
                }
                query_tokens.emplace_back(term,term_str,qry_tok.second);
                ++index;
            }
            query_t q(qry_id,query_tokens);
            return {true,q};
        }

        // error
        query_t q;
        return {false,q};
    }

    static std::vector<query_t> parse_queries(const std::string& collection_dir,
                                              const std::string& query_file,
                                              bool only_complete = false) {
        std::vector<query_t> queries;

        /* load the mapping */
        auto mapping = load_dictionary(collection_dir);
        /* parse queries */
        std::ifstream qfs(query_file); 
        if(!qfs.is_open()) {
            std::cerr << "cannot load query file.";
            exit(EXIT_FAILURE);
        }

        std::string query_str;
        while( std::getline(qfs,query_str) ) {
            auto parsed_qry = parse_query(mapping,query_str);
            if(parsed_qry.first) {
                queries.emplace_back(parsed_qry.second);
            }
        }

        return queries;
    }
};

#endif
