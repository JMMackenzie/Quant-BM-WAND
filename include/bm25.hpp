#ifndef _BM25_H
#define _BM25_H

#include "util.hpp"
#include "generic_rank.hpp"

struct rank_bm25 : public generic_rank {

  static constexpr double k1 = 0.9;
  static constexpr double b = 0.4;
  const double epsilon_score = 1e-6;
  size_t num_docs;
  size_t num_terms;
  double avg_doc_len;
  double min_doc_len;
  std::vector<uint64_t> doc_lengths;

  static std::string name() {
    return "bm25";
  }

  rank_bm25(){}

  rank_bm25& operator=(const rank_bm25&) = default;

  rank_bm25(std::vector<uint64_t> doc_len, 
          uint64_t terms) : rank_bm25(doc_len, terms, doc_len.size()) { }

  rank_bm25(std::vector<uint64_t> doc_len, 
          uint64_t terms, uint64_t numdocs) : num_docs(numdocs), 
          avg_doc_len((double)terms/(double)numdocs) {
    doc_lengths = std::move(doc_len); //Takes ownership of the vector!

    std::cerr<<"num_docs = "<<num_docs<<std::endl;
    std::cerr<<"avg_doc_len = "<<avg_doc_len<<std::endl;
  }

  uint64_t doc_length(const uint64_t doc_id) const {
    return doc_lengths[doc_id];
  }
  
  double calculate_docscore(const uint64_t f_dt,
                            const uint64_t f_t, const double W_d) const {
    double w_qt = std::max(epsilon_score, 
                  log((num_docs - (double)f_t + 0.5) / ((double)f_t+0.5)));
    double K_d = k1*((1-b) + (b*(W_d/avg_doc_len)));
    double w_dt = ((k1+1)*(double)f_dt) / (K_d + f_dt);
    return w_dt*w_qt;
  }

};

#endif
