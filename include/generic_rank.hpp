#ifndef _GENERIC_RANK_H
#define _GENERIC_RANK_H

#include "util.hpp"

struct generic_rank {
  virtual uint64_t doc_length(const uint64_t docid) const = 0;
	
  virtual double calculate_docscore(const uint64_t f_dt, 
                                    const uint64_t f_t, 
                                    const double w_d) const = 0;
  
  virtual ~generic_rank() {} // Avoid memory leaks: need virtual destruction
};

#endif
