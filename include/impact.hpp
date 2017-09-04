#ifndef _IMPACT_H
#define _IMPACT_H

#include "util.hpp"
#include "generic_rank.hpp"

struct rank_impact : public generic_rank {

	rank_impact(){}
	rank_impact& operator=(const rank_impact&) = default;
	rank_impact(std::vector<uint64_t> doc_len, uint64_t terms) { }
	rank_impact(std::vector<uint64_t> doc_len, uint64_t terms, uint64_t numdocs){ }

	static std::string name() {
		return "impact";
	}
	uint64_t doc_length(const uint64_t) const {
		return 0;
	}
	double calc_doc_weight(double ) const {
		return 0;
	}
	double calculate_docscore(const uint64_t f_dt, const uint64_t, const double) const {
		return f_dt;
	}
};

#endif
