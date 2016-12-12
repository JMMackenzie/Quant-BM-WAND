#ifndef _IMPACT_H
#define _IMPACT_H

#include <vector>

struct my_rank_impact {
	my_rank_impact(){}
	my_rank_impact& operator=(const my_rank_impact&) = default;
	my_rank_impact(std::vector<uint64_t> doc_len, uint64_t terms) { }
	my_rank_impact(std::vector<uint64_t> doc_len, uint64_t terms, uint64_t numdocs){ }

	static std::string name() {
		return "impact";
	}
	double doc_length(size_t ) const {
		return 0;
	}
	double calc_doc_weight(double ) const {
		return 0;
	}
	double calculate_docscore(const double , const double f_dt, const double , const double , bool) const {
		return f_dt;
	}
};

#endif
