#ifndef BLOCK_POSTINGS_LIST_H
#define BLOCK_POSTINGS_LIST_H

#include <limits>
#include <stdexcept>
#include <x86intrin.h>

#include "util.hpp"
#include "util.h"
#include "memutil.h"
#include "codecs.h"
#include "codecfactory.h"
#include "bitpacking.h"
#include "simdfastpfor.h"
#include "deltautil.h"

#include "compress_qmx.h"

#include "sdsl/int_vector.hpp"

using namespace sdsl;


template<uint64_t t_block_size>
class block_postings_list;

template<uint64_t t_block_size>
class plist_iterator
{
  public:
    typedef block_postings_list<t_block_size> list_type;
    typedef typename list_type::size_type     size_type;
    typedef uint64_t                          value_type;
  public: // default implementation used. not necessary to list here
    plist_iterator() = default;
    plist_iterator(const plist_iterator& pi) = default;
    plist_iterator(plist_iterator&& pi) = default;
    plist_iterator& operator=(const plist_iterator& pi) = default;
    plist_iterator& operator=(plist_iterator&& pi) = default;
  public:
    plist_iterator(const list_type& l,size_t pos);
    plist_iterator& operator++();
    bool operator ==(const plist_iterator& b) const;
    bool operator !=(const plist_iterator& b) const;
    uint64_t docid() const;
    uint64_t freq() const;
    void skip_to_id(const uint64_t id);
    void skip_to_block_with_id(const uint64_t id);
    uint64_t block_max() const;
    uint64_t block_max(const uint64_t id) const {
      return m_plist_ptr->m_block_maximums[id];
    }
    uint64_t block_rep() const { 
      return m_plist_ptr->block_rep(m_cur_block_id); 
    }
    uint64_t block_rep(const uint64_t id) const {
      return m_plist_ptr->block_rep(id);
    }
    const uint64_t block_containing_id(const uint64_t id);
    uint64_t num_blocks() const {
      return m_plist_ptr->num_blocks();
    }
    size_t size() const { return m_plist_ptr->size(); }
    size_t remaining() const { return size() - m_cur_pos; }
    size_t offset() const { return m_cur_pos; }
  private:
    void access_and_decode_cur_pos() const;
  private:
    size_type m_cur_pos = std::numeric_limits<uint64_t>::max();
    mutable size_type m_cur_block_id = std::numeric_limits<uint64_t>::max();
    mutable size_type m_last_accessed_block = 
            std::numeric_limits<uint64_t>::max()-1;
    mutable size_type m_last_accessed_id = 
            std::numeric_limits<uint64_t>::max()-1;
    mutable value_type m_cur_docid = 0;
    mutable value_type m_cur_freq = 0;
    const list_type* m_plist_ptr = nullptr;
    mutable std::vector<uint32_t, FastPForLib::cacheallocator> m_decoded_ids;
    mutable std::vector<uint32_t, FastPForLib::cacheallocator> m_decoded_freqs;
};

template<uint64_t t_block_size=128>
class block_postings_list {
	static_assert(t_block_size % 32 == 0,"blocksize must be multiple of 32.");
  public: // types
	  friend class plist_iterator<t_block_size>;
	  using comp_codec = ANT_compress_qmx;
	  using freq_codec = ANT_compress_qmx;
	  using size_type = sdsl::int_vector<>::size_type;
	  using const_iterator = plist_iterator<t_block_size>;
	  using pfor_data_type = std::vector<uint32_t, FastPForLib::cacheallocator>;
	  #pragma pack(push, 1)
	  struct block_data {
		  uint32_t max_block_id = 0;
		  uint32_t id_offset = 0;
		  uint32_t freq_offset = 0;
		  uint32_t id_bytes = -1;
		  uint32_t freq_bytes = -1;
	  };
	  #pragma pack(pop)
  public: // actual data
	  uint64_t m_size = 0;
	  uint64_t m_list_maximum = std::numeric_limits<uint64_t>::lowest();
	  std::vector<block_data> m_block_data;
    pfor_data_type m_docid_data;
    pfor_data_type m_freq_data;
    std::vector<uint64_t> m_block_maximums;
  public: // default 
    block_postings_list() {
    	m_block_data.resize(1);
    }
    uint64_t block_max(uint64_t bid) const {
      return m_block_maximums[bid];
    }

    block_postings_list(const block_postings_list& pl) = default;
    block_postings_list(block_postings_list&& pl) = default;
    block_postings_list& operator=(const block_postings_list& pi) = default;
    block_postings_list& operator=(block_postings_list&& pi) = default;
    uint64_t list_max_score() const { return m_list_maximum; };
public: // constructors
    block_postings_list(std::istream& in) {
      load(in);
    }

 
    block_postings_list(std::vector<std::pair<uint64_t,uint64_t>>& pre_sorted_data,
                       index_form index_type) {

    	m_size = pre_sorted_data.size();

	    // extract doc_ids and freqs
	    sdsl::int_vector<32> tmp_data(pre_sorted_data.size());
	    sdsl::int_vector<32> tmp_freq(pre_sorted_data.size());
	    for (size_type i=0; i<pre_sorted_data.size(); i++) {
	        tmp_data[i] = pre_sorted_data[i].first;
	        tmp_freq[i] = pre_sorted_data[i].second;
	    }

	    // create block max structure first
	    create_block_support(tmp_data);

  
	    // create rank support structure
	
      if (index_type == BMW) {
        // BMW specific
        create_rank_support_bmw(tmp_data,tmp_freq);
      }
      else {
        // Wand specific
        create_rank_support_wand(tmp_data,tmp_freq);
      }

	    // compress postings
	    compress_postings_data(tmp_data,tmp_freq);
    }
  
  private: // functions used during construction
	  void create_block_support(const sdsl::int_vector<32>& ids)
	  {
	    size_t num_blocks = ids.size() / t_block_size;
	    if (ids.size() % t_block_size != 0) num_blocks++;
	    m_block_data.resize(num_blocks);
	    size_t j = 0;
	    for (size_t i=t_block_size-1; i<ids.size(); i+=t_block_size) {
	      m_block_data[j++].max_block_id = ids[i];
	    }
	    if (ids.size() % t_block_size != 0) {
        m_block_data[j].max_block_id = ids[ids.size()-1];
      }
	  }

	  void create_rank_support_wand(const sdsl::int_vector<32>& ids,
							               const sdsl::int_vector<32>& freqs)
	  {
		  auto F_t = std::accumulate(freqs.begin(),freqs.end(),0);
		  auto f_t = ids.size();
      
      uint64_t max_score = 0;
      m_list_maximum = std::numeric_limits<uint64_t>::lowest();
  
	    for (size_t l=0; l<ids.size(); l++) {
	      auto id = ids[l];
	      uint64_t score = freqs[l];
	      max_score = std::max(max_score, score);
      }
      m_list_maximum = std::max(m_list_maximum, max_score);
	  }


	  void create_rank_support_bmw(const sdsl::int_vector<32>& ids,
							               const sdsl::int_vector<32>& freqs)
	  {
		  auto F_t = std::accumulate(freqs.begin(),freqs.end(),0);
		  auto f_t = ids.size();
      
      size_t num_blocks = ids.size() / t_block_size;
      if(ids.size() % t_block_size != 0)
        num_blocks++;

      m_block_maximums.resize(num_blocks);
      uint64_t max_score = 0;
      m_list_maximum = std::numeric_limits<uint64_t>::lowest();
      size_t i = 1;
      size_t j = 0;
  
	    for (size_t l=0; l<ids.size(); l++) {
	      auto id = ids[l];
	      uint64_t score = freqs[l];
	      max_score = std::max(max_score, score);
        //Block max support
        if(i % t_block_size == 0){
          m_block_maximums[j] = max_score;
          m_list_maximum = std::max(m_list_maximum, max_score);
          i = 0;
          max_score = 0;
          j++;
        }
        i++;
      }
      if (ids.size() % t_block_size != 0){
        m_block_maximums[num_blocks-1] = max_score;
      }
      m_list_maximum = std::max(m_list_maximum, max_score);
	  }

	  void compress_postings_data(const sdsl::int_vector<32>& ids,
	            					        sdsl::int_vector<32>& freqs)
	  {
		  static comp_codec c;
		  static freq_codec fc;

		  uint32_t *id_input = (uint32_t *)ids.data();
		  FastPForLib::Delta::fastDelta(id_input,ids.size());
		  uint32_t *freq_input = (uint32_t *)freqs.data();

		  m_docid_data.resize(2 * ids.size() + 1024);
		  m_freq_data.resize(2 * freqs.size() + 1024);

		  uint32_t *id_out = m_docid_data.data();
		  uint32_t *freq_out = m_freq_data.data();

		  size_type cur_block = 0;
		  uint64_t id_offset = 0;
		  uint64_t freq_offset = 0;
		  size_t encoded_id_size = 0;
		  size_t encoded_freq_size = 0;
		  uint64_t bytes_used = 0;
		  uint64_t freq_bytes_used = 0;

		  size_type n = t_block_size;
		  for (size_t i = 0; i < ids.size(); i+=t_block_size) {
			  if (i + t_block_size > ids.size())
				  n = ids.size() % t_block_size;

			  m_block_data[cur_block].id_offset = id_offset;
			  m_block_data[cur_block].freq_offset = freq_offset;
			  c.encodeArray(&id_input[i], n, &id_out[id_offset], &bytes_used);
			  fc.encodeArray(&freq_input[i], n, &freq_out[freq_offset], &freq_bytes_used);

			  id_offset += (bytes_used / sizeof(uint32_t));
			  freq_offset += (freq_bytes_used / sizeof(uint32_t));

			  if (bytes_used % sizeof(uint32_t) != 0){
				  id_offset++;
			  }
			  if (freq_bytes_used % sizeof(uint32_t) != 0){
				  freq_offset++;
			  }
			  auto alignment = sizeof(uint32_t);
			  if (id_offset % alignment != 0) {
				  id_offset += alignment - (id_offset % alignment);
			  }
			  if (freq_offset % alignment != 0) {
				  freq_offset += alignment - (freq_offset % alignment);
			  }
			  if (id_offset > m_docid_data.size() || freq_offset > m_freq_data.size()) {
				  std::cerr << "Run out of room encoding!" << std::endl;
				  exit(EXIT_FAILURE);
			  }
			  m_block_data[cur_block].id_bytes = bytes_used;
			  m_block_data[cur_block].freq_bytes = freq_bytes_used;

			  cur_block++;
		  }

		  m_docid_data.resize(id_offset);
		  m_docid_data.shrink_to_fit();
		  m_freq_data.resize(freq_offset);
		  m_freq_data.shrink_to_fit();
	  }
  public: // functions used during processing
	  
    void decompress_block(const size_t block_id,
	            					  pfor_data_type& id_data,
						              pfor_data_type& freq_data) const
	{
		static comp_codec c;
		static freq_codec fc;

		uint32_t delta_offset = 0;
		if (block_id != 0) {
			delta_offset = m_block_data[block_id - 1].max_block_id;
		}

		const uint32_t *id_start = m_docid_data.data() + m_block_data[block_id].id_offset;
		const uint32_t *freq_start = m_freq_data.data() + m_block_data[block_id].freq_offset;
		auto block_size = postings_in_block(block_id);

		if (id_data.size() != block_size) {
			id_data.resize(block_size);
			freq_data.resize(block_size);
		}

		c.decodeArray(id_start, m_block_data[block_id].id_bytes, id_data.data(), block_size);

		/* Extracted from: https:github.com/lemire/FastDifferentialCoding */
		__m128i prev = _mm_set1_epi32(delta_offset);
		size_t i = 0;
		for (; i  < block_size/4; i++) {
			__m128i curr = _mm_lddqu_si128((const __m128i *)id_data.data() + i);
			const __m128i _tmp1 = _mm_add_epi32(_mm_slli_si128(curr, 8), curr);
			const __m128i _tmp2 = _mm_add_epi32(_mm_slli_si128(_tmp1, 4), _tmp1);
			prev = _mm_add_epi32(_tmp2, _mm_shuffle_epi32(prev, 0xff));
			_mm_storeu_si128((__m128i *)id_data.data() + i,prev);
		}
		uint32_t lastprev = _mm_extract_epi32(prev, 3);
		for(i = 4 * i ; i < block_size; ++i) {
			lastprev = lastprev + id_data[i];
			id_data[i] = lastprev;
		}

		fc.decodeArray(freq_start, m_block_data[block_id].freq_bytes, freq_data.data(), block_size);
	}

	  size_type find_block_with_id(const uint64_t id, const size_t start_block) const {
	    size_t block_id = start_block;
	    size_t nblocks = m_block_data.size();
	    while (block_id < nblocks && m_block_data[block_id].max_block_id < id) {
	      block_id++;
	    }
	    return block_id;
	  }

	  size_type size() const {
		  return m_size;
	  }

	  uint64_t block_rep(const size_t bid) const {
		  return m_block_data[bid].max_block_id;
	  }

	  size_type num_blocks() const {
		  return m_block_data.size();
	  }

	  size_type postings_in_block(const size_type block_id) const {
		  size_type block_size = t_block_size;
		  size_type mod = m_size % t_block_size;
		  if (block_id == m_block_data.size()-1 && mod != 0) {
			  block_size = mod;
		  }
		  return block_size;
	  }

    const_iterator begin() const {
      return const_iterator(*this,0);
    }

    const_iterator end() const {
      return const_iterator(*this,m_size);
    }

    auto serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, 
    			         std::string name = "") const -> size_type 
	  {

	    size_type written_bytes = 0;

		  sdsl::structure_tree_node* child;
	    if (m_size <= t_block_size) { // only one block
        //std::cout << "single" << std::endl;
	    	child = sdsl::structure_tree::add_child(v,"single block list",
                                                sdsl::util::class_name(*this));
	    } else {
        //std::cout << "multi" << std::endl;
	    	child = sdsl::structure_tree::add_child(v, "multi block list",
                                                sdsl::util::class_name(*this));
	    }

	    written_bytes += sdsl::write_member(m_size,out,child,"size");

	    if (m_size <= t_block_size) { // only one block
	     	written_bytes += sdsl::write_member(m_block_data[0].max_block_id,out,
                                            child,"max block id");
			written_bytes += sdsl::write_member(m_block_data[0].id_bytes, out, child, "id bytes used");
			written_bytes += sdsl::write_member(m_block_data[0].freq_bytes, out, child, "freq bytes used");
	    } else {
	    	auto* blockdata = sdsl::structure_tree::add_child(child, "block data",
                                                          "block data");
	    	out.write((const char*)m_block_data.data(), 
                  m_block_data.size()*sizeof(block_data));
	    	written_bytes += m_block_data.size()*sizeof(block_data);
	    	sdsl::structure_tree::add_size(blockdata, 
                                       m_block_data.size()*sizeof(block_data));
	    }

      uint32_t docidu32 = m_docid_data.size();
      uint32_t frequ32 = m_freq_data.size();
      written_bytes += sdsl::write_member(docidu32,out,child,"docid u32s");
      written_bytes += sdsl::write_member(frequ32,out,child,"freq u32s");

    	auto* idchild = sdsl::structure_tree::add_child(child, "id data",
                                                      "delta compressed");
      out.write((const char*)m_docid_data.data(), 
                m_docid_data.size()*sizeof(uint32_t));
      sdsl::structure_tree::add_size(idchild, 
                                     m_docid_data.size()*sizeof(uint32_t));
      written_bytes +=  m_docid_data.size()*sizeof(uint32_t);

      auto* fchild = sdsl::structure_tree::add_child(child, "freq data", 
                                                     "compressed");
      out.write((const char*)m_freq_data.data(), 
                m_freq_data.size()*sizeof(uint32_t));
      written_bytes +=  m_freq_data.size()*sizeof(uint32_t);
    	sdsl::structure_tree::add_size(fchild, 
                                     m_freq_data.size()*sizeof(uint32_t));
     
      //Write blockmax data
      auto *bm_c = sdsl::structure_tree::add_child(child, "blockmax", 
                                                                    "blockmax");
      auto *bmm_c = sdsl::structure_tree::add_child(bm_c, "block maximums",
                                   "uint64_t");
      size_type bm_written_bytes = sdsl::write_member(m_block_maximums.size(), 
                                  out, bm_c, "num max_scores");

      out.write((const char *)m_block_maximums.data(), 
                                      m_block_maximums.size() * sizeof(uint64_t));
      bm_written_bytes += m_block_maximums.size() * sizeof(uint64_t);
      sdsl::structure_tree::add_size(bmm_c, 
                                   m_block_maximums.size() * sizeof(uint64_t));
      sdsl::structure_tree::add_size(bm_c, bm_written_bytes);
      written_bytes += bm_written_bytes; 
      

	    written_bytes += sdsl::write_member(m_list_maximum,out,
                                          child,"list max score");

	    sdsl::structure_tree::add_size(child, written_bytes);
	    return written_bytes;
	  }

	  void load(std::istream& in) {
		  read_member(m_size,in);
		  if (m_size <= t_block_size) { // only one block
			  uint32_t max_block_id;
			  uint32_t id_bytes_used;
			  uint32_t freq_bytes_used;
			  read_member(max_block_id,in);
			  read_member(id_bytes_used,in);
			  read_member(freq_bytes_used,in);
			  m_block_data.resize(1);
			  m_block_data[0].max_block_id = max_block_id;
			  m_block_data[0].id_bytes = id_bytes_used;
			  m_block_data[0].freq_bytes = freq_bytes_used;
		  } else {
			  uint64_t num_blocks = m_size / t_block_size;
			  if (m_size % t_block_size != 0) num_blocks++;
			  m_block_data.resize(num_blocks);
			  in.read((char*)m_block_data.data(),num_blocks*sizeof(block_data));
		  }

		  // load compressed data
      uint32_t docidu32;
      uint32_t frequ32;
      read_member(docidu32,in);
      read_member(frequ32,in);
      m_docid_data.resize(docidu32);
      m_freq_data.resize(frequ32);
      in.read((char*)m_docid_data.data(),docidu32*sizeof(uint32_t));
      in.read((char*)m_freq_data.data(),frequ32*sizeof(uint32_t));
      size_t num_block_max_scores;
      read_member(num_block_max_scores, in);
      m_block_maximums.resize(num_block_max_scores);
      in.read((char *)m_block_maximums.data(), 
                                         num_block_max_scores * sizeof(uint64_t));

      read_member(m_list_maximum,in);
	}
};


template<uint64_t t_bs>
plist_iterator<t_bs>::plist_iterator(const list_type& l,
                                     size_t pos) : plist_iterator()
{
  m_cur_pos = pos;
  m_plist_ptr = &l;
}

template<uint64_t t_bs>
plist_iterator<t_bs>& plist_iterator<t_bs>::operator++()
{
  if (m_cur_pos != size()) { // end?
    (*this).m_cur_pos++;
  } else {
    std::cerr << "ERROR: trying to advance plist iterator beyond list end.\n";
    throw std::out_of_range("trying to advance plist iterator beyond list end");
  }
  return (*this);
}

template<uint64_t t_bs>
bool plist_iterator<t_bs>::operator ==(const plist_iterator& b) const
{
  return ((*this).m_cur_pos == b.m_cur_pos) && 
          ((*this).m_plist_ptr == b.m_plist_ptr);
}

template<uint64_t t_bs>
bool plist_iterator<t_bs>::operator !=(const plist_iterator& b) const
{
  return !((*this)==b);
}

template<uint64_t t_bs>
typename plist_iterator<t_bs>::value_type plist_iterator<t_bs>::docid() const
{
  if (m_cur_pos == m_plist_ptr->size()) { // end?
    std::cerr << "ERROR: plist iterator dereferenced at list end.\n";
    throw std::out_of_range("plist iterator dereferenced at list end");
  }
  if (m_cur_pos == m_last_accessed_id) {
    return m_cur_docid;
  }
  access_and_decode_cur_pos();
  return m_cur_docid;
}

template<uint64_t t_bs>
uint64_t plist_iterator<t_bs>::block_max() const 
{
  return m_plist_ptr->block_max(m_cur_block_id);
}

template<uint64_t t_bs>
typename plist_iterator<t_bs>::value_type plist_iterator<t_bs>::freq() const
{
  if (m_cur_pos == m_plist_ptr->size()) { // end?
    std::cerr << "ERROR: plist iterator dereferenced at list end.\n";
    throw std::out_of_range("plist iterator dereferenced at list end");
  }
  if (m_cur_pos == m_last_accessed_id) {
    return m_cur_freq;
  }
  access_and_decode_cur_pos();
  return m_cur_freq;
}

template<uint64_t t_bs>
void plist_iterator<t_bs>::access_and_decode_cur_pos() const
{
  m_cur_block_id = m_cur_pos / t_bs;
  if (m_cur_block_id != m_last_accessed_block) {  // decompress block
    m_last_accessed_block = m_cur_block_id;
    m_plist_ptr->decompress_block(m_cur_block_id,m_decoded_ids,m_decoded_freqs);
  }
  size_t in_block_offset = m_cur_pos % t_bs;
  m_cur_docid = m_decoded_ids[in_block_offset];
  m_cur_freq = m_decoded_freqs[in_block_offset];
  m_last_accessed_id = m_cur_pos;
}

template<uint64_t t_bs>
const uint64_t plist_iterator<t_bs>::block_containing_id(const uint64_t id) {
  size_t block = m_plist_ptr->find_block_with_id(id, m_cur_block_id);
  m_cur_block_id = block;
  return block;
} 

template<uint64_t t_bs>
void plist_iterator<t_bs>::skip_to_block_with_id(const uint64_t id)
{
  size_t old_block = m_cur_block_id;
  m_cur_block_id = m_plist_ptr->find_block_with_id(id,m_cur_block_id);

  // we now go to the first id in the new block!
  if (old_block != m_cur_block_id) {
    m_cur_pos = m_cur_block_id*t_bs;
    if (m_cur_pos > m_plist_ptr->size()) { // don't go past the end!
      m_cur_pos = m_plist_ptr->size();
    }
  }
}

template<uint64_t t_bs>
void plist_iterator<t_bs>::skip_to_id(const uint64_t id)
{
  if (id == m_cur_docid) {
    return;
  }

  skip_to_block_with_id(id);
  // check if we reached list end!
  if (m_cur_block_id >= m_plist_ptr->num_blocks()) {
    m_cur_pos = m_plist_ptr->size();
    return;
  }
  if (m_last_accessed_block != m_cur_block_id) {
    m_last_accessed_block = m_cur_block_id;
    m_plist_ptr->decompress_block(m_cur_block_id,m_decoded_ids,m_decoded_freqs);
    auto block_itr = std::lower_bound(m_decoded_ids.begin(),
                                      m_decoded_ids.end(),id);
    m_cur_pos = (t_bs*m_cur_block_id) + 
                std::distance(m_decoded_ids.begin(),block_itr);
  } else {
    size_t in_block_offset = m_cur_pos % t_bs;
    auto block_itr = std::lower_bound(m_decoded_ids.begin()+in_block_offset,
                                      m_decoded_ids.end(),id);
    m_cur_pos = (t_bs*m_cur_block_id) + 
                std::distance(m_decoded_ids.begin(),block_itr);
  }
  size_t inblock_offset = m_cur_pos % t_bs;
  m_cur_docid = m_decoded_ids[inblock_offset];
  m_cur_freq = m_decoded_freqs[inblock_offset];
  m_last_accessed_id = m_cur_pos;
}

#endif
