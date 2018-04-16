/****
AC-DIAMOND: DNA-protein alignment tool
Copyright (C) 2018 Huijun Mai

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****/

#ifndef MASTER_THREAD_H_
#define MASTER_THREAD_H_

#include <iostream>
#include <boost/timer/timer.hpp>
#include "../data/reference.h"
#include "../data/queries.h"
#include "../basic/statistics.h"
#include "../output/join_blocks.h"
#include "../align/align_queries.h"
#include "../search/align.h"
#include "../basic/setup.h"
#include "../util/tinythread.h"
#include <pthread.h>
#include <string>
#include <algorithm>

using boost::timer::cpu_timer;
using boost::ptr_vector;

template<typename _locr, typename _locl>
bool cmp_locs(const pair<_locr, _locl> &lhs, const pair<_locr, _locl> &rhs){
  const uint64_t x = (uint64_t)lhs.first + (uint64_t)rhs.second;
  const uint64_t y = (uint64_t)rhs.first + (uint64_t)lhs.second;
  return (x < y) || (x == y && lhs.second < rhs.second);
}

unsigned binary_search(unsigned start, uint16_t cnt, uint16_t value, uint16_t *key16){
  if (cnt == 1)
    return key16[start] == value ? start : 0;
  
  int64_t mid, end;
  end = start + cnt - 1;
  
  while ((int64_t)start <= end){
    mid = (start + end) / 2;
    if (key16[mid] == value) return mid;
    else if (key16[mid] < value) start = mid + 1;
    else end = mid - 1;
  }
  
  return 0;
}

template<typename _val, typename _locr>
unsigned binary_search_upper(unsigned start, unsigned cnt, uint8_t value, unsigned pos, unsigned *sa, _locr part_len){
  int64_t s, e, mid;
  _val *subject;
  uint8_t ref_value;
  s = start;
  e = start + cnt;
  
  while (s < e){
    mid = (s + e) / 2;
    subject = ref_seqs<_val>::data_->data(part_len + sa[mid]);
    ref_value = Reduction<Amino_acid>::reduction((char)(*(subject+pos))) + 1;
    
    if (ref_value <= value) s = mid + 1;
    else e = mid;
  }
  
  return (unsigned)s;
}

template<typename _val, typename _locr>
unsigned binary_search_lower(unsigned start, unsigned cnt, uint8_t value, unsigned pos, unsigned *sa, _locr part_len){
  int64_t s, e, mid;
  _val *subject;
  uint8_t ref_value;
  s = start;
  e = start + cnt;

  while (s < e){
    mid = (s + e) / 2;
    subject = ref_seqs<_val>::data_->data(part_len + sa[mid]);
    ref_value = Reduction<Amino_acid>::reduction((char)(*(subject+pos))) + 1;
    
    if (ref_value >= value) e = mid;
    else s = mid + 1;
  }
  
  return (unsigned)s;
}

template<typename _val, typename _locr>
_locr binary_search_upper_long(_locr start, unsigned cnt, uint8_t value, unsigned pos, unsigned *sa_4bytes, uint8_t *sa_1byte){
  _locr s, e, mid, loc;
  _val *subject;
  uint8_t ref_value;
  s = start;
  e = start + cnt;

  while (s < e){
    mid = (s + e) / 2;
    loc = ((_locr)sa_4bytes[mid] << 8) | sa_1byte[mid];
    subject = ref_seqs<_val>::data_->data(loc);
    ref_value = Reduction<Amino_acid>::reduction((char)(*(subject+pos))) + 1;
    
    if (ref_value <= value) s = mid + 1;
    else e = mid;
  }
  
  return s;
}

template<typename _val, typename _locr>
_locr binary_search_lower_long(_locr start, unsigned cnt, uint8_t value, unsigned pos, unsigned *sa_4bytes, uint8_t *sa_1byte){
  _locr s, e, mid, loc;
  _val *subject;
  uint8_t ref_value;
  s = start;
  e = start + cnt;
  
  while (s < e){
    mid = (s + e) / 2;
    loc = ((_locr)sa_4bytes[mid] << 8) | sa_1byte[mid];
    subject = ref_seqs<_val>::data_->data(loc);
    ref_value = Reduction<Amino_acid>::reduction((char)(*(subject+pos))) + 1;

    if (ref_value >= value) e = mid;
    else s = mid + 1;
  }
  
  return s;
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void process(unsigned *sa,
	     unsigned *table,
	     uint16_t *count,
	     unsigned *ref_sal,
	     unsigned *ref_sar,
	     uint16_t *key16,
	     uint64_t mis_pattern,
	     vector<unsigned> &split,
	     unsigned partition,
	     unsigned max_seq_len,
	     _locr part_len,
	     unsigned thread_id)
{  
  typename Trace_pt_buffer<_locr,_locl>::Iterator* out = new typename Trace_pt_buffer<_locr,_locl>::Iterator (*Trace_pt_buffer<_locr,_locl>::instance, thread_id);

  const Sequence_set<_val> *query = query_seqs<_val>::data_;
  unsigned hashValue, m, n, i, seq_len, k1, k2, id, l, lower_bound, upper_bound, index;
  unsigned r, cnt_24bits, band;
  uint16_t value_16bit;
  uint64_t pattern_id64, pattern_id64_mask, pattern_id40, tmp;
  uint8_t value_8bit;
  _locr loc;
  _locl j;
  int score, k;
  _locr ref_loc;
  _locl query_loc;

  unsigned candidates;
  if (program_options::aligner_mode == program_options::sensitive){
    candidates = CANDIDATES_SENSITIVE;
  }
  else {
    candidates = CANDIDATES_FAST;
  }

  // added
  //candidates = 4096;

  vector<pair<_locr, _locl> > locs;
  
  uint64_t base[MAX_SEED_LEN];
  calculate_factor(base);
  
  uint8_t* query_reduced_aa = new uint8_t[max_seq_len];

  for (i = split[partition]; i < split[partition+1]; i++){
    const sequence<const _val> seq = (*query)[i];
    seq_len = seq.length();
    if (seq_len < MAX_SEED_LEN) continue;
    
    memset(query_reduced_aa, 0, sizeof(uint8_t) * max_seq_len);
    locs.clear();
    j = 0;
    
    while (j < seq_len - 1){
      if ((int)seq[j] == (int)Value_traits<Amino_acid>::MASK_CHAR){
	if ((int)seq[j+1] == (int)Value_traits<Amino_acid>::MASK_CHAR){
	  query_reduced_aa[j/2] = 0;
	}
	else {
	  k2 = Reduction<Amino_acid>::reduction((char)seq[j+1]) + 1;
	  query_reduced_aa[j/2] = k2 << 4;
	}
      }
      else {
	k1 = Reduction<Amino_acid>::reduction((char)seq[j]) + 1;
	if ((int)seq[j+1] == (int)Value_traits<Amino_acid>::MASK_CHAR){
	  query_reduced_aa[j/2] = k1;
	}
	else {
	  k2 = Reduction<Amino_acid>::reduction((char)seq[j+1]) + 1;
	  query_reduced_aa[j/2] =  (k2 << 4) | k1;
	}
      }
      
      j += 2;
    }
    
    if ((seq_len & 0x01) != 0){
      if ((int)seq[j] == (int)Value_traits<Amino_acid>::MASK_CHAR){
	query_reduced_aa[j/2] = 0;
      }
      else {
	k1 = Reduction<Amino_acid>::reduction((char)seq[j]) + 1;
	query_reduced_aa[j/2] = k1;
      }
    }

    j = 0;
    seq_len = seq_len - MAX_SEED_LEN + 1;
    
    while (j < seq_len){
      // get the pattern from each position of seq, except the last SEED_LEN-1 positions
      // for each pattern
      id = 0;
      
      if ((j & 0x01) == 0){
        pattern_id64 = (uint64_t)*((uint64_t*)(query_reduced_aa+j/2));
      }
      else {
        pattern_id64 = ((uint64_t)*((uint64_t*)(query_reduced_aa+j/2)) & 0xfffffffffffffff0) >> 4;
	tmp = ((uint64_t)*((uint64_t*)(query_reduced_aa+(j+1)/2)) & 0x0f00000000000000) << 4;
	pattern_id64 = pattern_id64 | tmp;
      }
      
      for (k = MAX_SEED_LEN-1; k >= 0; k--){
	if ((pattern_id64 & base[k]) == 0){
	  id = 1;
	  j = j + k + 1;
	  break;
	}
      }
      
      if (id == 1) continue;
      
      pattern_id64_mask = pattern_id64 & mis_pattern;
      pattern_id40 = pattern_id64_mask >> BS_BITS;
      
      // hash
      hashValue = (pattern_id40 >> HT_BITS) & 0xffffff;
      value_16bit = pattern_id40 & 0xffff;
      cnt_24bits = count[hashValue];
      
      /*
      index = 0;
      r = 0;
      if (cnt_24bits > 0){
	index = table[hashValue];
	l = ref_sal[index];
	for (unsigned a = 0; a < cnt_24bits; a++)
	  r += ref_sar[index + a];
      }
      
      if (index != 0){
	k = 36;
	while (r > candidates && k >= 0){
	  value_8bit = (pattern_id64_mask >> k) & 0x0f;
	  
	  if (value_8bit == 0x00){
	    k = k - 4;
	    continue;
	  }
	  
	  upper_bound = binary_search_upper<_val, _locr>(l, r, value_8bit, k/4, sa, part_len);
	  lower_bound = binary_search_lower<_val, _locr>(l, r, value_8bit, k/4, sa, part_len);
	  
	  if (upper_bound - lower_bound == 0){
	    // not found
	    break;
	  }
	  
	  r = upper_bound - lower_bound;
	  l = lower_bound;
	  
	  k = k - 4;
	}
	  
	r = min((unsigned)r, (unsigned)candidates);

	for (m = 0; m < r; m++){
	  loc = (_locr)sa[l+m] + part_len + (_locr)OFFSET;
	  score = align<_val,_locr>(&seq[j+OFFSET], loc);
	  if (score >= program_options::min_hit_score){
	    locs.push_back(make_pair(loc, j+OFFSET));
	  }
	}
      }
      */
      
      index = 0;
      //cout << cnt_24bits << endl;
      // check whether the high 24 bits (last 6 characters) of pattern_id64 in ref hash table
      // if cnt_24bits >= 0, then found on hash table
      if (cnt_24bits > 0){
	// found
	// check whether the last 16 bits of pattern_id64 exist in ref by binary search
	index = binary_search(table[hashValue], cnt_24bits, value_16bit, key16);
      }
      
      if (index != 0){
	// found
	// get SA range
	l = ref_sal[index];
	r = ref_sar[index];

        k = MAX_BS_CONST;
        while (r > candidates && k >= 0){
          value_8bit = (pattern_id64_mask >> k) & 0x0f;
	  
	  upper_bound = binary_search_upper<_val, _locr>(l, r, value_8bit, k/4, sa, part_len);
          lower_bound = binary_search_lower<_val, _locr>(l, r, value_8bit, k/4, sa, part_len);
	  
          if (upper_bound - lower_bound == 0){
            // not found
	    break;
          }
	  
          r = upper_bound - lower_bound;
          l = lower_bound;
	  
          k = k - 4;
        }
	
	r = min((unsigned)r, (unsigned)candidates);
	
	//r = (uint16_t)min((uint16_t)ref_sar[index], (uint16_t)CANDIDATES);
	
	for (m = 0; m < r; m++){
	  loc = (_locr)sa[l+m] + part_len + (_locr)OFFSET;
	  score = align<_val,_locr>(&seq[j+OFFSET], loc);
	  if (score >= program_options::min_hit_score){
	    locs.push_back(make_pair(loc, j+OFFSET));
	  }
	}	
      }

      j++;
    }
    
    if (locs.size() == 0) continue;
    
    std::sort(locs.begin(), locs.end(), cmp_locs<_locr, _locl>);
    
    band = program_options::read_padding<_val>(seq.length());
        
    (*out).push(hit<_locr,_locl> (i, locs[0].first, locs[0].second));
    
    for (n = 1; n < locs.size(); n++){
      ref_loc = locs[n].first;
      query_loc = locs[n].second;
      if ((ref_loc - query_loc) - (locs[n-1].first - locs[n-1].second) <= band){
	continue;
      }
      
      (*out).push(hit<_locr,_locl> (i, ref_loc, query_loc));
    }
    
  }

  delete []query_reduced_aa;
  delete out;
  return;
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void process_long(unsigned *sa_4bytes,
		  uint8_t *sa_1byte,
		  unsigned *table,
		  uint16_t *count,
		  unsigned *ref_sal_4bytes,
		  uint8_t *ref_sal_1byte,
		  unsigned *ref_sar,
		  uint16_t *key16,
		  uint64_t mis_pattern,
		  vector<unsigned> &split,
		  unsigned partition,
		  unsigned max_seq_len,
		  unsigned thread_id)
{
  typename Trace_pt_buffer<_locr,_locl>::Iterator* out = new typename Trace_pt_buffer<_locr,_locl>::Iterator (*Trace_pt_buffer<_locr,_locl>::instance, thread_id);

  const Sequence_set<_val> *query = query_seqs<_val>::data_;
  unsigned hashValue, m, n, seq_len, k1, k2, id, index;
  _locr i;
  unsigned r, cnt_24bits, band;
  uint16_t value_16bit;
  uint64_t pattern_id64, pattern_id64_mask, pattern_id40, tmp;
  uint8_t value_8bit;
  _locr loc, l, lower_bound, upper_bound;
  _locl j;
  int score, k;
  _locr ref_loc;
  _locl query_loc;
  
  unsigned candidates;
  if (program_options::aligner_mode == program_options::sensitive){
    candidates = CANDIDATES_SENSITIVE;
  }
  else {
    candidates = CANDIDATES_FAST;
  }
  
  vector<pair<_locr, _locl> > locs;

  uint64_t base[MAX_SEED_LEN];
  calculate_factor(base);

  uint8_t* query_reduced_aa = new uint8_t[max_seq_len];

  for (i = split[partition]; i < split[partition+1]; i++){
    const sequence<const _val> seq = (*query)[i];
    seq_len = seq.length();
    if (seq_len < MAX_SEED_LEN) continue;

    memset(query_reduced_aa, 0, sizeof(uint8_t) * max_seq_len);
    locs.clear();
    j = 0;

    while (j < seq_len - 1){
      if ((int)seq[j] == (int)Value_traits<Amino_acid>::MASK_CHAR){
        if ((int)seq[j+1] == (int)Value_traits<Amino_acid>::MASK_CHAR){
          query_reduced_aa[j/2] = 0;
        }
        else {
          k2 = Reduction<Amino_acid>::reduction((char)seq[j+1]) + 1;
	  query_reduced_aa[j/2] = k2 << 4;
        }
      }
      else {
        k1 = Reduction<Amino_acid>::reduction((char)seq[j]) + 1;
        if ((int)seq[j+1] == (int)Value_traits<Amino_acid>::MASK_CHAR){
          query_reduced_aa[j/2] = k1;
        }
        else {
          k2 = Reduction<Amino_acid>::reduction((char)seq[j+1]) + 1;
          query_reduced_aa[j/2] =  (k2 << 4) | k1;
        }
      }

      j += 2;
    }
    
    if ((seq_len & 0x01) != 0){
      if ((int)seq[j] == (int)Value_traits<Amino_acid>::MASK_CHAR){
        query_reduced_aa[j/2] = 0;
      }
      else {
        k1 = Reduction<Amino_acid>::reduction((char)seq[j]) + 1;
        query_reduced_aa[j/2] = k1;
      }
    }

    j = 0;
    seq_len = seq_len - MAX_SEED_LEN + 1;

    while (j < seq_len){
      id = 0;

      if ((j & 0x01) == 0){
        pattern_id64 = (uint64_t)*((uint64_t*)(query_reduced_aa+j/2));
      }
      else {
        pattern_id64 = ((uint64_t)*((uint64_t*)(query_reduced_aa+j/2)) & 0xfffffffffffffff0) >> 4;
        tmp = ((uint64_t)*((uint64_t*)(query_reduced_aa+(j+1)/2)) & 0x0f00000000000000) << 4;
        pattern_id64 = pattern_id64 | tmp;
      }

      for (k = MAX_SEED_LEN-1; k >= 0; k--){
        if ((pattern_id64 & base[k]) == 0){
          id = 1;
          j = j + k + 1;
          break;
        }
      }

      if (id == 1) continue;

      pattern_id64_mask = pattern_id64 & mis_pattern;
      pattern_id40 = pattern_id64_mask >> BS_BITS;
      
      // hash
      hashValue = (pattern_id40 >> HT_BITS) & 0xffffff;
      value_16bit = pattern_id40 & 0xffff;
      cnt_24bits = count[hashValue];

      index = 0;
      if (cnt_24bits > 0){
	index = binary_search(table[hashValue], cnt_24bits, value_16bit, key16);
      }

      if (index != 0){
	// found
        // get SA range
        l = ((_locr)ref_sal_4bytes[index] << 8) | ref_sal_1byte[index];
        r = ref_sar[index];
	
	k = MAX_BS_CONST;
        while (r > candidates && k >= 0){
          value_8bit = (pattern_id64_mask >> k) & 0x0f;
	  
          upper_bound = binary_search_upper_long<_val, _locr>(l, r, value_8bit, k/4, sa_4bytes, sa_1byte);
          lower_bound = binary_search_lower_long<_val, _locr>(l, r, value_8bit, k/4, sa_4bytes, sa_1byte);

          if (upper_bound - lower_bound == 0){
            // not found
            break;
          }
	  
          r = upper_bound - lower_bound;
          l = lower_bound;

          k = k - 4;
        }

        r = min((unsigned)r, (unsigned)candidates);
	
	for (m = 0; m < r; m++){
	  loc = (((_locr)sa_4bytes[l+m] << 8) | sa_1byte[l+m]) + OFFSET;
          score = align<_val,_locr>(&seq[j+OFFSET], loc);
          if (score >= program_options::min_hit_score){
            locs.push_back(make_pair(loc, j+OFFSET));
          }
        }
      }

      j++;
    }

    if (locs.size() == 0) continue;
    
    std::sort(locs.begin(), locs.end(), cmp_locs<_locr, _locl>);

    band = program_options::read_padding<_val>(seq.length());

    (*out).push(hit<_locr,_locl> (i, locs[0].first, locs[0].second));

    for (n = 1; n < locs.size(); n++){
      ref_loc = locs[n].first;
      query_loc = locs[n].second;
      if ((ref_loc - query_loc) - (locs[n-1].first - locs[n-1].second) <= band){
        continue;
      }

      (*out).push(hit<_locr,_locl> (i, ref_loc, query_loc));
    }

  }

  delete []query_reduced_aa;
  delete out;
  return;
}

template<typename _val, typename _locr, typename _locq, typename _locl>
struct Search_context
{
        Search_context(unsigned *sa,
		       unsigned *table,
		       uint16_t *count,
		       unsigned *ref_sal,
		       unsigned *ref_sar,
		       uint16_t *key16,
		       uint64_t mis_pattern,
		       vector<unsigned> &split,
		       unsigned max_seq_len,
		       _locr part_len):
                sa (sa),
		table (table),
		count (count),
		ref_sal (ref_sal),
		ref_sar (ref_sar),
		key16 (key16),
		mis_pattern (mis_pattern),
		split (split),
		max_seq_len (max_seq_len),
		part_len (part_len)
        { }
        void operator()(unsigned thread_id, unsigned partition) const
        {
	  process<_val, _locr, _locq, _locl>(sa, table, count, ref_sal, ref_sar, key16, mis_pattern, split, partition, max_seq_len, part_len, thread_id);
	}
        unsigned *sa;
        unsigned *table;
        uint16_t *count;
        unsigned *ref_sal;
        unsigned *ref_sar;
        uint16_t *key16;
        uint64_t mis_pattern;
        vector<unsigned> &split;
        unsigned max_seq_len;
        _locr part_len;
};

template<typename _val, typename _locr, typename _locq, typename _locl>
struct Search_context_long
{
        Search_context_long(unsigned *sa_4bytes,
			    uint8_t *sa_1byte,
			    unsigned *table,
			    uint16_t *count,
			    unsigned *ref_sal_4bytes,
			    uint8_t *ref_sal_1byte,
			    unsigned *ref_sar,
			    uint16_t *key16,
			    uint64_t mis_pattern,
			    vector<unsigned> &split,
			    unsigned max_seq_len):
                sa_4bytes (sa_4bytes),
		sa_1byte (sa_1byte),
		table (table),
		count (count),
		ref_sal_4bytes (ref_sal_4bytes),
		ref_sal_1byte (ref_sal_1byte),
		ref_sar (ref_sar),
		key16 (key16),
		mis_pattern (mis_pattern),
		split (split),
		max_seq_len (max_seq_len)
        { }
        void operator()(unsigned thread_id, unsigned partition) const
        {
	  process_long<_val, _locr, _locq, _locl>(sa_4bytes, sa_1byte, table, count, ref_sal_4bytes, ref_sal_1byte, ref_sar, key16, mis_pattern, split, partition, max_seq_len, thread_id);
	}
        unsigned *sa_4bytes;
        uint8_t *sa_1byte;
        unsigned *table;
        uint16_t *count;
        unsigned *ref_sal_4bytes;
        uint8_t *ref_sal_1byte;
        unsigned *ref_sar;
        uint16_t *key16;
        uint64_t mis_pattern;
        vector<unsigned> &split;
        unsigned max_seq_len;
};

template<typename _val, typename _locr, typename _locq, typename _locl>
void split_query(unsigned *sa,
		 unsigned *table,
                 uint16_t *count,
		 unsigned *ref_sal,
		 unsigned *ref_sar,
		 uint16_t *key16,
		 uint64_t mis_pattern,
		 cpu_timer &timer_mapping,
		 unsigned query_chunk,
		 vector<unsigned> &split,
		 unsigned max_seq_len,
		 _locr part_len)
{
  // search patterns in queries
  _locq n_query_seqs = query_seqs<_val>::data_->get_length();
  _locq cnt = n_query_seqs / Const::seedp;
  
  _locq j;
  
  for (unsigned i = 0; i < Const::seedp; i++){
    j = i * cnt;
    split.push_back(j);
  }
  split.push_back(n_query_seqs);
  
  task_timer timer ("Pattern searching", true);
  Search_context<_val, _locr, _locq, _locl> context (sa, table, count, ref_sal, ref_sar, key16, mis_pattern, split, max_seq_len, part_len);
  launch_scheduled_thread_pool(context, Const::seedp, program_options::threads());

  timer.finish();
  timer_mapping.stop();
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void split_query_long(unsigned *sa_4bytes,
		      uint8_t *sa_1byte,
		      unsigned *table,
		      uint16_t *count,
		      unsigned *ref_sal_4bytes,
		      uint8_t *ref_sal_1byte,
		      unsigned *ref_sar,
		      uint16_t *key16,
		      uint64_t mis_pattern,
		      cpu_timer &timer_mapping,
		      unsigned query_chunk,
		      vector<unsigned> &split,
		      unsigned max_seq_len)
{
  // search patterns in queries
  _locq n_query_seqs = query_seqs<_val>::data_->get_length();
  _locq cnt = n_query_seqs / Const::seedp;

  _locq j;

  for (unsigned i = 0; i < Const::seedp; i++){
    j = i * cnt;
    split.push_back(j);
  }
  split.push_back(n_query_seqs);
  
  task_timer timer ("Pattern searching", true);
  Search_context_long<_val, _locr, _locq, _locl> context (sa_4bytes, sa_1byte, table, count, ref_sal_4bytes, ref_sal_1byte, ref_sar, key16, mis_pattern, split, max_seq_len);
  launch_scheduled_thread_pool(context, Const::seedp, program_options::threads());
  
  timer.finish();
  timer_mapping.stop();
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void run_ref_chunk_fast(Database_file<_val> &db_file,
			cpu_timer &timer_mapping,
			cpu_timer &total_timer,
			unsigned query_chunk,
			pair<size_t,size_t> query_len_bounds,
			unsigned ref_chunk,
			DAA_output &master_out,
			vector<Temp_file> &tmp_file)
{
  verbose_stream << "Processing reference chunk " << ref_chunk << "." << endl;
  
  task_timer timer ("Loading reference sequences", true);
  ref_seqs<_val>::data_ = new Sequence_set<_val> (db_file);
  ref_ids::data_ = new String_set<char,0> (db_file);

  setup_search_params<_val>(query_len_bounds, ref_seqs<_val>::data_->letters());
  ref_map.init(ref_seqs<_val>::get().get_length());

  timer.go("Initializing temporary storage");
  timer_mapping.resume();
  Trace_pt_buffer<_locr,_locl>::instance = new Trace_pt_buffer<_locr,_locl> (query_seqs<_val>::data_->get_length()/query_contexts(), program_options::tmpdir, program_options::mem_buffered());
  
  timer.finish();
  timer_mapping.stop();
  
  unsigned max_seq_len = query_len_bounds.second;

  string file_prefix = program_options::database.substr(0,program_options::database.size()-5) + "_" + to_string(ref_chunk) + ".index";
  FILE *f = fopen(file_prefix.c_str(), "rb");

  unsigned *sa;
  unsigned *table;
  uint16_t *count;
  unsigned *ref_sal;
  unsigned *ref_sar;
  uint16_t *key16;
  unsigned sa_len, len, partN_len;
  _locr part_len = 0;
  
  for (unsigned i = 0; i < N_PART; i++){
    timer.go("Loading reference index");
    
    fread(&partN_len, sizeof(unsigned), 1, f);
    fread(&sa_len, sizeof(unsigned), 1, f);
    sa = (unsigned*)malloc(sa_len * sizeof(unsigned)); // suffix array
    fread(sa, sizeof(unsigned), sa_len, f);
    
    // perfect hash table: 24 bits per key
    table = (unsigned*)malloc(REF_TABLE_SIZE * sizeof(unsigned)); // record the first 7 characters of each pattern
    count = (uint16_t*)malloc(REF_TABLE_SIZE * sizeof(uint16_t)); // occurrences of first 7 characters of patterns
    fread(table, sizeof(unsigned), REF_TABLE_SIZE, f);
    fread(count, sizeof(uint16_t), REF_TABLE_SIZE, f);
    
    fread(&len, sizeof(unsigned), 1, f);
    //cout << "len " << len << endl;

    ref_sal = (unsigned*)malloc(len * sizeof(unsigned)); // SAL
    ref_sar = (unsigned*)malloc(len * sizeof(unsigned)); // SAR
    key16 = (uint16_t*)malloc(len * sizeof(uint16_t)); // last 4 characters of each pattern
    fread(ref_sal, sizeof(unsigned), len, f);
    fread(ref_sar, sizeof(unsigned), len, f);
    fread(key16, sizeof(uint16_t), len, f);
    
    timer.finish();
    
    uint64_t mis_pattern = MIS_PATTERN[0];
    vector<unsigned> split;
    split.clear();
    split_query<_val,_locr,_locq,_locl>(sa, table, count, ref_sal, ref_sar, key16, mis_pattern, timer_mapping, query_chunk, split, max_seq_len, part_len);
    
    part_len += (_locr)partN_len;

    timer.go("Deallocating reference index");
    free(sa);
    free(table);
    free(count);
    free(ref_sal);
    free(ref_sar);
    free(key16);
  }
  
  fclose(f);

  timer.go("Closing temporary storage");
  Trace_pt_buffer<_locr,_locl>::instance->close();
  exception_state.sync();
  
  timer_mapping.resume();
  Output_stream* out;
  if(ref_header.n_blocks > 1) {
    timer.go ("Opening temporary output file");
    tmp_file.push_back(Temp_file ());
    out = new Output_stream (tmp_file.back());
  } else
    out = &master_out.stream();
  
  timer.go("Computing alignments");
  align_queries<_val,_locr,_locl>(*Trace_pt_buffer<_locr,_locl>::instance, out);
  delete Trace_pt_buffer<_locr,_locl>::instance;
  
  if(ref_header.n_blocks > 1) {
    timer.go("Closing temporary output file");
    out->close();
    delete out;
  }
  timer_mapping.stop();
  
  timer.go("Deallocating reference sequences");
  delete ref_seqs<_val>::data_;
  delete ref_ids::data_;
  
  timer.finish();
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void run_ref_chunk_long(Database_file<_val> &db_file,
			cpu_timer &timer_mapping,
			cpu_timer &total_timer,
			unsigned query_chunk,
			pair<size_t,size_t> query_len_bounds,
			unsigned ref_chunk,
			DAA_output &master_out,
			vector<Temp_file> &tmp_file,
			unsigned count_patterns)
{
  verbose_stream << "Processing reference chunk " << ref_chunk << "." << endl;
  
  task_timer timer ("Loading reference sequences", true);
  ref_seqs<_val>::data_ = new Sequence_set<_val> (db_file);
  ref_ids::data_ = new String_set<char,0> (db_file);

  setup_search_params<_val>(query_len_bounds, ref_seqs<_val>::data_->letters());
  ref_map.init(ref_seqs<_val>::get().get_length());

  timer.go("Initializing temporary storage");
  timer_mapping.resume();
  Trace_pt_buffer<_locr,_locl>::instance = new Trace_pt_buffer<_locr,_locl> (query_seqs<_val>::data_->get_length()/query_contexts(), program_options::tmpdir, program_options::mem_buffered());

  timer.finish();
  timer_mapping.stop();

  unsigned max_seq_len = query_len_bounds.second;

  string file_prefix = program_options::database.substr(0,program_options::database.size()-5) + "_" + to_string(ref_chunk) + ".index";
  FILE *f = fopen(file_prefix.c_str(), "rb");
  
  unsigned *sa_4bytes;
  uint8_t *sa_1byte;
  unsigned *table;
  uint16_t *count;
  unsigned *ref_sal_4bytes;
  uint8_t *ref_sal_1byte;
  unsigned *ref_sar;
  uint16_t *key16;
  _locr sa_len, len;
  unsigned begin, i;
  
  if (count_patterns == 1) begin = 0;
  else begin = 1;
  
  for (i = begin; i < count_patterns; i++){
    timer.go("Loading reference index");
    cout << "------ Pattern " << i << " ------" << endl;

    fread(&sa_len, sizeof(_locr), 1, f);
    sa_4bytes = (unsigned*)malloc(sa_len * sizeof(unsigned));
    fread(sa_4bytes, sizeof(unsigned), sa_len, f);
    sa_1byte = (uint8_t*)malloc(sa_len * sizeof(uint8_t));
    fread(sa_1byte, sizeof(uint8_t), sa_len, f);
    
    // perfect hash table: 24 bits per key
    table = (unsigned*)malloc(REF_TABLE_SIZE * sizeof(unsigned)); // record the first 7 characters of each pattern
    count = (uint16_t*)malloc(REF_TABLE_SIZE * sizeof(uint16_t)); // occurrences of first 7 characters of patterns
    fread(table, sizeof(unsigned), REF_TABLE_SIZE, f);
    fread(count, sizeof(uint16_t), REF_TABLE_SIZE, f);
    
    fread(&len, sizeof(_locr), 1, f);
    ref_sal_4bytes = (unsigned*)malloc(len * sizeof(unsigned)); // SAL
    ref_sal_1byte = (uint8_t*)malloc(len * sizeof(uint8_t));
    ref_sar = (unsigned*)malloc(len * sizeof(unsigned)); // SAR
    key16 = (uint16_t*)malloc(len * sizeof(uint16_t)); // last 4 characters of each pattern
    fread(ref_sal_4bytes, sizeof(unsigned), len, f);
    fread(ref_sal_1byte, sizeof(uint8_t), len, f);
    fread(ref_sar, sizeof(unsigned), len, f);
    fread(key16, sizeof(uint16_t), len, f);
    
    timer.finish();
  
    uint64_t mis_pattern = MIS_PATTERN[i];
    vector<unsigned> split;
    split.clear();
    split_query_long<_val,_locr,_locq,_locl>(sa_4bytes, sa_1byte, table, count, ref_sal_4bytes, ref_sal_1byte, ref_sar, key16, mis_pattern, timer_mapping, query_chunk, split, max_seq_len);
    
    timer.go("Deallocating reference index");
    free(sa_4bytes);
    free(sa_1byte);
    free(table);
    free(count);
    free(ref_sal_4bytes);
    free(ref_sal_1byte);
    free(ref_sar);
    free(key16);
  }

  fclose(f);

  timer.go("Closing temporary storage");
  Trace_pt_buffer<_locr,_locl>::instance->close();
  exception_state.sync();

  timer_mapping.resume();
  Output_stream* out;
  if(ref_header.n_blocks > 1) {
    timer.go ("Opening temporary output file");
    tmp_file.push_back(Temp_file ());
    out = new Output_stream (tmp_file.back());
  } else
    out = &master_out.stream();
  
  timer.go("Computing alignments");
  align_queries<_val,_locr,_locl>(*Trace_pt_buffer<_locr,_locl>::instance, out);
  delete Trace_pt_buffer<_locr,_locl>::instance;
  
  if(ref_header.n_blocks > 1) {
    timer.go("Closing temporary output file");
    out->close();
    delete out;
  }
  timer_mapping.stop();
  
  timer.go("Deallocating reference sequences");
  delete ref_seqs<_val>::data_;
  delete ref_ids::data_;
  
  timer.finish();
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void run_ref_chunk_sensitive(Database_file<_val> &db_file,
			     cpu_timer &timer_mapping,
			     cpu_timer &total_timer,
			     unsigned query_chunk,
			     pair<size_t,size_t> query_len_bounds,
			     unsigned ref_chunk,
			     DAA_output &master_out,
			     vector<Temp_file> &tmp_file)
{
  verbose_stream << "Processing reference chunk " << ref_chunk << "." << endl;
  
  task_timer timer ("Loading reference sequences", true);
  ref_seqs<_val>::data_ = new Sequence_set<_val> (db_file);
  ref_ids::data_ = new String_set<char,0> (db_file);

  setup_search_params<_val>(query_len_bounds, ref_seqs<_val>::data_->letters());
  ref_map.init(ref_seqs<_val>::get().get_length());
  
  timer.go("Initializing temporary storage");
  timer_mapping.resume();
  Trace_pt_buffer<_locr,_locl>::instance = new Trace_pt_buffer<_locr,_locl> (query_seqs<_val>::data_->get_length()/query_contexts(), program_options::tmpdir, program_options::mem_buffered());
  timer.finish();
  timer_mapping.stop();
  
  unsigned max_seq_len = query_len_bounds.second;

  string file_prefix = program_options::database.substr(0,program_options::database.size()-5) +	"_" + to_string(ref_chunk) + ".index";
  FILE *f = fopen(file_prefix.c_str(), "rb");

  unsigned *sa;
  unsigned *table;
  uint16_t *count;
  unsigned *ref_sal;
  unsigned *ref_sar;
  uint16_t *key16;
  unsigned sa_len, len;

  for (unsigned i = 1; i < N_PATTERN; i++){
    timer.go("Loading reference index");
    cout << "------ Pattern " << i << " ------" << endl;

    fread(&sa_len, sizeof(unsigned), 1, f);
    sa = (unsigned*)malloc(sa_len * sizeof(unsigned)); // suffix array
    fread(sa, sizeof(unsigned), sa_len, f);
    
    // perfect hash table: 24 bits per key
    table = (unsigned*)malloc(REF_TABLE_SIZE * sizeof(unsigned)); // record the first 7 characters of each pattern
    count = (uint16_t*)malloc(REF_TABLE_SIZE * sizeof(uint16_t)); // occurrences of first 7 characters of patterns
    fread(table, sizeof(unsigned), REF_TABLE_SIZE, f);
    fread(count, sizeof(uint16_t), REF_TABLE_SIZE, f);
    
    fread(&len, sizeof(unsigned), 1, f);
    //cout << "len " << len << endl;
    ref_sal = (unsigned*)malloc(len * sizeof(unsigned)); // SAL
    ref_sar = (unsigned*)malloc(len * sizeof(unsigned)); // SAR
    key16 = (uint16_t*)malloc(len * sizeof(uint16_t)); // last 4 characters of each pattern
    fread(ref_sal, sizeof(unsigned), len, f);
    fread(ref_sar, sizeof(unsigned), len, f);
    fread(key16, sizeof(uint16_t), len, f);
    
    timer.finish();

    uint64_t mis_pattern = MIS_PATTERN[i];
    vector<unsigned> split;
    split.clear();
    split_query<_val,_locr,_locq,_locl>(sa, table, count, ref_sal, ref_sar, key16, mis_pattern, timer_mapping, query_chunk, split, max_seq_len, 0);
    
    timer.go("Deallocating reference index");
    free(sa);
    free(table);
    free(count);
    free(ref_sal);
    free(ref_sar);
    free(key16);
  }

  fclose(f);

  timer.go("Closing temporary storage");
  Trace_pt_buffer<_locr,_locl>::instance->close();
  exception_state.sync();

  timer_mapping.resume();
  Output_stream* out;
  if(ref_header.n_blocks > 1) {
    timer.go ("Opening temporary output file");
    tmp_file.push_back(Temp_file ());
    out = new Output_stream (tmp_file.back());
  } else
    out = &master_out.stream();
  
  timer.go("Computing alignments");
  align_queries<_val,_locr,_locl>(*Trace_pt_buffer<_locr,_locl>::instance, out);
  delete Trace_pt_buffer<_locr,_locl>::instance;

  if(ref_header.n_blocks > 1) {
    timer.go("Closing temporary output file");
    out->close();
    delete out;
  }
  timer_mapping.stop();
  
  timer.go("Deallocating reference sequences");
  delete ref_seqs<_val>::data_;
  delete ref_ids::data_;

  timer.finish();
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void run_query_chunk(Database_file<_val> &db_file,
		     cpu_timer &timer_mapping,
		     cpu_timer &total_timer,
		     unsigned query_chunk,
		     pair<size_t,size_t> query_len_bounds,
		     DAA_output &master_out)
{
  task_timer timer ("Allocating buffers", true);
  vector<Temp_file> tmp_file;
  timer.finish();
  unsigned count_patterns;

  if (program_options::aligner_mode == program_options::sensitive) {
    count_patterns = N_PATTERN;
    if (ref_header.long_addressing){
      for (current_ref_block=0;current_ref_block<ref_header.n_blocks;current_ref_block++){
        run_ref_chunk_long<_val,_locr,_locq,_locl>(db_file, timer_mapping, total_timer, query_chunk, query_len_bounds, current_ref_block, master_out, tmp_file, count_patterns);
      }
    }
    else {
      for (current_ref_block=0;current_ref_block<ref_header.n_blocks;current_ref_block++){
	run_ref_chunk_sensitive<_val,_locr,_locq,_locl>(db_file, timer_mapping, total_timer, query_chunk, query_len_bounds, current_ref_block, master_out, tmp_file);
      }
    }
  }
  else {
    if (ref_header.long_addressing){
      count_patterns = 1;
      for (current_ref_block=0;current_ref_block<ref_header.n_blocks;current_ref_block++){
        run_ref_chunk_long<_val,_locr,_locq,_locl>(db_file, timer_mapping, total_timer, query_chunk, query_len_bounds, current_ref_block, master_out, tmp_file, count_patterns);
      }
    }
    else {
      for (current_ref_block=0;current_ref_block<ref_header.n_blocks;current_ref_block++){
	run_ref_chunk_fast<_val,_locr,_locq,_locl>(db_file, timer_mapping, total_timer, query_chunk, query_len_bounds, current_ref_block, master_out, tmp_file);
      }
    }
  }

  timer.go("Deallocating buffers");
  timer_mapping.resume();
  
  if(ref_header.n_blocks > 1) {
    timer.go("Joining output blocks");
    join_blocks<_val>(ref_header.n_blocks, master_out, tmp_file);
  }
  
  timer.go("Deallocating queries");
  delete query_seqs<_val>::data_;
  delete query_ids::data_;
  delete query_source_seqs::data_;
  timer_mapping.stop();
}

template<typename _val, typename _locr>
void master_thread(Database_file<_val> &db_file, cpu_timer &timer_mapping, cpu_timer &total_timer)
{
  //shape_config::instance = shape_config (program_options::index_mode, _val ());
  
  task_timer timer ("Opening the input file", true);
  timer_mapping.resume();
  const Sequence_file_format<Nucleotide> *format_n (guess_format<Nucleotide>(program_options::query_file));
  const Sequence_file_format<Amino_acid> *format_a (guess_format<Amino_acid>(program_options::query_file));
  Input_stream query_file (program_options::query_file, true);
  current_query_chunk=0;
  
  timer.go("Opening the output file");
  DAA_output master_out;
  timer_mapping.stop();
  timer.finish();
    
  for(;;++current_query_chunk) {
    db_file.rewind();
    task_timer timer ("Loading query sequences", true);
    timer_mapping.resume();
    size_t n_query_seqs;
    if(input_sequence_type() == nucleotide)
      n_query_seqs = load_seqs<Nucleotide,_val,Double_strand>(query_file, *format_n, query_seqs<_val>::data_, query_ids::data_, query_source_seqs::data_, (size_t)(program_options::query_chunk_size * 1e9));
    else
      n_query_seqs = load_seqs<Amino_acid,_val,Single_strand>(query_file, *format_a, query_seqs<_val>::data_, query_ids::data_, query_source_seqs::data_, (size_t)(program_options::query_chunk_size * 1e9));
    if(n_query_seqs == 0)
      break;
    timer.finish();
    query_seqs<_val>::data_->print_stats();
    
    if(sequence_type() == amino_acid && program_options::seg == "yes") {
      timer.go("Running complexity filter");
      Complexity_filter<_val>::get().run(*query_seqs<_val>::data_);
    }
    
    verbose_stream << "Processing query chunk " << current_query_chunk << "." << endl;

    const pair<size_t,size_t> query_len_bounds = query_seqs<_val>::data_->len_bounds(10);

    timer_mapping.stop();
    timer.finish();
    const bool long_addressing_query = query_seqs<_val>::data_->raw_len() > (size_t)std::numeric_limits<uint32_t>::max();
    
    if (query_len_bounds.second <= (size_t)std::numeric_limits<uint8_t>::max()){
      if (long_addressing_query)
	run_query_chunk<_val,_locr,uint64_t,uint8_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
      else
	run_query_chunk<_val,_locr,uint32_t,uint8_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
    } else if (query_len_bounds.second <= (size_t)std::numeric_limits<uint16_t>::max()){
      if (long_addressing_query)
	run_query_chunk<_val,_locr,uint64_t,uint16_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
      else
	run_query_chunk<_val,_locr,uint32_t,uint16_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
    } else {
      if (long_addressing_query)
	run_query_chunk<_val,_locr,uint64_t,uint32_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
      else
	run_query_chunk<_val,_locr,uint32_t,uint32_t>(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, master_out);
    }
  }

  timer.go("Closing the output file");
  timer_mapping.resume();
  master_out.finish();
  timer_mapping.stop();
  
  timer.go("Closing the database file");
  db_file.close();
  
  timer.finish();
  verbose_stream << "Total time = " << boost::timer::format(total_timer.elapsed(), 1, "%ws\n");
  verbose_stream << "Mapping time = " << boost::timer::format(timer_mapping.elapsed(), 1, "%ws\n");
  statistics.print();
}

template<typename _val>
void master_thread()
{
  cpu_timer timer2, timer_mapping;
  timer_mapping.stop();
  
  if(!check_dir(program_options::tmpdir))
    throw std::runtime_error("Temporary directory " + program_options::tmpdir + " does not exist or is not a directory. Please use option -t to specify a different directory.");
  
  task_timer timer ("Opening the database", 1);
  Database_file<_val> db_file;
  timer.finish();
  program_options::set_options<_val>(ref_header.block_size);
  verbose_stream << "Reference = " << program_options::database << endl;
  verbose_stream << "Sequences = " << ref_header.sequences << endl;
  verbose_stream << "Letters = " << ref_header.letters << endl;
  verbose_stream << "Block size = " << (size_t)(ref_header.block_size * 1e9) << endl;
  
  if (ref_header.block_size > 4)
    master_thread<_val,uint64_t>(db_file, timer_mapping, timer2);
  else
    master_thread<_val,unsigned>(db_file, timer_mapping, timer2);
}

#endif /* MASTER_THREAD_H_ */
