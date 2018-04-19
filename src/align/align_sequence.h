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

/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef ALIGN_SEQUENCE_H_
#define ALIGN_SEQUENCE_H_

#include <vector>
#include "../dp/floating_sw2.h"

using std::vector;
using std::max;

template<typename _val>
bool cmp_seed_offset(const local_match<_val> &lhs, const local_match<_val> &rhs){
  return lhs.query_anchor_ < rhs.query_anchor_;
}

typedef score_vector2<int8_t> sv;

template<typename _val>
void align_group(vector<sequence<const _val> > &ref_suffix_set,
		 vector<sequence<const _val> > &ref_prefix_set,
		 const sequence<const _val> query,
		 int begin,
		 int end,
		 int band,
		 vector<local_match<_val> > &segments,
		 sv *traceback,
		 sv *Mv,
		 sv *Mh,
		 vector<char> &transcript_buf)
{
        int query_suffix_len, query_prefix_len, ref_suffix_len, ref_prefix_len, query_offset, i;
	const int query_len = query.length();
	query_offset = segments[begin].query_anchor_;
	query_suffix_len = query_len - query_offset;
	query_prefix_len = query_offset + 1;
	ref_suffix_len = query_suffix_len + band;
	ref_prefix_len = query_prefix_len + band;
	
	sequence<const _val> query_suffix (&query[query_offset], query_suffix_len, 0);
        sequence<const _val> query_prefix (&query[0], query_prefix_len, 0);
	ref_suffix_set.clear();
        ref_prefix_set.clear();
	
	for (i = begin; i <= end; i++){
	  sequence<const _val> ref_suffix (ref_seqs<_val>::data_->window_infix_dp_right(segments[i].subject_position_, ref_suffix_len));
	  sequence<const _val> ref_prefix (ref_seqs<_val>::data_->window_infix_dp_left(segments[i].subject_position_, ref_prefix_len));
	  ref_suffix_set.push_back(ref_suffix);
	  ref_prefix_set.push_back(ref_prefix);
	}
	
	floating_sw2(ref_suffix_set,
		     ref_prefix_set,
		     query_suffix,
		     query_prefix,
		     score_matrix::get().rawscore(program_options::gapped_xdrop),
		     program_options::gap_open + program_options::gap_extend,
		     program_options::gap_extend,
		     begin,
		     band,
		     segments,
		     traceback,
		     Mv,
		     Mh,
		     transcript_buf,
		     int8_t());
}

template<typename _val, typename _locr, typename _locl>
void align_sequence(vector<Segment<_val> > &matches,
		    Statistics &stat,
		    vector<local_match<_val> > &local,
		    unsigned *padding,
		    size_t db_letters,
		    unsigned dna_len,
		    typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &begin,
		    typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &end,
		    vector<char> &transcript_buf)
{
        const unsigned q_num (begin->query_);
	const sequence<const _val> query (query_seqs<_val>::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const int query_len = query.length();
	padding[frame] = program_options::read_padding<_val>(query_len);
	const Sequence_set<_val> *ref = ref_seqs<_val>::data_;

	int GAP = -program_options::gap_open - program_options::gap_extend;
	unsigned band, max_slen, width;
	int k, local_temp_size, start;
	band = padding[frame];
	width = band * 2 + 2;

	vector<local_match<_val> > local_temp;
	local_temp.clear();
	
	if (program_options::aligner_mode == program_options::sensitive){
	  std::sort(begin, end, hit<_locr,_locl>::cmp_normalized_subject);
	  
	  for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i) {
	    if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) {
	      stat.inc(Statistics::DUPLICATES);
	      continue;
	    }
	  
	    local_temp.push_back(local_match<_val> (i->seed_offset_, ref->data(i->subject_), i->subject_));
	  }
	}
	else {
	  for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i) {
	    local_temp.push_back(local_match<_val> (i->seed_offset_, ref->data(i->subject_), i->subject_));
          }
	}

	// sort for grouping
	std::sort(local_temp.begin(), local_temp.end(), cmp_seed_offset<_val>);
	
	// bottom local_match, which is useless
	local_temp.push_back(local_match<_val> (-1, ref->data(begin->subject_)));
	
	vector<sequence<const _val> > ref_suffix_set;
	vector<sequence<const _val> > ref_prefix_set;
	
	local_temp_size = local_temp.size() - 2;
	
	max_slen = max(query_len-local_temp[0].query_anchor_+band, local_temp[local_temp_size].query_anchor_+band+1);                                                                   
        max_slen += 2;

	sv *traceback = (sv*) malloc(max_slen * width * sizeof(sv));
        sv *Mv = (sv*)malloc(width * 3 * sizeof(sv));
        sv *Mh = (sv*)malloc(max_slen * 2 * sizeof(sv));
	
	unsigned a;
	// Mh
	memset(Mh, GAP, max_slen * 2 * sizeof(sv));
        Mh[0] = sv (0);   // mh
	// Dh
	for (a = 2; a <= band; a++){
          Mh[a*2] = sv (-1);
        }

        // traceback
	memset(traceback, BASE, max_slen * width * sizeof(sv));
	
	start = 0;
	for (k = 0; k <= local_temp_size; k++){
	  if (local_temp[k].query_anchor_ == local_temp[k+1].query_anchor_){
	    continue;
	  }
	  
	  // a group of hits range from local_temp[start] to local_temp[k]
	  align_group<_val>(ref_suffix_set, ref_prefix_set, query, start, k, band, local_temp, traceback, Mv, Mh, transcript_buf);
	  //count++;
	  start = k + 1;
	}
	//output_mutex.lock();
	//cout << "ha: " << count << endl;
	//output_mutex.unlock();

	for (k = 0; k <= local_temp_size; k++){
	  const int score = local_temp[k].score_;
          
	  local.push_back(local_temp[k]);
	  std::pair<size_t, size_t> l = ref_seqs<_val>::data_->local_position(local.back().subject_position_);
	  matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
	  anchored_transform(local.back(), l.second, local.back().query_anchor_);
	  stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);
	  
	  //local.back().print(query, ref_seqs<_val>::get()[l.first], transcript_buf);
	  //cout << local.back().subject_begin_ << " " << local.back().subject_len_ << " " << local.back().query_begin_ << " " << local.back().query_len_ << endl;
	  
	  to_source_space(local.back(), frame, dna_len);
	  stat.inc(Statistics::SCORE_TOTAL, score);
	  stat.inc(Statistics::OUT_HITS);
	}

	
	free(traceback);
	free(Mv);
	free(Mh);
	
	//output_mutex.unlock();
}

#endif /* ALIGN_SEQUENCE_H_ */
