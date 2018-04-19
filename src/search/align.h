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

#ifndef ALIGN_H_
#define ALIGN_H_

#include "../data/reference.h"
#include "../basic/statistics.h"
#include "../basic/score_matrix.h"
#include "../search/sse_dist.h"
#include "../search/align_ungapped.h"
#include "../dp/score_vector2.h"
#include "../dp/score_profile2.h"

	//#define likely(x) __builtin_expect((x), 1)
	//#define unlikely(x) __builtin_expect((x), 0)

template<typename _val, typename _locr>
int align(const _val *query,
	  _locr s)
{
  const _val* subject = ref_seqs<_val>::data_->data(s);
  
  if (match_block(query-8, subject-8) + match_block(query+8, subject+8) < program_options::min_identities) return -1;
  
  int score = xdrop_ungapped<_val,_locr>(query, subject);  
  return score;
}

/*
template<typename _val, typename _locr, typename _locl>
void align_ungap_group(const sequence<const _val> query,
			 vector<_locr> &refs,
		       vector<pair<_locr, _locl> > &locs,
		       int suffix_len,
		       int prefix_len,
		       int query_offset)
{
  typedef score_vector2<int8_t> sv;
  int i, j;
  vector<sequence<const _val> > ref_suffix_set;
  vector<sequence<const _val> > ref_prefix_set;
  sequence_stream2 dseq_ref;
  
  suffix_len = min((int)max(program_options::window, SEED_LEN2 - Const::seed_anchor) + Const::seed_anchor, suffix_len);
  prefix_len = min((int)max(program_options::window, (unsigned)Const::seed_anchor) - Const::seed_anchor, prefix_len);
  sequence<const _val> query_suffix (&query[query_offset], suffix_len, 0);
  sequence<const _val> query_prefix (&query[query_offset-prefix_len], prefix_len, 0);
  
  for (i = 0; i < (int)refs.size(); i++){
    sequence<const _val> ref_suffix (ref_seqs<_val>::data_->window_infix_dp_right(refs[i], suffix_len));
    sequence<const _val> ref_prefix (ref_seqs<_val>::data_->window_infix_dp_left(refs[i]-1, prefix_len));
    ref_suffix_set.push_back(ref_suffix);
    ref_prefix_set.push_back(ref_prefix);
  }
  
  sv xdrop (0-program_options::xdrop);
  sv min_score (program_options::min_hit_score);
  sv zero (0);
  sv score, st, best, tmp, stop;
  __m128i seq;
  int n_ref, n;

  n = 0;
  typename vector<sequence<const _val> >::const_iterator ref_suffix_begin (ref_suffix_set.begin());
  typename vector<sequence<const _val> >::const_iterator ref_prefix_begin (ref_prefix_set.begin());
  while (ref_suffix_begin < ref_suffix_set.end()){
    n_ref = min((int)SIXTEEN, (int)(ref_suffix_set.end() - ref_suffix_begin));
    typename vector<sequence<const _val> >::const_iterator ref_suffix_end (ref_suffix_begin + n_ref);
    typename vector<sequence<const _val> >::const_iterator ref_prefix_end (ref_prefix_begin + n_ref);
    
    dseq_ref.reset();
    
    st = zero;
    best = zero;
    stop = sv (0xff);
    
    // seed
    for (i = 0; i < (int)SEED_LEN2; i++){
      seq = dseq_ref.get_right<_val>(ref_suffix_begin, ref_suffix_end, i, int8_t());
      score = sv ((int)query_suffix[i], seq);
      best = best + score;
    }
    
    // suffix
    st = best;
    for (i = 0; i < suffix_len; i++){
      seq = dseq_ref.get_right<_val>(ref_suffix_begin, ref_suffix_end, i, int8_t());
      
      // score
      score = sv ((int)query_suffix[i], seq);
      score.mask_stop(stop);
      st = st + score;
      best.max(st);
      // 0 means stop (best score >= min_score or st - best < xdrop)
      tmp = min_score.cmpgt(best);
      stop = (st - best).cmpgt(xdrop).mask_stop(tmp);
      
      if (stop.cmpgt_count() == 0) break;
    }

    // prefix
    dseq_ref.reset();
    st = best;
    stop = min_score.cmpgt(best);
    for (i = 0; i < prefix_len; i++){
      seq = dseq_ref.get_left<_val>(ref_prefix_begin, ref_prefix_end, i, int8_t());
      
      // score
      score = sv ((int)query_prefix[prefix_len-1-i], seq);
      score.mask_stop(stop);
      st = st + score;
      best.max(st);
      tmp = min_score.cmpgt(best);
      stop = (st - best).cmpgt(xdrop).mask_stop(tmp);
      
      if (stop.cmpgt_count() == 0) break;
    }
    
    for (i = 0; i < n_ref; i++){
      if ((int)best[i] >= program_options::min_hit_score){
	locs.push_back(make_pair(refs[i+n], query_offset));
      }
    }
    
    ref_suffix_begin += SIXTEEN;
    ref_prefix_begin += SIXTEEN;
    n += SIXTEEN;
  }

}
*/

#endif /* ALIGN_H_ */
