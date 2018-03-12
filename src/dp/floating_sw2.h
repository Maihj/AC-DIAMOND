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

#ifndef FLOATING_SW2_H_
#define FLOATING_SW2_H_

#include "../basic/match.h"
#include "dp_matrix2.h"
#include "score_vector2.h"
#include "score_profile2.h"
#include "../util/tinythread.h"
#include "traceback.h"

#define INF -64
	//#define GAP -12
	//#define GAP_EXTEND -1

using boost::thread_specific_ptr;

using std::max;
using std::min;

tthread::mutex output_mutex;

template<typename _val, typename _score>
void floating_sw_suffix(vector<sequence<const _val> > &ref_suffix_set,
			sequence<const _val> query_suffix,
			int xdrop,
			int gap_open,
			int gap_extend,
			int begin,
			int band,
			vector<local_match<_val> > &segments,
			score_vector2<_score> *traceback,
			score_vector2<_score> *Mv,
			score_vector2<_score> *Mh,
			vector<char> &transcript_buf)
{
  typedef score_vector2<_score> sv;
  
  int slen (ref_suffix_set[0].length());
  int qlen (query_suffix.length());
  //cout << "qlen slen " << qlen << " " << slen << endl;
  
  sv open_penalty (gap_open);
  sv extend_penalty (gap_extend);
  sv GAP (0-gap_open);
  sv GAP_EXTEND (0-gap_extend);
  sv zero (0);
  sv one (1);
  sv drop_score (0-xdrop);
  int i, j, n, pos, pos2, query_end, it;
  int width = band * 2 + 2;
  int width2 = width * 2;
  sv base (BASE);
  sv V_I, V_D, scores, direction;
  sv tmp, tmp2;
  int cnt, delta, best_row_index, best_column_index;
  sv mh, dh, mv, delta_col_current, delta_col_prev;
  n = 0;
  
  sequence_stream2 dseq_ref;
  score_profile2<_score> profile;
  
  typename vector<sequence<const _val> >::const_iterator ref_begin (ref_suffix_set.begin());
  
  while (ref_begin < ref_suffix_set.end()){
    const int n_ref (min((int)score_traits2<_score>::channels, (int)(ref_suffix_set.end() - ref_begin)));
    typename vector<sequence<const _val> >::const_iterator ref_end (ref_begin + n_ref);
    
    dseq_ref.reset();
    Mv[0] = zero;
    Mv[1] = GAP;
    Mv[width] = GAP;
    Mv[width+1] = GAP;
    Mv[width2] = base;
    Mv[width2+1] = base + GAP;

    for (i = 2; i < width; i++){
      Mv[i] = GAP_EXTEND;
      Mv[width+i] = GAP;
      Mv[width2+i] = Mv[width2+i-1] - one;
    }
    
    sv best (0), best_row (0), best_column1 (0), best_column2 (0), cell, row, row_max, column1, column2;
    sv mask (0xff);
    __m128i thisseq;

    traceback[0] = zero;

    //dseq_ref.reset2();
    // 0 < j <= band
    for (j = 0; j <= band; j++){
      it = 0;
      query_end = min(j + band + 1, qlen);
      
      delta_col_prev = traceback[j*width];
      delta_col_current = zero;
      
      mh = Mh[j*2];
      dh = Mh[j*2+1];
      
      row_max = zero;
      
      pos = (j+1)*width + 1;

      thisseq = dseq_ref.get_right<_val>(ref_begin, ref_end, j, _score());
      profile.template set<_val> (thisseq);
    
      while (it < query_end){
	mv = Mv[it];
	V_I = (mv - open_penalty).max(Mv[it + width] - extend_penalty);
	V_D = (mh - open_penalty).max(dh - extend_penalty);

	if ((it * j) == 0 && (it ^ j) != 0){
	  scores = sv (INF);
	  if (j == 0) V_I = sv (INF);
	  else V_D = sv (INF);
	}
	else {
	  scores = profile.get(query_suffix[it]);
	}
	
	scores.update_score(scores, V_I, V_D, direction);
	
	Mv[it + width] = V_I - mh;
	Mv[it] = scores - mh;
	dh = V_D - mv;
	mh = scores - mv;
	
	cell = mh - delta_col_prev;
	tmp2 = Mv[width2 + it];

	row = sv (j - it);
	
	delta_col_current.update_row(tmp2, cell, row, row_max);
	
	Mv[width2 + it] = cell;
	traceback[pos] = direction;
	
	it++;
	pos++;
      }
      
      Mv[width2 + it] = (Mv[width2 + it - 1] - one).max(zero);
      delta_col_current -= base;
      
      column1 = sv (j >> 7);
      column2 = sv (j & 0x7f);
      
      delta_col_current.mask_stop(mask);
      
      best.update_best(delta_col_current, row_max, column1, column2, best_row, best_column1, best_column2);
      
      traceback[(j+1)*width] = delta_col_current;
      
      mask = best.cmpgt(drop_score);
    }

    // j > band
    //dseq_ref.reset2();
    for (j = band + 1; j < slen; j++){
      it = j - band;
      delta = 1;
      query_end = min(j + band + 1, qlen);
      
      delta_col_prev = traceback[j*width];
      delta_col_current = zero;
      
      mh = Mh[j*2];
      dh = Mh[j*2+1];

      row_max = zero;
      //sv row_max2 (0);

      pos = (j+1)*width + 1;
      pos2 = 0;
      thisseq = dseq_ref.get_right<_val>(ref_begin, ref_end, j, _score());
      profile.template set<_val> (thisseq);
      
      while (it < query_end){
	mv = Mv[pos2 + delta];
        V_I = (mv - open_penalty).max(Mv[pos2 + width + delta] - extend_penalty);
	V_D = (mh - open_penalty).max(dh - extend_penalty);

	scores = profile.get(query_suffix[it]);
	scores.update_score(scores, V_I, V_D, direction);
	
	Mv[pos2 + width] = V_I - mh;
        Mv[pos2] = scores - mh;
        dh = V_D - mv;
	mh = scores - mv;
	
	cell = mh - delta_col_prev;
	tmp2 = Mv[width2 + pos2 + delta];
	
	row = sv (j - it);
        //sv row2 (it & 0x7f);

        delta_col_current.update_row(tmp2, cell, row, row_max);
	
	Mv[width2 + pos2] = cell;
	traceback[pos] = direction;
	
        it++;
	pos++;
        pos2++;
      }

      Mv[width2 + pos2] = (Mv[width2 + pos2 - 1] - one).max(zero);

      delta_col_current -= base;
      
      column1 = sv (j >> 7);
      column2 = sv (j & 0x7f);
      
      delta_col_current.mask_stop(mask);
      
      best.update_best(delta_col_current, row_max, column1, column2, best_row, best_column1, best_column2);
            
      traceback[(j+1)*width] = delta_col_current;
      
      mask = best.cmpgt(drop_score);
      cnt = mask.cmpgt_count();
      if (cnt == 0) break;
    }
    
    /*__m128i best_score1 = _mm_set1_epi16(-gap_open);
    __m128i best_score2 = best_score1;
    __m128i st1 = _mm_setzero_si128();
    __m128i st2 = _mm_setzero_si128();
    int max_col = -1;
    for (i = 0; i < n_ref; i++){
      max_col = max(max_col, (int)best_column1[i] << 7 | best_column2[i]);
    }
    
    for (i = 1; i <= max_col + 1; i++){
      cell = traceback[i*width];
      st1 = _mm_add_epi16(st1, _mm_srai_epi16(_mm_unpacklo_epi8(zero.data_, cell.data_), 8));
      st2 = _mm_add_epi16(st2, _mm_srai_epi16(_mm_unpackhi_epi8(zero.data_, cell.data_), 8));
      best_score1 = _mm_max_epi16(best_score1, st1);
      best_score2 = _mm_max_epi16(best_score2, st2);
    }
    */
    for (i = 0; i < n_ref; i++){
      best_column_index = (int)best_column1[i] << 7 | best_column2[i];
      best_row_index = best_column_index - (int)best_row[i];
      //int best_score = i < 8 ? *(((int16_t*)&best_score1)+i) : *(((int16_t*)&best_score2)+(i-8));
      dp_traceback_right<_val, _score>(segments[begin], ref_suffix_set[i+n], query_suffix, traceback, i, best_row_index, best_column_index, band, transcript_buf);
      begin++;
    }
        
    ref_begin += score_traits2<_score>::channels;
    n += score_traits2<_score>::channels;
  }
  
}

template<typename _val, typename _score>
void floating_sw_prefix(vector<sequence<const _val> > &ref_prefix_set,
			sequence<const _val> query_prefix,
			int xdrop,
			int gap_open,
			int gap_extend,
			int begin,
			int band,
			vector<local_match<_val> > &segments,
			score_vector2<_score> *traceback,
			score_vector2<_score> *Mv,
			score_vector2<_score> *Mh,
			vector<char> &transcript_buf)
{
  typedef score_vector2<_score> sv;
  
  int slen (ref_prefix_set[0].length());
  int qlen (query_prefix.length());
  
  sv open_penalty (gap_open);
  sv extend_penalty (gap_extend);
  sv GAP (0-gap_open);
  sv GAP_EXTEND (0-gap_extend);
  sv zero (0);
  sv one (1);
  sv drop_score (0-xdrop);
  int i, j, n, pos, pos2, query_end, it;
  int width = band * 2 + 2;
  int width2 = width * 2;
  sv base (BASE);
  sv V_I, V_D, scores, direction;
  sv tmp, tmp2;
  int cnt, delta, best_row_index, best_column_index;
  sv mh, dh, mv, delta_col_current, delta_col_prev;
  n = 0;
  sequence_stream2 dseq_ref;
  score_profile2<_score> profile;
  
  typename vector<sequence<const _val> >::const_iterator ref_begin (ref_prefix_set.begin());
  
  while (ref_begin < ref_prefix_set.end()){
    const int n_ref (min((int)score_traits2<_score>::channels, (int)(ref_prefix_set.end() - ref_begin)));
    typename vector<sequence<const _val> >::const_iterator ref_end (ref_begin + n_ref);
    
    dseq_ref.reset();
    Mv[0] = sv (0);
    Mv[1] = GAP;
    Mv[width] = GAP;
    Mv[width+1] = GAP;
    Mv[width2] = base;
    Mv[width2+1] = base + GAP;

    for (i = 2; i < width; i++){
      Mv[i] = GAP_EXTEND;
      Mv[width+i] = GAP;
      Mv[width2+i] = Mv[width2+i-1] - one;
    }
    
    sv best (0), best_row (0), best_column1 (0), best_column2 (0), cell, row_max, row, column1, column2;
    sv mask (0xff);
    __m128i thisseq;

    traceback[0] = zero;

    // 0 <= j <= band
    for (j = 0; j <= band; j++){
      it = 0;
      query_end = min(j + band + 1, qlen);
      
      delta_col_prev = traceback[j*width];
      delta_col_current = zero;
      
      mh = Mh[j*2];
      dh = Mh[j*2+1];
      
      row_max = zero;
            
      pos = (j+1)*width + 1;
      
      thisseq = dseq_ref.get_left<_val>(ref_begin, ref_end, j, _score());
      profile.template set<_val> (thisseq);
      
      while (it < query_end){
	mv = Mv[it];
	V_I = (mv - open_penalty).max(Mv[it + width] - extend_penalty);
        V_D = (mh - open_penalty).max(dh - extend_penalty);
	
	if ((it * j) == 0 && (it ^ j) != 0){
	  scores = sv (INF);
	  if (j == 0) V_I = sv (INF);
	  else V_D = sv (INF);
	}
	else {
	  scores = profile.get(query_prefix[qlen-it-1]);
	}
	
	scores.update_score(scores, V_I, V_D, direction);
	
	Mv[it + width] = V_I - mh;
        Mv[it] = scores - mh;
        dh = V_D - mv;
        mh = scores - mv;

	cell = mh - delta_col_prev;
        tmp2 = Mv[width2 + it];
	
	row = sv (j - it);
	
	delta_col_current.update_row(tmp2, cell, row, row_max);

	Mv[width2 + it] = cell;
        traceback[pos] = direction;
	
	it++;
	pos++;
      }

      Mv[width2 + it] = (Mv[width2 + it - 1] - one).max(zero);
      
      delta_col_current -= base;
      
      column1 = sv (j >> 7);
      column2 = sv (j & 0x7f);
      
      delta_col_current.mask_stop(mask);
      
      best.update_best(delta_col_current, row_max, column1, column2, best_row, best_column1, best_column2);

      traceback[(j+1)*width] = delta_col_current;
      
      mask = best.cmpgt(drop_score);
    }

    //dseq_ref.reset2();
    for (j = band + 1; j < slen; j++){
      it = j - band;
      delta = 1;
      query_end = min(j + band + 1, qlen);

      delta_col_prev = traceback[j*width];
      delta_col_current = zero;
      
      mh = Mh[j*2];
      dh = Mh[j*2+1];

      row_max = zero;

      pos = (j+1)*width + 1;
      pos2 = 0;
      
      thisseq = dseq_ref.get_left<_val>(ref_begin, ref_end, j, _score());
      profile.template set<_val> (thisseq);
      
      while (it < query_end){
        mv = Mv[pos2 + delta];
        V_I = (mv - open_penalty).max(Mv[pos2 + width + delta] - extend_penalty);
        V_D = (mh - open_penalty).max(dh - extend_penalty);

        scores = profile.get(query_prefix[qlen-it-1]);
        scores.update_score(scores, V_I, V_D, direction);
	
        Mv[pos2 + width] = V_I - mh;
        Mv[pos2] = scores - mh;
        dh = V_D - mv;
        mh = scores - mv;
	
	cell = mh - delta_col_prev;
        tmp2 = Mv[width2 + pos2 + delta];
	
	row = sv (j - it);
        
        delta_col_current.update_row(tmp2, cell, row, row_max);

	Mv[width2 + pos2] = cell;
        traceback[pos] = direction;

        it++;
        pos++;
        pos2++;
      }

      Mv[width2 + pos2] = (Mv[width2 + pos2 - 1] - one).max(zero);
      
      delta_col_current -= base;
      
      column1 = sv (j >> 7);
      column2 = sv (j & 0x7f);
      
      delta_col_current.mask_stop(mask);
      
      best.update_best(delta_col_current, row_max, column1, column2, best_row, best_column1, best_column2);
      
      traceback[(j+1)*width] = delta_col_current;
      
      mask = best.cmpgt(drop_score);
      cnt = mask.cmpgt_count();
      if (cnt == 0) break;
    }
    /*
    __m128i best_score1 = _mm_set1_epi16(-gap_open);
    __m128i best_score2 = best_score1;
    __m128i st1 = _mm_setzero_si128();
    __m128i st2 = _mm_setzero_si128();
    int max_col = -1;
    for (i = 0; i < n_ref; i++){
      max_col = max(max_col, (int)best_column1[i] << 7 | best_column2[i]);
    }
    
    for (i = 1; i <= max_col + 1; i++){
      cell = traceback[i*width];
      st1 = _mm_add_epi16(st1, _mm_srai_epi16(_mm_unpacklo_epi8(zero.data_, cell.data_), 8));
      st2 = _mm_add_epi16(st2, _mm_srai_epi16(_mm_unpackhi_epi8(zero.data_, cell.data_), 8));
      best_score1 = _mm_max_epi16(best_score1, st1);
      best_score2 = _mm_max_epi16(best_score2, st2);
    }
    */

    for (i = 0; i < n_ref; i++){
      best_column_index = (int)best_column1[i] << 7 | best_column2[i];
      best_row_index = best_column_index - (int)best_row[i];
      dp_traceback_left<_val, _score>(segments[begin], ref_prefix_set[i+n], query_prefix, traceback, i, best_row_index, best_column_index, band, transcript_buf);      
      begin++;
    }

    ref_begin += score_traits2<_score>::channels;
    n += score_traits2<_score>::channels;
  }
}

template<typename _val, typename _score>
void floating_sw2(vector<sequence<const _val> > &ref_suffix_set,
		  vector<sequence<const _val> > &ref_prefix_set,
		  sequence<const _val> query_suffix,
		  sequence<const _val> query_prefix,
		  int xdrop,
		  int gap_open,
		  int gap_extend,
		  int begin,
		  int band,
		  vector<local_match<_val> > &segments,
		  score_vector2<_score> *traceback,
		  score_vector2<_score> *Mv,
		  score_vector2<_score> *Mh,
		  vector<char> &transcript_buf,
		  const _score&)
{
  // suffixes
  //output_mutex.lock();
  floating_sw_suffix<_val, _score>(ref_suffix_set, query_suffix, xdrop, gap_open, gap_extend, begin, band, segments, traceback, Mv, Mh, transcript_buf);
  //output_mutex.unlock();
  // prefixes
  floating_sw_prefix<_val, _score>(ref_prefix_set, query_prefix, xdrop, gap_open, gap_extend, begin, band, segments, traceback, Mv, Mh, transcript_buf);
}

#endif /* FLOATING_SW2_H_ */
