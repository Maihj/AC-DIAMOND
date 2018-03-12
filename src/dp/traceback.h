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

#ifndef TRACEBACK_H_
#define TRACEBACK_H_

#include "score_vector2.h"

#define BASE 100
#define MATCH 0
#define HGAP 1
#define VGAP 3

template<typename _score>
struct Traceback_matrix
{
        Traceback_matrix(score_vector2<_score> *scores, int band, int num):
                scores_ (scores),
		band_ (band+1),
		width_ (2*band+2),
		num_ (num)
	{ }
        int operator()(int col, int row) const
        {
	  int top_delta (col > band_ ? col - band_ : 0);
	  int pos = col*width_+row-top_delta;
	  return (int)scores_[pos][num_];
	}
  /*
        bool in_band(int col, int row) const
        {
	  return abs(col - row) < band_ && row > 0 && col > 0;
	}
  */
        /*void print(int col, int row) const
	{
		for(unsigned j=0;j<=row;++j) {
			for(unsigned i=0;i<=col;++i)
				printf("%4i", in_band(i, j) ? this->operator()(i, j) : 0);
			printf("\n");
		}
		}*/
private:
	const score_vector2<_score> *scores_;
	const int band_;
        const int width_;
        const int num_;
};

template<typename _val, typename _score>
void dp_traceback_right(local_match<_val> &segment,
			sequence<const _val> &ref_suffix,
			sequence<const _val> &query_suffix,
			score_vector2<_score> *relevant_scores,
			int num,
			int i,
			int j,
			int band,
			vector<char> &transcript_buf)
{
  //if (j < 0 || i < 0) return;

	Traceback_matrix<_score> dp (relevant_scores, band, num);

	segment.query_len_ = i + 1;
  	segment.subject_len_ = j + 1;
  	segment.query_begin_ = 0;
  	segment.subject_begin_ = 0;
	
	int gap_len, width, best_score;
	Edit_transcript transcript (transcript_buf);
	width = band * 2 + 2;
	best_score = 0;
	i++;
	j++;
	
	while (i > 0 && j > 0){
	  //cout << "dp(j,i): (" << j << ", " << i << ") " << (int)dp(j,i) << endl;
	  if (dp(j, i) == MATCH){
	    //printf("i=%i j=%i subject=%c query=%c \n",i,j,Value_traits<_val>::ALPHABET[mask_critical(ref_suffix[j-1])],Value_traits<_val>::ALPHABET[query_suffix[i-1]]);
	    best_score += (int)relevant_scores[j*width][num];
	    if(query_suffix[i-1] == mask_critical(ref_suffix[j-1]))
	      ++segment.identities_;
	    else
	      ++segment.mismatches_;
	    --i;
	    --j;
	    ++segment.len_;
	    transcript_buf.push_back(op_match);
	  }
	  else if (dp(j, i) == HGAP){
	    // hgap
	    gap_len = 1;
	    
	    best_score += (int)relevant_scores[j*width][num];
	    j--;
	    while (dp(j, i) == HGAP){
	      best_score += (int)relevant_scores[j*width][num];
	      j--;
	      gap_len++;
	    }
	    
	    ++segment.gap_openings_;
	    segment.len_ += gap_len;
	    transcript_buf.insert(transcript_buf.end(), gap_len, op_deletion);
	  }
	  else {
	    // vgap
	    gap_len = 1;
	    
	    i--;
	    while (dp(j, i) == VGAP){
	      i--;
	      gap_len++;
	    }
	    
	    ++segment.gap_openings_;
	    segment.len_ += gap_len;
	    transcript_buf.insert(transcript_buf.end(), gap_len, op_insertion);
	  }
	}
	
	segment.score_ = best_score;
	segment.transcript_right_ = transcript.set_end(transcript_buf);
}

template<typename _val, typename _score>
void dp_traceback_left(local_match<_val> &segment,
		       sequence<const _val> &ref_prefix,
		       sequence<const _val> &query_prefix,
		       score_vector2<_score> *relevant_scores,
		       int num,
		       int i,
		       int j,
		       int band,
		       vector<char> &transcript_buf)
{
  //if (j < 0 || i < 0) return;
	
        int slen = ref_prefix.length();
        int qlen = query_prefix.length();
	//int score = segment.score_ + best_score - score_matrix::get().letter_score(query_prefix[qlen-1], mask_critical(ref_prefix[slen-1]));
	//if (score < min_score) return;

	Traceback_matrix<_score> dp (relevant_scores, band, num);
	
	local_match<_val> segment_left;
	segment_left.query_len_ = i + 1;
	segment_left.subject_len_ = j + 1;
	segment_left.query_begin_ = 0;
	segment_left.subject_begin_ = 0;
		
	int gap_len, width, best_score;
	Edit_transcript transcript (transcript_buf);
	width = band * 2 + 2;
	best_score = 0;
        i++;
        j++;
	
	while (i > 0 && j > 0){
	  //cout << "dp(j,i): (" << j << ", " << i << ") " << (int)dp(j,i) << endl;
	  if (dp(j, i) == MATCH){
	    best_score += (int)relevant_scores[j*width][num];
	    if(query_prefix[qlen-i] == mask_critical(ref_prefix[slen-j]))
	      ++segment_left.identities_;
	    else
	      ++segment_left.mismatches_;
	    --i;
	    --j;
	    ++segment_left.len_;
	    transcript_buf.push_back(op_match);
	  }
	  else if (dp(j, i) == HGAP){
	    // hgap
	    gap_len = 1;
	    
	    best_score += (int)relevant_scores[j*width][num];
	    j--;
	    while (dp(j, i) == HGAP){
	      best_score += (int)relevant_scores[j*width][num];
	      j--;
	      gap_len++;
	    }
	    
	    ++segment_left.gap_openings_;
	    segment_left.len_ += gap_len;
	    transcript_buf.insert(transcript_buf.end(), gap_len, op_deletion);
	  }
	  else {
	    // vgap
	    gap_len = 1;
	    
	    i--;
	    while (dp(j, i) == VGAP){
	      i--;
	      gap_len++;
	    }
	    
	    ++segment_left.gap_openings_;
	    segment_left.len_ += gap_len;
	    transcript_buf.insert(transcript_buf.end(), gap_len, op_insertion);
	  }
	}
	
	segment_left.score_ = best_score;
       	segment_left.transcript_right_ = transcript.set_end(transcript_buf);
	
	segment -= segment_left;
	segment.query_begin_--;
	segment.subject_begin_--;
	segment.score_ -= score_matrix::get().letter_score(query_prefix[qlen-1], mask_critical(ref_prefix[slen-1]));
	if (query_prefix[qlen-1] == mask_critical(ref_prefix[slen-1]))
	  segment.identities_--;
	else
	  segment.mismatches_--;
	segment.len_--;
	segment.subject_len_--;
	segment.query_len_--;
	
}

/*template<typename _val, typename _score>
void dp_traceback(const _score &scores,
	     int gap_open,
	     int gap_extend,
	     int i,
	     int j,
	     int score)
{ return local_match<_val> (score); }
*/
#endif /* TRACEBACK_H_ */
