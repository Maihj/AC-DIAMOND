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

#ifndef DP_MATRIX2_H_
#define DP_MATRIX2_H_

#include <vector>
#include <boost/thread/tss.hpp>
#include "score_vector2.h" //added

using std::vector;
using boost::thread_specific_ptr;

template<typename _score>
void array_clear2(score_vector2<_score> *v, unsigned n)
{
        score_vector2<_score> *end (v+n);
	while(v < end)
	        *(v++) = score_vector2<_score> ();
}

template<typename _score>
struct DP_matrix2
{

	typedef score_vector2<_score> sv;

	struct Column_iterator
	{
	        Column_iterator(int column, sv* hgap_front, sv* score_front, int row_pos, int row_end, unsigned delta):
			row_pos_ (row_pos),
			row_end_ (row_end),
			delta_ (delta),
			hgap_ptr_ (hgap_front),
			score_ptr_ (score_front),
			d_ (delta == 0 ? sv ((int) 0) : *score_front)
		{ }
		inline bool at_end() const
		{ return row_pos_ >= row_end_; }
		inline void operator++()
		{ ++row_pos_; ++hgap_ptr_; ++score_ptr_; }
		inline sv hgap() const
		{ return *(hgap_ptr_+delta_); }
		inline sv diag() const
		{ return d_; }
		inline void set_hgap(const sv& x)
		{ *hgap_ptr_ = x; }
		inline void set_score(const sv& x)
		{ d_ = *(score_ptr_+delta_); *score_ptr_ = x; }
	        int row_pos_, row_end_, delta_;
	        sv *hgap_ptr_, *score_ptr_, d_;
	};

        DP_matrix2(int columns, int rows, int band):
		rows_ (rows),
		band_ (band),
		scores_ (&s2),
		hgap_ (&h2)
	{
		scores_->resize(2*band+1);
		hgap_->resize(2*band+2);
		hgap_front_ = &hgap_->front();
		score_front_ = &scores_->front();
	}

	inline void clear2()
	{
	        array_clear2(hgap_front_, 2*band_+2);
		array_clear2(score_front_, 2*band_+1);
		//cout << "hh:" << hgap_front_[0][0] << endl;
	}
  
  /*
	inline void band_range(unsigned column, unsigned& begin, unsigned& end)
	{
		if(column >= rows_) {
			begin = 0;
			end = band_;
		} else if(column >= 0) {
			unsigned pj (column);
			unsigned top_delta (pj >= band_ ? 0 : band_ - pj);
			unsigned query_start (pj >= band_ ? pj - band_ : 0);
			unsigned query_end (std::min(pj+band_+1, rows_));
			begin = top_delta;
			end = begin + query_end - query_start;
		}
	}*/

	inline Column_iterator begin(int column)
	{
	  if(column == 0) {
	    int query_end (rows_ >= band_+1 ? band_+1 : rows_);
	    return Column_iterator (column, hgap_front_+band_, score_front_+band_, 0, query_end, 0);
	  }
	  /*else if(column > rows_) {
	    unsigned query_start (rows_ >= band_ ? rows_-band_ : 0);
	    return Column_iterator (column, hgap_front_, score_front_, query_start, rows_, 1);
	    }*/
	  /*else if(column > band_){
	    unsigned pj (column);
            unsigned top_delta (0);
            unsigned query_start (pj - band_);
            unsigned query_end (std::min(pj+band_+1, rows_));
            return Column_iterator (column, hgap_front_+top_delta, score_front_+top_delta, query_start, query_end, 0, 1);
	    }*/
	  else {
	    int pj (column);
	    int top_delta (pj > band_ ? 0 : band_ - pj);
	    int query_start (pj > band_ ? pj - band_ : 0);
	    int query_end (std::min(pj+band_+1, rows_));
	    return Column_iterator (column, hgap_front_+top_delta, score_front_+top_delta, query_start, query_end, 1); //top_dalta: there are top_dalta vectors in hgap and score not be used
	  }
	}

	inline void sub_all(sv *ptr, const sv *end, const sv& x)
	{
		while(ptr < end)
			*(ptr++) -= x;
	}

	inline sv get_min(const sv *ptr, const sv *end) const
	{
		sv x (*(ptr++));
		while(ptr < end)
			x = x.min(*(ptr++));
		return x;
	}

        vector<sv> score_buffer2() const
        {
	  return *scores_;
	}

private:

	static thread_specific_ptr<vector<sv> > scores_ptr;
	static thread_specific_ptr<vector<sv> > hgap_ptr;

	const int rows_, band_;
	sv *hgap_front_, *score_front_;
	//Tls<vector<sv> > scores_, hgap_;
	vector<sv>* scores_, *hgap_;
	vector<sv> s2,h2;

};

template<typename _score> thread_specific_ptr<vector<score_vector2<_score> > > DP_matrix2<_score>::scores_ptr;
template<typename _score> thread_specific_ptr<vector<score_vector2<_score> > > DP_matrix2<_score>::hgap_ptr;

#endif /* DP_MATRIX2_H_ */
