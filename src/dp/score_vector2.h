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

#ifndef SCORE_VECTOR2_H_
#define SCORE_VECTOR2_H_

#ifdef __SSSE3__
#include <tmmintrin.h>
#include <smmintrin.h>
#include <emmintrin.h>
#endif

#include "../basic/score_matrix.h"

template<typename _score>
struct score_traits2
{
	static const unsigned channels = 1;
        enum { zero = 0, inf = -127, byte_size = 4 };
	typedef bool Mask;
};

template<>
struct score_traits2<int8_t>
{
        enum { channels = 16, zero = 0x00, inf = 0x81, byte_size = 1 };
	typedef int8_t Mask;
};

template<typename _score>
struct score_vector2
{ };

template<>
struct score_vector2<int8_t>
{

        score_vector2()
	{
	  data_ = _mm_set1_epi8(score_traits2<int8_t>::inf);
	}

	explicit score_vector2(int x):
		data_ (_mm_set1_epi8(x))
	{ }

	explicit score_vector2(__m128i data):
		data_ (data)
	{ }

	explicit score_vector2(int a, const __m128i &seq)
	{
	        if(program_options::have_ssse3) {
#ifdef __SSSE3__
		  set_ssse3(a, seq);
#else
		  set_generic(a, seq);
#endif
		}
		else
		  set_generic(a, seq);
	}

	void set_ssse3(int a, const __m128i &seq)
	{
#ifdef __SSSE3__
	        const __m128i *row = reinterpret_cast<const __m128i*>(&score_matrix::get().matrix8()[a << 5]);
		__m128i high_mask = _mm_slli_epi16(_mm_and_si128(seq, _mm_set1_epi8(0x10)), 3);
                __m128i seq_low = _mm_or_si128(seq, high_mask);
                __m128i seq_high = _mm_or_si128(seq, _mm_xor_si128(high_mask, _mm_set1_epi8(0x80)));
		__m128i r1 = _mm_load_si128(row);
		__m128i r2 = _mm_load_si128(row+1);
		__m128i s1 = _mm_shuffle_epi8(r1, seq_low);
		__m128i s2 = _mm_shuffle_epi8(r2, seq_high);
		data_ = _mm_or_si128(s1, s2);
#endif
	}

	void set_generic(unsigned a, const __m128i &seq)
        {
	        const int8_t* row (&score_matrix::get().matrix8()[a<<5]);
		const int8_t* seq_ptr (reinterpret_cast<const int8_t*>(&seq));
		int8_t* dest (reinterpret_cast<int8_t*>(&data_));
		for(unsigned i=0;i<16;i++)
		  *(dest++) = row[*(seq_ptr++)];
	}
	
	score_vector2(const int8_t* s):
		data_ (_mm_loadu_si128(reinterpret_cast<const __m128i*>(s)))
	{ }

	score_vector2 operator+(const score_vector2 &rhs) const
	{
		return score_vector2 (_mm_add_epi8(data_, rhs.data_));
	}

	score_vector2 operator-(const score_vector2 &rhs) const
	{
		return score_vector2 (_mm_sub_epi8(data_, rhs.data_));
	}

	score_vector2& operator-=(const score_vector2 &rhs)
	{
		data_ = _mm_sub_epi8(data_, rhs.data_);
		return *this;
	}

	//void unbias(const score_vector2 &bias)
	//{ this->operator -=(bias); }

	void p()
	{
	  //return *(((int16_t*)&data_)+i);
	        print2(data_);
	}

	int operator [](unsigned i) const
        {
	  /*int value;
	  int16_t *pt = (int16_t*)&data_;
	  for(unsigned j=0; j<i; j++)
	    value = (int)*(pt++);
	  return value;*/
	  return *(((int8_t*)&data_)+i);
	}

	void set(unsigned i, int8_t v)
	{
	        *(((int8_t*)&data_)+i) = v;
	}

	score_vector2& max(const score_vector2 &rhs)
	{
		data_ = _mm_max_epi8(data_, rhs.data_);
		return *this;
	}

	score_vector2& min(const score_vector2 &rhs)
	{
		data_ = _mm_min_epi8(data_, rhs.data_);
		return *this;
	}

	score_vector2& mask_highbit(const score_vector2 &rhs)
	{
	  __m128i tmp = _mm_set1_epi8(0x7f);
	  data_ = _mm_and_si128(rhs.data_, tmp);
	  return *this;
	}

	score_vector2& mask_stop(const score_vector2 &rhs)
	{
	  data_ = _mm_and_si128(data_, rhs.data_);
	  return *this;
	}

	score_vector2& update_score(score_vector2 &scores, score_vector2 &V_I, score_vector2 &V_D, score_vector2 &direction)
	{
	  __m128i select_V_I = _mm_cmpgt_epi8(V_I.data_, scores.data_);
	  direction.data_ = _mm_and_si128(select_V_I, _mm_set1_epi8(0x01));
	  __m128i temp = _mm_max_epi8(scores.data_, V_I.data_);
          __m128i select_V_D = _mm_cmpgt_epi8(V_D.data_, temp);
          direction.data_ = _mm_or_si128(direction.data_, _mm_and_si128(select_V_D, _mm_set1_epi8(0x03)));
	  
          data_ = _mm_max_epi8(temp, V_D.data_);
          return *this;
	}
	
	score_vector2& update_row(score_vector2 &rhs, score_vector2 &cell, const score_vector2 &row, score_vector2 &row_best)
	{
	  // select those rhs.data_ = 0 and do not want them
	  __m128i zero = _mm_setzero_si128();
	  __m128i temp = _mm_cmpgt_epi8(rhs.data_, zero);
	  __m128i select_cell = _mm_and_si128(cell.data_, temp);
	  cell.data_ = _mm_max_epi8(_mm_add_epi8(rhs.data_, select_cell), zero);
	  
	  // choose those relevant_scores on row.data_
	  temp = _mm_cmpgt_epi8(cell.data_, data_);
	  __m128i select_row = _mm_and_si128(row.data_, temp);
	  temp = _mm_xor_si128(temp, _mm_set1_epi8(0xff));
	  __m128i origin_row = _mm_and_si128(row_best.data_, temp);
	  row_best.data_ = _mm_or_si128(select_row, origin_row);
	  
	  data_ = _mm_max_epi8(data_, cell.data_);
	  return *this;
	}
	
	score_vector2& update_best(const score_vector2 &delta, const score_vector2 &row, const score_vector2 &column1, const score_vector2 &column2, score_vector2 &best_row, score_vector2 &best_column1, score_vector2 &best_column2)
	{
	  __m128i zero = _mm_setzero_si128();
	  __m128i new_best = _mm_add_epi8(data_, delta.data_);
	  __m128i temp = _mm_cmpgt_epi8(new_best, zero);
	  __m128i select_row = _mm_and_si128(row.data_, temp);
	  __m128i select_column1 = _mm_and_si128(column1.data_, temp);
	  __m128i select_column2 = _mm_and_si128(column2.data_, temp);
	  temp = _mm_xor_si128(temp, _mm_set1_epi8(0xff)); // select_best (0)
	  __m128i origin_row = _mm_and_si128(best_row.data_, temp);
	  __m128i origin_column1 = _mm_and_si128(best_column1.data_, temp);
	  __m128i origin_column2 = _mm_and_si128(best_column2.data_, temp);
	  best_row.data_ = _mm_or_si128(select_row, origin_row);
	  best_column1.data_ = _mm_or_si128(select_column1, origin_column1);
	  best_column2.data_ = _mm_or_si128(select_column2, origin_column2);
	  data_ = _mm_and_si128(new_best, temp);
	  return *this;
	}

	/*
	  score_vector2& update(const score_vector2 &rhs, const score_vector2 &row, const score_vector2 &column, score_vector2 &row_best, score_vector2 &column_best)
	{
	        __m128i temp = _mm_cmplt_epi16(data_, rhs.data_);
		__m128i select_row = _mm_and_si128(row.data_, temp);
		__m128i select_column = _mm_and_si128(column.data_, temp);
		__m128i temp2 = _mm_xor_si128(temp, _mm_set1_epi16(0xffff));
		__m128i origin_row = _mm_and_si128(row_best.data_, temp2);
		__m128i origin_column = _mm_and_si128(column_best.data_, temp2);
		row_best.data_ = _mm_or_si128(select_row, origin_row);
		column_best.data_ = _mm_or_si128(select_column, origin_column);
		data_ = _mm_max_epi16(data_, rhs.data_);
		return *this;
	}
	*/

	friend score_vector2 max(const score_vector2& lhs, const score_vector2 &rhs)
	{
		return score_vector2 (_mm_max_epi8(lhs.data_, rhs.data_));
	}
	
	friend score_vector2 min(const score_vector2& lhs, const score_vector2 &rhs)
	{
		return score_vector2 (_mm_min_epi8(lhs.data_, rhs.data_));
	}

	/*int16_t cmpeq(const score_vector2 &rhs) const
	{
		return _mm_movemask_epi8(_mm_cmpeq_epi16(data_, rhs.data_));
	}

	__m128i cmpeq2(const score_vector2 &rhs) const
	{
		return _mm_cmpeq_epi16(data_, rhs.data_);
	}

	int16_t cmpgt(const score_vector2 &rhs) const
	{
		return _mm_movemask_epi8(_mm_cmpgt_epi16(data_, rhs.data_));
	}*/

	int cmpgt_count()
	{
	  return _mm_movemask_epi8(data_);
	}

	score_vector2 cmpeq(const score_vector2 &rhs) const
        {
	  return score_vector2 (_mm_cmpeq_epi8(data_, rhs.data_));
        }

	score_vector2 cmpgt(const score_vector2 &rhs) const
        {
	        return score_vector2 (_mm_cmpgt_epi8(data_, rhs.data_));
        }

	__m128i data_;

};

#endif /* SCORE_VECTOR2_H_ */
