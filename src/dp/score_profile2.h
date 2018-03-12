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

#ifndef SCORE_PROFILE2_H_
#define SCORE_PROFILE2_H_

#include <vector>
#include "../basic/sequence.h"
#include "score_vector2.h"

#define MASK 24

using std::vector;

struct sequence_stream2
{
	sequence_stream2():
		next (buffer_len),
		mask (0)
	{ }
	inline void reset()
	{
		next = buffer_len;
		mask = 0;
	}
  
        // right
        template<typename _val, typename _score>
	inline const __m128i& get_right(const typename vector<sequence<const _val> >::const_iterator &begin,
					const typename vector<sequence<const _val> >::const_iterator &end,
					unsigned pos,
					const _score&)
        {
		if(next == buffer_len)
		      fill_right<_val,_score>(begin, end, pos);
		return data_[next++];
	}
  
        // left
        template<typename _val, typename _score>
	inline const __m128i& get_left(const typename vector<sequence<const _val> >::const_iterator &begin,
				       const typename vector<sequence<const _val> >::const_iterator &end,
				       unsigned pos,
				       const _score&)
        {
	        if(next == buffer_len)
		  fill_left<_val,_score>(begin, end, pos);
		return data_[next++];
	}

        // right
	template<typename _val, typename _score>
	inline void fill_right(const typename vector<sequence<const _val> >::const_iterator &begin,
			       const typename vector<sequence<const _val> >::const_iterator &end,
			       unsigned pos)
        {
	  //memset(data_, Value_traits<_val>::MASK_CHAR, buffer_len*8);
	  for (unsigned i = 0; i < buffer_len; i++){
	    data_[i] = _mm_set1_epi8(MASK);
	  }
	  
		unsigned n = 0;
		typename vector<sequence<const _val> >::const_iterator it (begin);
	
		assert(pos < it->length());
		const unsigned read_len (std::min(unsigned(buffer_len), static_cast<unsigned>(it->length())-pos));
		while(it < end) {
		        const uint8_t *src (reinterpret_cast<const uint8_t*>(it->data()));
			src += pos;
			_score *dest (reinterpret_cast<_score*>(data_) + n);
			
			if((mask & (1 << n)) == 0) {
			          if(copy_char_right(src, dest, mask, n))
				  if(read_len > 1 && copy_char_right(src, dest, mask, n))
				  if(read_len > 2 && copy_char_right(src, dest, mask, n))
				  if(read_len > 3) copy_char_right(src, dest, mask, n);
			}
			++it;
			++n;
		}
		next = 0;
	}

        // left
        template<typename _val, typename _score>
	inline void fill_left(const typename vector<sequence<const _val> >::const_iterator &begin,
			      const typename vector<sequence<const _val> >::const_iterator &end,
			      unsigned pos)
        {
	  //memset(data_, Value_traits<_val>::MASK_CHAR, buffer_len*8);
	  for (unsigned i = 0; i < buffer_len; i++){
	    data_[i] = _mm_set1_epi8(MASK);
	  }
	       unsigned n = 0;
	       typename vector<sequence<const _val> >::const_iterator it (begin);
	       assert(pos < it->length());
	       unsigned pos2 = it->length() - pos - 1;
	       const unsigned read_len (std::min(unsigned(buffer_len), static_cast<unsigned>(pos2+1)));
	       while(it < end) {
		       const uint8_t *src (reinterpret_cast<const uint8_t*>(it->data()));
		       src += pos2;
		       _score *dest (reinterpret_cast<_score*>(data_) + n);

		       if((mask & (1 << n)) == 0) {
			         if(copy_char_left(src, dest, mask, n))
				 if(read_len > 1 && copy_char_left(src, dest, mask, n))
				 if(read_len > 2 && copy_char_left(src, dest, mask, n))
				 if(read_len > 3) copy_char_left(src, dest, mask, n);
		       }
		       ++it;
		       ++n;
	       }
	       next = 0;
	}
  
        // right
        template<typename _score>
	static inline bool copy_char_right(const uint8_t*& src, _score*& dest, unsigned &mask, unsigned n)
        {
	        if(*src == 0xff) {
		  /*
		    mask |= 1 << n;
		    *dest = 24;
		    dest += 16/sizeof(_score);
		    return true;
		  */
		  mask |= 1 << n;
		  *dest = 24;
		  dest += 16/sizeof(_score);
		  return false;
		}
		*dest = *src & 0x7f;
		++src;
		dest += 16/sizeof(_score);
		return true;
	}

        // left
        template<typename _score>
        static inline bool copy_char_left(const uint8_t*& src, _score*& dest, unsigned &mask, unsigned n)
        { 
	        if(*src == 0xff) {
		  /*
		  mask |= 1 << n;
		  *dest = 24;
		  dest += 16/sizeof(_score);
		  return true;
		  */
		  mask |= 1 << n;
		  *dest = 24;
		  dest += 16/sizeof(_score);
		  return false;
		}
		*dest = *src & 0x7f;
		--src;
		dest += 16/sizeof(_score);
		return true;
	}
  
        static const unsigned buffer_len = 4;
	__m128i data_[buffer_len];
	unsigned next;
	unsigned mask;
};

template<typename _score>
struct score_profile2
{

        template<typename _val>
	inline void set(const __m128i &seq)
        {
		assert(sizeof(data_)/sizeof(score_vector2<_score>) >= Value_traits<_val>::ALPHABET_SIZE);
		/*unsigned j = 0;
		do {
			data_[j] = score_vector<_score> (j, seq);
			++j;
			data_[j] = score_vector<_score> (j, seq);
			++j;
			data_[j] = score_vector<_score> (j, seq);
			++j;
			data_[j] = score_vector<_score> (j, seq);
			++j;
		} while(j<24);
		data_[j] = score_vector<_score> (j, seq);
		assert(j+1 == Value_traits<_val>::ALPHABET_SIZE);*/
		score_vector2<_score> temp;
		for(int j=0;j<Value_traits<_val>::ALPHABET_SIZE;++j){
		        data_[j] = score_vector2<_score> (j, seq);
		}
	}

	template<typename _val>
	inline const score_vector2<_score>& get(_val i) const
	{ return data_[(int)i]; }

        score_vector2<_score> data_[25];
  
};

#endif /* SCORE_PROFILE2_H_ */
