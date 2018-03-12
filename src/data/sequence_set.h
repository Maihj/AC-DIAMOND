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

#ifndef SEQUENCE_SET_H_
#define SEQUENCE_SET_H_

#include <iostream>
#include <string>
#include "../basic/sequence.h"
#include "string_set.h"

using std::cout;
using std::endl;
using std::pair;

template<typename _val>
struct Sequence_set : public String_set<_val>
{

	Sequence_set()
	{ }

	Sequence_set(Input_stream &file):
		String_set<_val> (file)
	{ }

	void print_stats() const
	{ verbose_stream << "Sequences = " << this->get_length() << ", letters = " << this->letters() << endl; }

	pair<size_t,size_t> len_bounds(size_t min_len) const
	{
		const size_t l (this->get_length());
		size_t max = 0, min = std::numeric_limits<size_t>::max();
		for(size_t i=0;i<l;++i) {
			max = std::max(this->length(i), max);
			min = this->length(i) >= min_len ? std::min(this->length(i), min) : min;
		}
		return pair<size_t,size_t> (min, max);
	}

	sequence<const _val> window_infix(size_t offset, unsigned &left) const
	{
		const _val* begin (this->data(offset));
		unsigned n (0);
		while(*begin != String_set<_val>::PADDING_CHAR && n <= program_options::window) {
			--begin;
			++n;
		}
		++begin;
		left = program_options::window + 1 - n;
		const _val* end (this->data(offset));
		n = 0;
		while(*end != String_set<_val>::PADDING_CHAR && n < program_options::window) {
			++end;
			++n;
		}
		return sequence<const _val> (begin, end - begin);
	}
	
	sequence<const _val> window_infix_query(size_t offset, unsigned &left) const
	{
	    const _val* begin (this->data(offset));
	    unsigned n (0);
	    while(*begin != String_set<_val>::PADDING_CHAR && n <= program_options::window) {
	      --begin;
	      ++n;
	    }
	    ++begin;
	    left = program_options::window + 1 - n;
	    const _val* end (this->data(offset));
	    n = 0;
	    while(*end != String_set<_val>::PADDING_CHAR && n < program_options::window) {
	      ++end;
	      ++n;
	    }
	    return sequence<const _val> (begin, end - begin);
	  }

	sequence<const _val> fixed_window_infix(size_t offset) const
	{
		const _val* begin (this->data(offset));
		unsigned n (0);
		while(*begin != String_set<_val>::PADDING_CHAR && n <= program_options::window) {
			--begin;
			++n;
		}
		++begin;
		const _val* s (this->data(offset - program_options::window));
		return sequence<const _val> (s, 2*program_options::window, begin - s);
	}
	
	sequence<const _val> window_infix_dp_right(size_t offset, unsigned len) const
	{
	        const _val* begin (this->data(offset));
		return sequence<const _val> (begin, len, 0);
	}

	sequence<const _val> window_infix_dp_left(size_t offset, unsigned len) const
	{
	        const _val* begin (this->data(offset));
		unsigned n (0);
		while(*begin != String_set<_val>::PADDING_CHAR && n < len) {
		        --begin;
			++n;
		}
		++begin;
		const _val* s (this->data(offset - len + 1));
		return sequence<const _val> (s, len, begin - s);
	}

	vector<size_t> partition() const
	{
	        vector<size_t> v;
		const size_t l = (this->letters()+Const::seqp-1) / Const::seqp;
		v.push_back(0);
		for(unsigned i=0;i<this->get_length();) {
			size_t n = 0;
			while(i<this->get_length() && n < l)
				n += this->length(i++);
			v.push_back(i);
		}
		for(unsigned i=v.size();i<Const::seqp+1;++i)
			v.push_back(this->get_length());
		return v;
	}

	size_t reverse_translated_len(size_t i) const
	{
		const size_t j (i - i%6);
		const size_t l (this->length(j));
		if(this->length(j+2) == l)
			return l*3 + 2;
		else if(this->length(j+1) == l)
			return l*3 + 1;
		else
			return l*3;
	}

};

#endif /* SEQUENCE_SET_H_ */
