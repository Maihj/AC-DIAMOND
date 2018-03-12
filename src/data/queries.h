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

#ifndef QUERIES_H_
#define QUERIES_H_

#include "../basic/translate.h"
#include "../util/complexity_filter.h"
#include "../basic/statistics.h"
#include "sequence_set.h"

unsigned current_query_chunk;

struct query_source_seqs
{
	static const Sequence_set<Nucleotide>& get()
	{ return *data_; }
	static Sequence_set<Nucleotide> *data_;
};

Sequence_set<Nucleotide>* query_source_seqs::data_ = 0;

template<typename _val>
struct query_seqs
{
	static const Sequence_set<_val>& get()
	{ return *data_; }
	static Sequence_set<_val> *data_;
};

template<typename _val> Sequence_set<_val>* query_seqs<_val>::data_ = 0;

struct query_ids
{
	static const String_set<char,0>& get()
	{ return *data_; }
	static String_set<char,0> *data_;
};

String_set<char,0>* query_ids::data_ = 0;

#endif /* QUERIES_H_ */
