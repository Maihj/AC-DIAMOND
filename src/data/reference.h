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

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include <memory>
#include <string>
#include <numeric>
#include "../util/binary_file.h"
#include "../basic/statistics.h"
#include "../util/hash_table.h"
#include "../util/hash_function.h"
#include "../util/util.h"
#include "sequence_set.h"
#include "boost/ptr_container/ptr_vector.hpp"

using std::auto_ptr;
using boost::ptr_vector;

struct Reference_header
{
	Reference_header():
		unique_id (0x24af8a415ee186dllu),
		build (Const::build_version),
		long_addressing (false),
		sequences (0),
		letters (0)
	{ }
	uint64_t unique_id;
	uint32_t build;
	bool long_addressing;
	unsigned n_blocks;
	size_t sequences, letters;
	double block_size;
#ifdef EXTRA
	Sequence_type sequence_type;
#endif
} ref_header;

struct Database_format_exception : public exception
{
	virtual const char* what() const throw()
	{ return "Database file is not a AC-DIAMOND database."; }
};

template<typename _val>
struct Database_file : public Input_stream
{
	Database_file():
		Input_stream (program_options::database)
	{
		if(this->read(&ref_header, 1) != 1)
			throw Database_format_exception ();
		if(ref_header.unique_id != Reference_header ().unique_id)
			throw Database_format_exception ();
		if(ref_header.build > Const::build_version || ref_header.build < Const::build_compatibility)
			throw invalid_database_version_exception();
#ifdef EXTRA
		if(sequence_type(_val()) != ref_header.sequence_type)
			throw std::runtime_error("Database has incorrect sequence type for current alignment mode.");
#endif
	}
	void rewind()
	{ this->seek(sizeof(Reference_header)); }
};

template<typename _val>
struct ref_seqs
{
	static const Sequence_set<_val>& get()
	{ return *data_; }
	static Sequence_set<_val> *data_;
};

template<typename _val> Sequence_set<_val>* ref_seqs<_val>::data_ = 0;

struct ref_ids
{
	static const String_set<char,0>& get()
	{ return *data_; }
	static String_set<char,0> *data_;
};

String_set<char,0>* ref_ids::data_ = 0;

unsigned current_ref_block;

struct Ref_map
{
	Ref_map():
		next_ (0)
	{ }
	void init(unsigned ref_count)
	{
		const unsigned block = current_ref_block;
		if(data_.size() < block+1) {
			data_.resize(block+1);
			data_[block].insert(data_[block].end(), ref_count, std::numeric_limits<unsigned>::max());
		}
	}
	template<typename _val>
	uint32_t get(unsigned block, unsigned i)
	{
		uint32_t n = data_[block][i];
		if(n != std::numeric_limits<unsigned>::max())
			return n;
		else {
			tthread::lock_guard<tthread::mutex> lock (mtx_);
			n = data_[block][i];
			if(n != std::numeric_limits<uint32_t>::max())
				return n;
			n = next_++;
			data_[block][i] = n;
			len_.push_back(ref_seqs<_val>::get().length(i));
			name_.push_back(get_str(ref_ids::get()[i].c_str(), Const::id_delimiters));
			return n;
		}
	}
	/*template<typename _val>
	void finish()
	{
		vector<pair<unsigned,unsigned> > v;
		for(unsigned i=0;i<data_.size();++i)
			if(data_[i] != std::numeric_limits<unsigned>::max())
				v.push_back(pair<unsigned,unsigned> (data_[i], i));
		std::sort(v.begin(), v.end());
		for(vector<pair<unsigned,unsigned> >::const_iterator i = v.begin(); i!=v.end(); ++i) {
			const char* s = ref_ids::get()[i->second].c_str();
			buf_ << (uint32_t)(strlen(s)+1);
			buf_.write_c_str(s);
			buf_ << (uint32_t)ref_seqs<_val>::get()[i->second].length();
		}
	}*/
private:
	tthread::mutex mtx_;
	vector<vector<uint32_t> > data_;
	vector<uint32_t> len_;
	ptr_vector<string> name_;
	uint32_t next_;
	friend struct DAA_output;
} ref_map;

#endif /* REFERENCE_H_ */
