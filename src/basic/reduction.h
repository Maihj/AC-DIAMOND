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

#ifndef REDUCTION_H_
#define REDUCTION_H_

using std::string;
using std::vector;

#include "value.h"

template<typename _val>
struct Reduction
{

    Reduction(const char *definition_string)
	{
		memset(map_, 0, sizeof(map_));
		memset(map8_, 0, sizeof(map8_));
		const vector<string> tokens (tokenize(definition_string, " "));
		size_ = tokens.size();
		for(unsigned i=0;i<size_;++i)
		  for(unsigned j=0;j<tokens[i].length();++j) {
		    const char ch = tokens[i][j];
		    map_[(long) Value_traits<_val>::from_char(ch)] =  i;
		    map8_[(long) Value_traits<_val>::from_char(ch)] =  i;
		  }
	}

	unsigned size() const
	{ return size_; }

	unsigned operator()(_val a) const
	{ return map_[(long)a]; }

	const char* map8() const
	{ return map8_; }

	static const Reduction reduction;

private:

	unsigned map_[256];
	char map8_[256];
	unsigned size_;

};

// alphabet 11
template<> const Reduction<Amino_acid> Reduction<Amino_acid>::reduction ("KREDQN C G H M F Y ILV W P STA BJZXUO");
// alphabet 13
//template<> const Reduction<Amino_acid> Reduction<Amino_acid>::reduction ("C G H M F Y W P KR ED QN ILV STA BJZXUO");
template<> const Reduction<Nucleotide> Reduction<Nucleotide>::reduction ("A C G T");

#ifdef EXTRA
#include "../../../extra/reduction.h"
#endif

#endif /* REDUCTION_H_ */
