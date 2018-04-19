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

#ifndef ALIGN_UNGAPPED_H_
#define ALIGN_UNGAPPED_H_

#define SEED_LEN2 (unsigned)10

using std::max;
using std::min;

template<typename _val, typename _locr>
int xdrop_ungapped(const _val *query, const _val *subject)
{
        int score, st;
	unsigned left (0), right (0);
	const _val *q;
	const _val *s;
	
	score = 0;
	for(unsigned i = 0; i < SEED_LEN2; i++)
	  score += score_matrix::get().letter_score(query[i], mask_critical(subject[i]));

	st = score;
	
	q = query + SEED_LEN2;
	s = subject + SEED_LEN2;
	
	const unsigned window_right = max(program_options::window, SEED_LEN2 - Const::seed_anchor) - SEED_LEN2 + Const::seed_anchor;
	while(score - st < program_options::xdrop
	      && right < window_right
	      && *q != String_set<_val>::PADDING_CHAR
	      && *s != String_set<_val>::PADDING_CHAR)
	  {
	        st += score_matrix::get().letter_score(*q, mask_critical(*s));
		score = max(score, st);
		++q;
		++s;
		++right;
	}
	
	if (score >= program_options::min_hit_score) return score;
	
	q = query - 1;
	s = subject - 1;
	
        const unsigned window_left = max(program_options::window, (unsigned)Const::seed_anchor) - Const::seed_anchor;
	
	st = score;
	
	while(score - st < program_options::xdrop
	      && left < window_left
	      && *q != String_set<_val>::PADDING_CHAR
              && *s != String_set<_val>::PADDING_CHAR)
	{
	        st += score_matrix::get().letter_score(*q, mask_critical(*s));
		score = max(score, st);
		--q;
		--s;
		++left;
	}
	
	return score;
}

#endif /* ALIGN_UNGAPPED_H_ */
