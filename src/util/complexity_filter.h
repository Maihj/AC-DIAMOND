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

#ifndef COMPLEXITY_FILTER_H_
#define COMPLEXITY_FILTER_H_

#include "../algo/blast/core/blast_seg.h"
#include "../algo/blast/core/blast_filter.h"
#include "../basic/value.h"
#include "thread.h"

template<class _val>
struct Complexity_filter
{
	unsigned filter(vector<_val> &seq) const
	{ return 0; }
	static const Complexity_filter& get()
	{ return instance; }
	void run(String_set<_val> &seqs) const
	{ }
	static const Complexity_filter instance;
};

template<>
struct Complexity_filter<Amino_acid>
{

	Complexity_filter()
	{ blast_seg_ = SegParametersNewAa(); }

	~Complexity_filter()
	{ SegParametersFree(blast_seg_); }

	unsigned filter(sequence<Amino_acid> seq) const
	{
		BlastSeqLoc *seg_locs;
		SeqBufferSeg ((uint8_t*) seq.data(), seq.length(), 0, blast_seg_, &seg_locs);
		unsigned nMasked = 0;

		if(seg_locs) {
			BlastSeqLoc *l = seg_locs;
			do {
				for(signed i=l->ssr->left;i<=l->ssr->right;i++) {
					nMasked++;
					seq[i] = Value_traits<Amino_acid>::MASK_CHAR;
				}
			} while((l=l->next) != 0);
			BlastSeqLocFree(seg_locs);
		}
		return nMasked;
	}

	static const Complexity_filter& get()
	{ return instance; }

	void run(String_set<Amino_acid> &seqs) const
	{
		Filter_context context (seqs, *this);
		launch_scheduled_thread_pool(context, seqs.get_length(), program_options::threads());
	}

private:

	struct Filter_context
	{
		Filter_context(String_set<Amino_acid> &seqs, const Complexity_filter &filter):
			seqs (seqs),
			filter (filter)
		{ }
		void operator()(unsigned thread_id, unsigned i)
		{
			filter.filter(seqs[i]);
		}
		String_set<Amino_acid> &seqs;
		const Complexity_filter &filter;
	};

	SegParameters *blast_seg_;

	static const Complexity_filter instance;

};

const Complexity_filter<Amino_acid> Complexity_filter<Amino_acid>::instance;
#ifdef NDEBUG
template<> const Complexity_filter<Nucleotide> Complexity_filter<Nucleotide>::instance;
#else
template<typename _val> const Complexity_filter<_val> Complexity_filter<_val>::instance;
#endif


#endif /* COMPLEXITY_FILTER_H_ */
