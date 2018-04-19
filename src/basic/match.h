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

#ifndef MATCH_H_
#define MATCH_H_

#include "sequence.h"
#include "../util/async_buffer.h"
#include "edit_transcript.h"

enum Strand { FORWARD, REVERSE };

interval normalized_range(unsigned pos, int len, Strand strand)
{
	return strand == FORWARD
			? interval (pos, pos + len)
			: interval (pos + 1 + len, pos + 1);
}

template<typename _locr, typename _locl>
struct hit
{
	typedef typename packed_sequence_location<_locr>::type packed_loc;

	unsigned	query_;
	packed_loc	subject_;
	_locl		seed_offset_;
	hit():
		query_ (),
		subject_ (),
		seed_offset_ ()
	{ }
	hit(unsigned query, _locr subject, _locl seed_offset):
		query_ (query),
		subject_ (subject),
		seed_offset_ (seed_offset)
	{ }
	bool operator<(const hit &rhs) const
	{ return query_ < rhs.query_; }
	bool blank() const
	{ return subject_ == 0; }
	unsigned operator%(unsigned i) const
	{ return (query_/6) % i; }
	unsigned operator/(unsigned i) const
	{ return (query_/6)/i; }
	int64_t global_diagonal() const
	{ return (int64_t)subject_ - (int64_t)seed_offset_; }
	template<unsigned _d>
	static unsigned query_id(const hit& x)
	{ return x.query_/_d; }
	template<unsigned _d>
	struct Query_id
	{
		unsigned operator()(const hit& x) const
		{ return query_id<_d>(x); }
	};
	static bool cmp_subject(const hit &lhs, const hit &rhs)
	{ return lhs.subject_ < rhs.subject_; }
	static bool cmp_normalized_subject(const hit &lhs, const hit &rhs)
	{
		const uint64_t x = (uint64_t)lhs.subject_ + (uint64_t)rhs.seed_offset_, y = (uint64_t)rhs.subject_ + (uint64_t)lhs.seed_offset_;
		return x < y || (x == y && lhs.seed_offset_ < rhs.seed_offset_);
	}
	friend std::ostream& operator<<(std::ostream &s, const hit &me)
	{
		s << me.query_ << '\t' << me.subject_ << '\t' << me.seed_offset_ << '\n';
		return s;
	}
} __attribute__((packed));

template<typename _val>
struct local_match
{
	typedef typename vector<local_match>::iterator iterator;
	local_match():
		len_ (),
		query_begin_ (),
		subject_len_ (),
		gap_openings_ (),
		identities_ (),
		mismatches_ (),
		subject_begin_ (),
		score_ (),
		query_len_ (),
		query_anchor_ (0),
		subject_ (0),
		subject_position_ (0)
	{ }
	local_match(int score):
                len_ (0),
		query_begin_ (0),
		subject_len_ (0),
		gap_openings_ (0),
		identities_ (0),
		mismatches_ (0),
		subject_begin_ (0),
		score_ (score),
		query_len_ (0),
		query_anchor_ (0),
		subject_ (0),
		subject_position_ (0)
	{ }
        local_match(int query_anchor, const _val *subject):
                len_ (0),
		query_begin_ (0),
		subject_len_ (0),
		gap_openings_ (0),
		identities_ (0),
		mismatches_ (0),
		subject_begin_ (0),
		score_ (0),
		query_len_ (0),
		query_anchor_ (query_anchor),
		subject_ (subject),
		subject_position_ (0)
	{ }
        local_match(int query_anchor, const _val *subject, unsigned long position):
                len_ (0),
		query_begin_ (0),
		subject_len_ (0),
	        gap_openings_ (0),
		identities_ (0),
		mismatches_ (0),
		subject_begin_ (0),
		score_ (0),
		query_len_ (0),
		query_anchor_ (query_anchor),
		subject_ (subject),
		subject_position_ (position)
        { }
        local_match(unsigned len, unsigned query_begin, unsigned query_len, unsigned subject_len, unsigned gap_openings, unsigned identities, unsigned mismatches, signed subject_begin, signed score, int query_anchor, unsigned long position):
		len_ (len),
		query_begin_ (query_begin),
		subject_len_ (subject_len),
		gap_openings_ (gap_openings),
		identities_ (identities),
		mismatches_ (mismatches),
		subject_begin_ (subject_begin),
		score_ (score),
		query_len_ (query_len),
		query_anchor_ (query_anchor),
		subject_ (0),
		subject_position_ (position)
        { }
	local_match& operator+=(const local_match& rhs)
	{
		add(rhs);
		transcript_right_ = rhs.transcript_right_;
		return *this;
	}
	local_match& operator-=(const local_match& rhs)
	{
		add(rhs);
		query_begin_ = rhs.query_len_;
		subject_begin_ = rhs.subject_len_;
		transcript_left_ = rhs.transcript_right_;
		return *this;
	}
	void add(const local_match &rhs)
	{
		len_ += rhs.len_;
		subject_len_ += rhs.subject_len_;
		gap_openings_ += rhs.gap_openings_;
		identities_ += rhs.identities_;
		mismatches_ += rhs.mismatches_;
		score_ += rhs.score_;
		query_len_ += rhs.query_len_;
	}
	interval query_range(Strand strand) const
	{ return normalized_range(query_begin_, query_len_, strand); }
	interval subject_range() const
	{ return normalized_range(subject_begin_, subject_len_, FORWARD); }
	void print(const sequence<const _val> &query, const sequence<const _val> &subject, const vector<char> &buf) const
	{
		cout << "Score = " << score_ << endl;
		::print(cout, &query[query_begin_], &subject[subject_begin_], transcript_right_, transcript_left_, buf);
	}
	/*friend std::ostream& operator<<(std::ostream &os, const local_match &x)
	{
		os << "(sbj=" << x.subject_range() << " score=" << Scoring<_val>::bitscore(x.score_) << ")";
		return os;
	}*/
	unsigned len_, query_begin_, subject_len_, gap_openings_, identities_, mismatches_;
	signed subject_begin_, score_, query_len_, query_anchor_;
	const _val *subject_;
        unsigned long subject_position_;
	Edit_transcript transcript_right_, transcript_left_;
};

template<typename _val>
struct Segment
{
	Segment(int score,
		unsigned frame,
		local_match<_val> *traceback = 0,
		unsigned subject_id = std::numeric_limits<unsigned>::max()):
		score_ (score),
		frame_ (frame),
		traceback_ (traceback),
		subject_id_ (subject_id),
		next_ (0),
		top_score_ (0)
	{ }
	Strand strand() const
	{ return frame_ < 3 ? FORWARD : REVERSE; }
	interval query_range() const
	{ return traceback_->query_range(strand()); }
	interval subject_range() const
	{ return traceback_->subject_range(); }
	bool operator<(const Segment &rhs) const
	{ return top_score_ > rhs.top_score_
			|| (top_score_ == rhs.top_score_
			&& (subject_id_ < rhs.subject_id_ || (subject_id_ == rhs.subject_id_ && (score_ > rhs.score_ || (score_ == rhs.score_ && traceback_->score_ > rhs.traceback_->score_))))); }
	static bool comp_subject(const Segment& lhs, const Segment &rhs)
	{ return lhs.subject_id_ < rhs.subject_id_ || (lhs.subject_id_ == rhs.subject_id_ && lhs.score_ > rhs.score_); }
	struct Subject
	{
		unsigned operator()(const Segment& x) const
		{ return x.subject_id_; }
	};
        int						score_;
	unsigned				frame_;
        local_match<_val>	   *traceback_;
	unsigned				subject_id_;
	Segment				   *next_;
	int						top_score_;
};

#endif /* MATCH_H_ */
