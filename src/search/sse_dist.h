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

#ifndef SSE_DIST_H_
#define SSE_DIST_H_

#ifdef __SSSE3__
#include <tmmintrin.h>
#include <smmintrin.h>
#include <emmintrin.h>
#endif

#include "../basic/reduction.h"

unsigned POPCOUNT_4bit[16] __attribute__ ((aligned)) = {
   /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
   /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
   /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
   /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4
 };

unsigned POPCOUNT_8bit[256] __attribute__ ((aligned)) = {
   /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
   /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
   /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
   /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4,
   /* 10 */ 1, /* 11 */ 2, /* 12 */ 2, /* 13 */ 3,
   /* 14 */ 2, /* 15 */ 3, /* 16 */ 3, /* 17 */ 4,
   /* 18 */ 2, /* 19 */ 3, /* 1a */ 3, /* 1b */ 4,
   /* 1c */ 3, /* 1d */ 4, /* 1e */ 4, /* 1f */ 5,
   /* 20 */ 1, /* 21 */ 2, /* 22 */ 2, /* 23 */ 3,
   /* 24 */ 2, /* 25 */ 3, /* 26 */ 3, /* 27 */ 4,
   /* 28 */ 2, /* 29 */ 3, /* 2a */ 3, /* 2b */ 4,
   /* 2c */ 3, /* 2d */ 4, /* 2e */ 4, /* 2f */ 5,
   /* 30 */ 2, /* 31 */ 3, /* 32 */ 3, /* 33 */ 4,
   /* 34 */ 3, /* 35 */ 4, /* 36 */ 4, /* 37 */ 5,
   /* 38 */ 3, /* 39 */ 4, /* 3a */ 4, /* 3b */ 5,
   /* 3c */ 4, /* 3d */ 5, /* 3e */ 5, /* 3f */ 6,
   /* 40 */ 1, /* 41 */ 2, /* 42 */ 2, /* 43 */ 3,
   /* 44 */ 2, /* 45 */ 3, /* 46 */ 3, /* 47 */ 4,
   /* 48 */ 2, /* 49 */ 3, /* 4a */ 3, /* 4b */ 4,
   /* 4c */ 3, /* 4d */ 4, /* 4e */ 4, /* 4f */ 5,
   /* 50 */ 2, /* 51 */ 3, /* 52 */ 3, /* 53 */ 4,
   /* 54 */ 3, /* 55 */ 4, /* 56 */ 4, /* 57 */ 5,
   /* 58 */ 3, /* 59 */ 4, /* 5a */ 4, /* 5b */ 5,
   /* 5c */ 4, /* 5d */ 5, /* 5e */ 5, /* 5f */ 6,
   /* 60 */ 2, /* 61 */ 3, /* 62 */ 3, /* 63 */ 4,
   /* 64 */ 3, /* 65 */ 4, /* 66 */ 4, /* 67 */ 5,
   /* 68 */ 3, /* 69 */ 4, /* 6a */ 4, /* 6b */ 5,
   /* 6c */ 4, /* 6d */ 5, /* 6e */ 5, /* 6f */ 6,
   /* 70 */ 3, /* 71 */ 4, /* 72 */ 4, /* 73 */ 5,
   /* 74 */ 4, /* 75 */ 5, /* 76 */ 5, /* 77 */ 6,
   /* 78 */ 4, /* 79 */ 5, /* 7a */ 5, /* 7b */ 6,
   /* 7c */ 5, /* 7d */ 6, /* 7e */ 6, /* 7f */ 7,
   /* 80 */ 1, /* 81 */ 2, /* 82 */ 2, /* 83 */ 3,
   /* 84 */ 2, /* 85 */ 3, /* 86 */ 3, /* 87 */ 4,
   /* 88 */ 2, /* 89 */ 3, /* 8a */ 3, /* 8b */ 4,
   /* 8c */ 3, /* 8d */ 4, /* 8e */ 4, /* 8f */ 5,
   /* 90 */ 2, /* 91 */ 3, /* 92 */ 3, /* 93 */ 4,
   /* 94 */ 3, /* 95 */ 4, /* 96 */ 4, /* 97 */ 5,
   /* 98 */ 3, /* 99 */ 4, /* 9a */ 4, /* 9b */ 5,
   /* 9c */ 4, /* 9d */ 5, /* 9e */ 5, /* 9f */ 6,
   /* a0 */ 2, /* a1 */ 3, /* a2 */ 3, /* a3 */ 4,
   /* a4 */ 3, /* a5 */ 4, /* a6 */ 4, /* a7 */ 5,
   /* a8 */ 3, /* a9 */ 4, /* aa */ 4, /* ab */ 5,
   /* ac */ 4, /* ad */ 5, /* ae */ 5, /* af */ 6,
   /* b0 */ 3, /* b1 */ 4, /* b2 */ 4, /* b3 */ 5,
   /* b4 */ 4, /* b5 */ 5, /* b6 */ 5, /* b7 */ 6,
   /* b8 */ 4, /* b9 */ 5, /* ba */ 5, /* bb */ 6,
   /* bc */ 5, /* bd */ 6, /* be */ 6, /* bf */ 7,
   /* c0 */ 2, /* c1 */ 3, /* c2 */ 3, /* c3 */ 4,
   /* c4 */ 3, /* c5 */ 4, /* c6 */ 4, /* c7 */ 5,
   /* c8 */ 3, /* c9 */ 4, /* ca */ 4, /* cb */ 5,
   /* cc */ 4, /* cd */ 5, /* ce */ 5, /* cf */ 6,
   /* d0 */ 3, /* d1 */ 4, /* d2 */ 4, /* d3 */ 5,
   /* d4 */ 4, /* d5 */ 5, /* d6 */ 5, /* d7 */ 6,
   /* d8 */ 4, /* d9 */ 5, /* da */ 5, /* db */ 6,
   /* dc */ 5, /* dd */ 6, /* de */ 6, /* df */ 7,
   /* e0 */ 3, /* e1 */ 4, /* e2 */ 4, /* e3 */ 5,
   /* e4 */ 4, /* e5 */ 5, /* e6 */ 5, /* e7 */ 6,
   /* e8 */ 4, /* e9 */ 5, /* ea */ 5, /* eb */ 6,
   /* ec */ 5, /* ed */ 6, /* ee */ 6, /* ef */ 7,
   /* f0 */ 4, /* f1 */ 5, /* f2 */ 5, /* f3 */ 6,
   /* f4 */ 5, /* f5 */ 6, /* f6 */ 6, /* f7 */ 7,
   /* f8 */ 5, /* f9 */ 6, /* fa */ 6, /* fb */ 7,
   /* fc */ 6, /* fd */ 7, /* fe */ 7, /* ff */ 8
 };

unsigned popcount_3(uint64_t x)
{
	const uint64_t m1  = 0x5555555555555555; //binary: 0101...
	const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
	const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
	const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits
    return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
}

template<typename _val>
unsigned match_block(const _val *x, const _val *y)
{
        static const __m128i mask = _mm_set1_epi8(0x7F);
	__m128i r1 = _mm_load_si128 ((__m128i const*)(x));
	__m128i r2 = _mm_load_si128 ((__m128i const*)(y));
	r2 =  _mm_and_si128(r2, mask);
	return __builtin_popcount(_mm_movemask_epi8(_mm_cmpeq_epi8(r1, r2)));
}

template<typename _val>
unsigned fast_match(const _val *q, const _val *s)
{ return __builtin_popcount(match_block(q-8, s-8) << 16 | match_block(q+8, s+8)); }
//{ return popcount_3(match_block(q-8, s-8)<<16 | match_block(q+8, s+8)); }

template<typename _val>
__m128i reduce_seq_ssse3(const __m128i &seq)
{
#ifdef __SSSE3__
	const __m128i *row = reinterpret_cast<const __m128i*>(Reduction<_val>::reduction.map8());
	__m128i high_mask = _mm_slli_epi16(_mm_and_si128(seq, _mm_set1_epi8(0x10)), 3);
	__m128i seq_low = _mm_or_si128(seq, high_mask);
	__m128i seq_high = _mm_or_si128(seq, _mm_xor_si128(high_mask, _mm_set1_epi8(0x80)));

	__m128i r1 = _mm_load_si128(row);
	__m128i r2 = _mm_load_si128(row+1);
	__m128i s1 = _mm_shuffle_epi8(r1, seq_low);
	__m128i s2 = _mm_shuffle_epi8(r2, seq_high);
	return _mm_or_si128(s1, s2);
#endif
}

template<typename _val>
__m128i reduce_seq_generic(const __m128i &seq)
{
	__m128i r;
	uint8_t* s = (uint8_t*)&seq;
	uint8_t* d = (uint8_t*)&r;
	for(unsigned i=0;i<16;++i){
	  if (*s == 0x7F) *(d++) = 0x7F;
	  else *(d++) = Reduction<_val>::reduction(*s);
	  s++;
	}
	return r;
}

template<typename _val>
__m128i reduce_seq(const __m128i &seq)
{
	if(program_options::have_ssse3) {
#ifdef __SSSE3__
		return reduce_seq_ssse3<_val>(seq);
#else
		return reduce_seq_generic<_val>(seq);
#endif
	} else
		return reduce_seq_generic<_val>(seq);
}

template<typename _val>
unsigned match_block_reduced(const _val *x, const _val *y)
{
        static const __m128i mask = _mm_set1_epi8(0x7F);
	__m128i r1 = _mm_load_si128 ((__m128i const*)(x));
	__m128i r2 = _mm_load_si128 ((__m128i const*)(y));
	r1 = _mm_and_si128(r1, mask);
	r2 = _mm_and_si128(r2, mask);
	r1 = reduce_seq_generic<_val>(r1);
	r2 = reduce_seq_generic<_val>(r2);
	return __builtin_popcount(_mm_movemask_epi8(_mm_cmpeq_epi8(r1, r2)));
}

template<typename _val>
uint64_t reduced_match32(const _val* q, const _val *s, unsigned len)
{
	uint64_t x = match_block_reduced(q+16, s+16)<<16 | match_block_reduced(q,s);
	if(len < 32)
		x &= (1 << len) - 1;
	return x;
}

#endif /* SSE_DIST_H_ */
