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

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "options.h"
#include "value_type.h"

namespace program_options {
string		input_ref_file;
uint32_t	threads_;
string		database;
string		query_file;
uint32_t	merge_seq_treshold;
uint32_t	block_size;
uint32_t	hit_cap;
int			min_ungapped_raw_score;
uint32_t	shapes;
uint32_t	index_mode;
uint64_t	max_alignments;
string		output_file;
string		match_file1;
string		match_file2;
int			padding;
uint32_t	output_threads;
uint32_t	compression;
double		chunk_size;
double          query_chunk_size;
unsigned	min_identities;
unsigned	min_identities2;
int			xdrop;
unsigned	window;
int			min_hit_score;
int			hit_band;
unsigned	min_compressed_identities;
int			min_seed_score;
unsigned	seed_signatures;
double		min_bit_score;
unsigned	run_len;
bool		alignment_traceback;
double		max_seed_freq;
string		tmpdir;
bool		long_mode;
int			gapped_xdrop;
double		max_evalue;
string		sam_output;
string		kegg_file;
int			gap_open;
int			gap_extend;
string		matrix;
string		seg;
bool		verbose;
bool		debug_log;
bool		have_ssse3;
bool		salltitles;
int			reward;
int			penalty;
string		db_type;
double		min_id;
unsigned	compress_temp;
double		toppercent;
string		daa_file;
string		output_format;
bool		forwardonly;
unsigned	fetch_size;
bool		single_domain;

Aligner_mode aligner_mode;
Command command;

template<typename _val>
void set_options(double block_size)
{
  	if(aligner_mode == sensitive) {
		set_option(seed_signatures, 1u);
		set_option(index_mode, 2u);
	} else if (aligner_mode == fast) {
		set_option(seed_signatures, 1u);
		set_option(index_mode, 1u);
	}

	set_option(chunk_size, block_size);

}

string get_temp_file()
{
	if(strlen(getenv("TMPDIR")) > 0)
		return string(getenv("TMPDIR")) + "/ac-diamond.tmp";
	else {
		std::cerr << "Warning: TMPDIR environment variable not set - using output directory for temporary storage.\n";
		return output_file + ".tmp";
	}
}

template<typename _val>
unsigned read_padding(size_t len)
{
	if(padding == 0) {
		if(len<=255)
			return 10;
		else
			return 32;
	} else
		return padding;
}

template<>
unsigned read_padding<Amino_acid>(size_t len)
{
	if(padding == 0) {
		if(len<=35)
			return 5;
		else if(len<=55)
			return 16;
		else
			return 32;
	} else
		return padding;
}

template void set_options<Amino_acid>(double block_size);
template void set_options<Nucleotide>(double block_size);
template unsigned read_padding<Nucleotide>(size_t len);
template unsigned read_padding<Amino_acid>(size_t len);

bool mem_buffered()
{ return tmpdir == "/dev/shm"; }

}
