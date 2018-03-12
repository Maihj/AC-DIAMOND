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

#ifndef CONST_H_
#define CONST_H_

struct Const
{

	enum {
		build_version = 1,
		build_compatibility = 1,
		seedp_bits = 10,
                seedp = 1<<seedp_bits,
		seqp_bits = 8,
                seqp = 1<<seqp_bits,
		daa_version = 0,
		max_shapes = 16,
		index_modes = 2,
		min_shape_len = 8,
		max_shape_len = 32,
		seed_anchor = 8
	};

	static const char* version_string;
	static const char* program_name;
	static const char* id_delimiters;

};

const char* Const::version_string = "v1";
const char* Const::program_name = "ac-diamond";
const char* Const::id_delimiters = " \a\b\f\n\r\t\v";

#endif /* CONST_H_ */
