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

#ifndef MAKE_DB_H_
#define MAKE_DB_H_

#include <iostream>
#include "../basic/options.h"
#include "../data/reference.h"
#include "../basic/exceptions.h"
#include "../basic/statistics.h"
#include "../basic/reduction.h"
#include "../data/load_seqs.h"
#include "../util/seq_file_format.h"

#define SIZE 11
#define SEED_LEN 10
#define MAX_SEED_LEN 16
#define CANDIDATES_FAST 200
#define CANDIDATES_SENSITIVE 300
#define REF_TABLE_SIZE 16777216
#define MAX_COUNT 65535
#define MAX_BS_CONST 20
#define BS_BITS 24
#define	HT_BITS 16
#define N_PART 3
#define N_PATTERN 9
#define OFFSET 6

//using namespace __gnu_cxx;
using boost::timer::cpu_timer;
using std::cout;
using std::endl;
using std::max;
using std::min;
using std::make_pair;

// fast mode: the first pattern; sensitive mode: 2nd~9th pattern (totally 8 patterns)
uint64_t MIS_PATTERN[N_PATTERN] = {0xfff0ffff0fffffff, 0xfff00f0ff0ffffff, 0xff0fff0f00ffffff, 0xf00ff0fff0ffffff, 0xff0ff0f0f0ffffff, 0xff0f0f000fffffff, 0xff0ff00f00ffffff, 0xf0ff0ff000ffffff, 0xf0000f0fffffffff };

bool cmp_keys(const pair<uint64_t, pair<unsigned, uint8_t> > &lhs, const pair<uint64_t, pair<unsigned, uint8_t> > &rhs){
  return lhs.first < rhs.first;
}

template<typename _locr>
class Compare {
  
 public:
  Compare (uint64_t p) : mis_pattern(p){ }
  bool operator() (const pair<uint64_t, _locr> &lhs, const pair<uint64_t, _locr> &rhs)
  {
    uint64_t l, r;
    l = lhs.first & mis_pattern;
    r = rhs.first & mis_pattern;
    return (l < r) || (l == r && lhs.first < rhs.first);
  }
  
  uint64_t mis_pattern;
};

class Compare_long {

 public:
  Compare_long (uint64_t p) : mis_pattern(p){ }
  bool operator() (const pair<uint64_t, pair<unsigned, uint8_t> > &lhs, const pair<uint64_t, pair<unsigned, uint8_t> > &rhs)
  {
    uint64_t l, r;
    l = lhs.first & mis_pattern;
    r = rhs.first & mis_pattern;
    return (l < r) || (l == r && lhs.first < rhs.first);
  }
  
  uint64_t mis_pattern;
};

void calculate_factor(uint64_t base[MAX_SEED_LEN]){
  unsigned i;
  uint64_t init;
  init = 0xf;
  for(i = 0; i < MAX_SEED_LEN; i++){
    base[i] = init;
    init = init << 4;
  }
}

template<typename _locr>
void store_to_file_long(unsigned id,
			unsigned *msa_4bytes,
			uint8_t *msa_1byte,
			_locr msa_len,
			vector<pair<uint64_t, pair<pair<unsigned, uint8_t>, unsigned> > > &pattern_ids64,
			const char *filename){
  FILE *f;
  if (id == 0){
    f = fopen(filename, "wb");
  }
  else {
    f = fopen(filename, "ab");
  }
  
  task_timer timer ("Building the hash table", true);

  _locr i;
  _locr len = pattern_ids64.size();
  unsigned hashValue;
  
  unsigned *table = (unsigned*)malloc(REF_TABLE_SIZE * sizeof(unsigned)); // 24 bits
  uint16_t *count = (uint16_t*)malloc(REF_TABLE_SIZE * sizeof(uint16_t));
  memset(table, 0, REF_TABLE_SIZE * sizeof(unsigned));
  memset(count, 0, REF_TABLE_SIZE * sizeof(uint16_t));
  unsigned *ref_sal_4bytes = (unsigned*)malloc(len * sizeof(unsigned));
  uint8_t *ref_sal_1byte = (uint8_t*)malloc(len * sizeof(uint8_t));
  unsigned *ref_sar = (unsigned*)malloc(len * sizeof(unsigned));
  uint16_t *key16 = (uint16_t*)malloc(len * sizeof(uint16_t));
  
  // hash table
  hashValue = (pattern_ids64[0].first >> HT_BITS) & 0xffffff;
  table[hashValue] = 0;
  count[hashValue] = 1;
  ref_sal_4bytes[0] = pattern_ids64[0].second.first.first;
  ref_sal_1byte[0] = pattern_ids64[0].second.first.second;
  ref_sar[0] = pattern_ids64[0].second.second;
  key16[0] = (uint16_t)(pattern_ids64[0].first & 0xffff);
  
  for (i = 1; i < len; i++){
    if (((pattern_ids64[i].first >> HT_BITS) & 0xffffff) != hashValue){
      hashValue = (pattern_ids64[i].first >> HT_BITS) & 0xffffff;
      table[hashValue] = i;
      count[hashValue] = 1;
    }
    else {
      count[hashValue]++;
    }
    
    ref_sal_4bytes[i] = pattern_ids64[i].second.first.first;
    ref_sal_1byte[i] = pattern_ids64[i].second.first.second;
    ref_sar[i] = pattern_ids64[i].second.second;
    key16[i] = (uint16_t)(pattern_ids64[i].first & 0xffff);
  }

  timer.go("Writing the index to disk");
  
  fwrite(&msa_len, sizeof(_locr), 1, f);
  fwrite(msa_4bytes, sizeof(unsigned), msa_len, f);
  fwrite(msa_1byte, sizeof(uint8_t), msa_len, f);
  fwrite(table, sizeof(unsigned), REF_TABLE_SIZE, f);
  fwrite(count, sizeof(uint16_t), REF_TABLE_SIZE, f);

  fwrite(&len, sizeof(_locr), 1, f);
  fwrite(ref_sal_4bytes, sizeof(unsigned), len, f);
  fwrite(ref_sal_1byte, sizeof(uint8_t), len, f);
  fwrite(ref_sar, sizeof(unsigned), len, f);
  fwrite(key16, sizeof(uint16_t), len, f);

  free(table);
  free(count);
  free(ref_sal_4bytes);
  free(ref_sal_1byte);
  free(ref_sar);
  free(key16);

  fclose(f);
}

template<typename _locr>
void store_to_file(unsigned id,
		   unsigned *msa,
		   _locr part_len,
		   unsigned msa_len,
		   vector<pair<uint64_t, pair<_locr, unsigned> > > &pattern_ids64,
		   const char *filename){
  FILE *f;
  if (id == 0){
    f = fopen(filename, "wb");
  }
  else {
    f = fopen(filename, "ab");
  }
  
  task_timer timer ("Building the hash table", true);

  unsigned i;
  unsigned len = pattern_ids64.size();
  unsigned hashValue;
  
  unsigned *table = (unsigned*)malloc(REF_TABLE_SIZE * sizeof(unsigned)); // 24 bits
  uint16_t *count = (uint16_t*)malloc(REF_TABLE_SIZE * sizeof(uint16_t));
  memset(table, 0, REF_TABLE_SIZE * sizeof(unsigned));
  memset(count, 0, REF_TABLE_SIZE * sizeof(uint16_t));
  unsigned *ref_sal = (unsigned*)malloc(len * sizeof(unsigned));
  unsigned *ref_sar = (unsigned*)malloc(len * sizeof(unsigned));
  uint16_t *key16 = (uint16_t*)malloc(len * sizeof(uint16_t));
    
  // hash table
  hashValue = (pattern_ids64[0].first >> HT_BITS) & 0xffffff;
  table[hashValue] = 0;
  count[hashValue] = 1;
  ref_sal[0] = pattern_ids64[0].second.first;
  ref_sar[0] = pattern_ids64[0].second.second;
  key16[0] = (uint16_t)(pattern_ids64[0].first & 0xffff);
  
  for (i = 1; i < len; i++){
    if (((pattern_ids64[i].first >> HT_BITS) & 0xffffff) != hashValue){
      hashValue = (pattern_ids64[i].first >> HT_BITS) & 0xffffff;
      table[hashValue] = i;
      count[hashValue] = 1;
    }
    else {
      count[hashValue]++;
    }
    
    ref_sal[i] = pattern_ids64[i].second.first;
    ref_sar[i] = pattern_ids64[i].second.second;
    key16[i] = (uint16_t)(pattern_ids64[i].first & 0xffff);
  }

  /*
  cout << "stat:" << endl;
  unsigned cnt[9] = {0};
  for (i = 0; i < len; i++){
    if (ref_sar[i] <= 128) cnt[0]++;
    else if (ref_sar[i] <= 256) cnt[1]++;
    else if (ref_sar[i] <= 512) cnt[2]++;
    else if (ref_sar[i] <= 1024) cnt[3]++;
    else if (ref_sar[i] <= 2048) cnt[4]++;
    else if (ref_sar[i] <= 4096) cnt[5]++;
    else if (ref_sar[i] <= 8192) cnt[6]++;
    else if (ref_sar[i] <= 65535) cnt[7]++;
    else cnt[8]++;
  }
  for (i = 0; i < 9; i++){
    cout << cnt[i] << endl;
  }
  */

  timer.go("Writing the index to disk");

  if (program_options::aligner_mode != program_options::sensitive){
    fwrite(&part_len, sizeof(_locr), 1, f);
  }

  fwrite(&msa_len, sizeof(unsigned), 1, f);
  fwrite(msa, sizeof(unsigned), msa_len, f);

  fwrite(table, sizeof(unsigned), REF_TABLE_SIZE, f);
  fwrite(count, sizeof(uint16_t), REF_TABLE_SIZE, f);
  
  fwrite(&len, sizeof(unsigned), 1, f);
  fwrite(ref_sal, sizeof(unsigned), len, f);
  fwrite(ref_sar, sizeof(unsigned), len, f);
  fwrite(key16, sizeof(uint16_t), len, f);
  
  free(table);
  free(count);
  
  free(ref_sal);
  free(ref_sar);
  free(key16);
  
  fclose(f);
}

template<typename _val>
void make_db()
{
	verbose_stream << "Database file = " << program_options::input_ref_file << endl;
	
	cpu_timer total;
	task_timer timer ("Opening the database file", true);
	Input_stream db_file (program_options::input_ref_file, true);
	timer.finish();
	
	if ((unsigned)program_options::chunk_size <= 4){
	  if (program_options::aligner_mode == program_options::sensitive)
	    ref_header.block_size = program_options::chunk_size;
	  else
	    ref_header.block_size = (unsigned)program_options::chunk_size * N_PART;
	}
	else
	  ref_header.block_size = program_options::chunk_size;
	
#ifdef EXTRA
	ref_header.sequence_type = sequence_type(_val ());
#endif
	size_t chunk = 0;
	Output_stream main(program_options::database);
	main.write(&ref_header, 1);

	if (program_options::aligner_mode == program_options::sensitive){ // sensitive mode
	  //cout << "sensitive mode" << endl;
	  for(;;++chunk) {
	    timer.go("Loading sequences");
            
	    Sequence_set<Nucleotide>* ss;
            size_t n_seq = load_seqs<_val,_val,Single_strand>(db_file, FASTA_format<_val> (), ref_seqs<_val>::data_, ref_ids::data_, ss, (size_t)(program_options::chunk_size * 1e9));
            log_stream << "load_seqs n=" << n_seq << endl;
	    if(n_seq == 0)
              break;
	    ref_header.letters += ref_seqs<_val>::data_->letters();
            ref_header.sequences += n_seq;
	    
            string file_prefix = program_options::database.substr(0, program_options::database.size()-5) + "_" + to_string(chunk) + ".index"; // cut out ".data"                                              
	    
            const bool long_addressing = (ref_seqs<_val>::data_->letters() + n_seq) > (size_t)std::numeric_limits<uint32_t>::max();
            ref_header.long_addressing = ref_header.long_addressing == true ? true : long_addressing;
	    
            timer.finish();
            ref_seqs<_val>::data_->print_stats();
	    
	    if (long_addressing){ // block size > 4G amino aicids
	      uint64_t len = ref_seqs<_val>::data_->letters() + n_seq;
              char *seqs = new char[len];
              memset(seqs, 'A', sizeof(len));
              uint64_t i, j, k, pattern_id64, msa_len, last, begin, p, id;
              unsigned r;
              uint8_t *pat = new uint8_t[8];

              vector<pair<uint64_t, pair<unsigned, uint8_t> > > pattern_ids64;
              vector<pair<uint64_t, pair<pair<unsigned, uint8_t>, unsigned> > > pattern_ids64_final;

              k = 0;
              memset(seqs, 'A', sizeof(len));
              timer.go("Converting sequences");
	      const Sequence_set<_val> *ref = ref_seqs<_val>::data_;
              for (i = 0; i < n_seq; i++){
                const sequence<const _val> seq = (*ref)[i];

                for (j = 0; j < seq.length(); j++){
                  seqs[k] = Reduction<_val>::reduction(seq[j]) + 'B';
                  k++;
                }
                seqs[k] = 'A';
                k++;
              }

              timer.go("Finding patterns");
              pattern_ids64.clear();
              pattern_ids64.reserve(len);

              i = 0;
	      while (i < len - MAX_SEED_LEN){
                id = 0;

                for (j = 0; j < MAX_SEED_LEN; j++){
                  if (seqs[i+j] == 'A' || seqs[i+j] == 'M'){ // alphabet11
                    //if (seqs[i+j] == 'A' || seqs[i+j] == 'O'){ // alphabet13
                    i = i + j + 1;
                    id = 1;
                    break;
                  }
                }

                if (id == 1) continue;

                // the pattern is valid
                memset(pat, 0, 8);
                pattern_id64 = 0;
                for (j = 0; j < MAX_SEED_LEN; j += 2){
                  pat[j/2] = ((seqs[i+j+1] - 'A') << 4) | (seqs[i+j] - 'A');
                }
		pattern_id64 = (uint64_t)*((uint64_t*)pat);
                pattern_ids64.push_back(make_pair(pattern_id64, make_pair((i >> 8) & 0xffffffff, i & 0xff)));
                i++;
              }
              cout << "pattern_len = 16, # of patterns: " << pattern_ids64.size() << endl;
	      
              msa_len = pattern_ids64.size();
              unsigned *msa_4bytes = new unsigned[msa_len];
              uint8_t *msa_1byte = new uint8_t[msa_len];
	      
	      for (unsigned part = 1; part < N_PATTERN; part++){
		cout << "------ Pattern " << part << " ------" << endl;
		timer.go("Sorting patterns");
		
		std::sort(pattern_ids64.begin(), pattern_ids64.end(), Compare_long(MIS_PATTERN[part]));
		
		timer.go("Constructing the modified suffix array");
		pattern_ids64_final.clear();
		memset(msa_4bytes, 0, msa_len * sizeof(unsigned));
		memset(msa_1byte, 0, msa_len * sizeof(uint8_t));
		
		p = (pattern_ids64[0].first & MIS_PATTERN[0]) >> BS_BITS;
		begin = 0;
		last = 0;
		
		for (i = 0; i < pattern_ids64.size(); i++){
		  if (((pattern_ids64[i].first & MIS_PATTERN[0]) >> BS_BITS) != p){
		    r = (unsigned)(last - begin);
		    pattern_ids64_final.push_back(make_pair(p, make_pair(make_pair((begin >> 8) & 0xffffffff, begin & 0xff), r)));
		    
		    begin = last;
		    p = (pattern_ids64[i].first & MIS_PATTERN[0]) >> BS_BITS;
		  }
		  
		  uint64_t loc = (((uint64_t)pattern_ids64[i].second.first << 8) | pattern_ids64[i].second.second) + String_set<_val>::PERIMETER_PADDING;
		  msa_4bytes[last] = (loc >> 8) & 0xffffffff;
		  msa_1byte[last] = loc & 0xff;
		  last++;
		}
		
		cout << "# of unique patterns(40 bits): " << pattern_ids64_final.size() << endl;
		
		timer.go("Saving the index to disk");
		store_to_file_long<uint64_t>(part-1, msa_4bytes, msa_1byte, msa_len, pattern_ids64_final, file_prefix.c_str());
	      }
	      
	      delete []msa_4bytes;
	      delete []msa_1byte;
		
	      delete []seqs;
	      pattern_ids64.clear();
	      pattern_ids64_final.clear();
	      
	      timer.finish();
	      timer.go("Saving the references to disk");
	      ref_seqs<_val>::data_->save(main);
	      ref_ids::get().save(main);
	      
	      timer.go("Deallocating sequences");
              delete ref_seqs<_val>::data_;
              delete ref_ids::data_;
              delete ss;
	    }
	    else {
	      uint64_t len = ref_seqs<_val>::data_->letters() + n_seq;
	      char *seqs = new char[len];
	      memset(seqs, 'A', sizeof(len));
	      uint64_t i, j, k;
	      unsigned r;
	      uint64_t pattern_id64;
	      unsigned id;
	      uint8_t *pat = new uint8_t[8];
	      unsigned msa_len;
	      unsigned last, begin;
	      uint64_t p;
	      
	      vector<pair<uint64_t, unsigned> > pattern_ids64;
	      vector<pair<uint64_t, pair<unsigned, unsigned> > > pattern_ids64_final;
	      
	      k = 0;
	      timer.go("Converting sequences");
	      const Sequence_set<_val> *ref = ref_seqs<_val>::data_;
	      for (i = 0; i < n_seq; i++){
		const sequence<const _val> seq = (*ref)[i];
		
		for (j = 0; j < seq.length(); j++){
		  seqs[k] = Reduction<_val>::reduction(seq[j]) + 'B';
		  k++;
		}
		seqs[k] = 'A';
		k++;
	      }
	      
	      timer.go("Finding patterns");
	      pattern_ids64.clear();
	      pattern_ids64.reserve(len);
	      
	      i = 0;
	      while (i < len - MAX_SEED_LEN){
		id = 0;
		
		for (j = 0; j < MAX_SEED_LEN; j++){
		  if (seqs[i+j] == 'A' || seqs[i+j] == 'M'){ // alphabet11
		    //if (seqs[i+j] == 'A' || seqs[i+j] == 'O'){ // alphabet13
		    i = i + j + 1;
		    id = 1;
		    break;
		  }
		}
		
		if (id == 1) continue;
		
		memset(pat, 0, 8);
		pattern_id64 = 0;
		for (j = 0; j < MAX_SEED_LEN; j += 2){
		  pat[j/2] = ((seqs[i+j+1] - 'A') << 4) | (seqs[i+j] - 'A');
		}
		pattern_id64 = (uint64_t)*((uint64_t*)pat);
		pattern_ids64.push_back(make_pair(pattern_id64, i));
		i++;
	      }
	      cout << "pattern_len = 16, # of patterns: " << pattern_ids64.size() << endl;
	      msa_len = pattern_ids64.size();
	      unsigned *msa = new unsigned[msa_len];
	      
	      for (unsigned part = 1; part < N_PATTERN; part++){
		cout << "------ Pattern " << part << " ------" << endl;
		timer.go("Sorting patterns");
		std::sort(pattern_ids64.begin(), pattern_ids64.end(), Compare<unsigned>(MIS_PATTERN[part]));
		
		timer.go("Constructing the modified suffix array");
		pattern_ids64_final.clear();
		memset(msa, 0, msa_len * sizeof(unsigned));
		
		p = (pattern_ids64[0].first & MIS_PATTERN[part]) >> BS_BITS;
		begin = 0;
		last = 0;
		
		for (i = 0; i < pattern_ids64.size(); i++){
		  if (((pattern_ids64[i].first & MIS_PATTERN[part]) >> BS_BITS) != p){
		    r = last - begin;
		    pattern_ids64_final.push_back(make_pair(p, make_pair(begin, r)));
		    
		    begin = last;
		    p = (pattern_ids64[i].first & MIS_PATTERN[part]) >> BS_BITS;
		  }
		  
		  msa[last] = pattern_ids64[i].second + String_set<_val>::PERIMETER_PADDING;
		  last++;
		}
		
		cout << "# of unique patterns(40 bits): " << pattern_ids64_final.size() << endl;
		
		timer.go("Saving the index to disk");
		store_to_file<unsigned>(part-1, msa, 0, msa_len, pattern_ids64_final, file_prefix.c_str());
	      }
	      
	      delete []msa;
	      delete []seqs;
	      pattern_ids64.clear();
	      pattern_ids64_final.clear();
	      
	      timer.finish();
	      timer.go("Saving the references to disk");
	      ref_seqs<_val>::data_->save(main);
	      ref_ids::get().save(main);
	      
	      timer.go("Deallocating sequences");
	      delete ref_seqs<_val>::data_;
	      delete ref_ids::data_;
	      delete ss;
	    }
	  }
	}
	else { // fast mode
	  for(;;++chunk) {
	    timer.go("Loading sequences");
	    
	    const bool long_addressing = (unsigned)program_options::chunk_size > 4;
	    unsigned chunkSize;
	    if (long_addressing) chunkSize = (unsigned)program_options::chunk_size;
	    else chunkSize = (unsigned)program_options::chunk_size * N_PART;
	    
	    Sequence_set<Nucleotide>* ss;
	    size_t n_seq = load_seqs<_val,_val,Single_strand>(db_file, FASTA_format<_val> (), ref_seqs<_val>::data_, ref_ids::data_, ss, (size_t)(chunkSize * 1e9));
	    log_stream << "load_seqs n=" << n_seq << endl;
	    if(n_seq == 0)
	      break;
	    ref_header.letters += ref_seqs<_val>::data_->letters();
	    ref_header.sequences += n_seq;
	    
	    string file_prefix = program_options::database.substr(0,program_options::database.size()-5) + "_" + to_string(chunk) + ".index"; // cut out ".data"
	    
	    //const bool long_addressing = (ref_seqs<_val>::data_->letters() + n_seq) > (size_t)std::numeric_limits<uint32_t>::max();
	    ref_header.long_addressing = ref_header.long_addressing == true ? true : long_addressing;
	    
	    timer.finish();
	    ref_seqs<_val>::data_->print_stats();
	    
	    if (long_addressing){ // block size > 4G amino acids
	      //cout << "long_addressing: " << (unsigned)ref_header.long_addressing << endl;
	      uint64_t len = ref_seqs<_val>::data_->letters() + n_seq;
	      char *seqs = new char[len];
	      memset(seqs, 'A', sizeof(len));
	      uint64_t i, j, k, pattern_id64, msa_len, last, begin, p, id;
	      unsigned r;
	      uint8_t *pat = new uint8_t[8];
	      
	      vector<pair<uint64_t, pair<unsigned, uint8_t> > > pattern_ids64;
	      vector<pair<uint64_t, pair<pair<unsigned, uint8_t>, unsigned> > > pattern_ids64_final;
	      
	      k = 0;
	      memset(seqs, 'A', sizeof(len));
	      timer.go("Converting sequences");
	      const Sequence_set<_val> *ref = ref_seqs<_val>::data_;
	      for (i = 0; i < n_seq; i++){
		const sequence<const _val> seq = (*ref)[i];
	
		for (j = 0; j < seq.length(); j++){
		  seqs[k] = Reduction<_val>::reduction(seq[j]) + 'B';
		  k++;
		}
		seqs[k] = 'A';
		k++;
	      }
	      
	      timer.go("Finding patterns");
	      pattern_ids64.clear();
	      pattern_ids64.reserve(len);
	      
	      i = 0;
	      while (i < len - MAX_SEED_LEN){
		id = 0;
		
		for (j = 0; j < MAX_SEED_LEN; j++){
		  if (seqs[i+j] == 'A' || seqs[i+j] == 'M'){ // alphabet11
		    //if (seqs[i+j] == 'A' || seqs[i+j] == 'O'){ // alphabet13
		    i = i + j + 1;
		    id = 1;
		    break;
		  }
		}
		
		if (id == 1) continue;
		
		// the pattern is valid
		memset(pat, 0, 8);
		pattern_id64 = 0;
		for (j = 0; j < MAX_SEED_LEN; j += 2){
		  pat[j/2] = ((seqs[i+j+1] - 'A') << 4) | (seqs[i+j] - 'A');
		}
		pattern_id64 = (uint64_t)*((uint64_t*)pat);
		pattern_ids64.push_back(make_pair(pattern_id64, make_pair((i >> 8) & 0xffffffff, i & 0xff)));
		i++;
	      }
	      cout << "pattern_len = 16, # of patterns: " << pattern_ids64.size() << endl;
	      
	      msa_len = pattern_ids64.size();
	      unsigned *msa_4bytes = new unsigned[msa_len];
	      uint8_t *msa_1byte = new uint8_t[msa_len];

	      timer.go("Sorting patterns");
	      std::sort(pattern_ids64.begin(), pattern_ids64.end(), Compare_long(MIS_PATTERN[0]));
	      
	      timer.go("Constructing the modified suffix array");
	      pattern_ids64_final.clear();
	      memset(msa_4bytes, 0, msa_len * sizeof(unsigned));
	      memset(msa_1byte, 0, msa_len * sizeof(uint8_t));

	      p = (pattern_ids64[0].first & MIS_PATTERN[0]) >> BS_BITS;
	      begin = 0;
	      last = 0;
	      
	      for (i = 0; i < pattern_ids64.size(); i++){
		if (((pattern_ids64[i].first & MIS_PATTERN[0]) >> BS_BITS) != p){
		  r = (unsigned)(last - begin);
		  pattern_ids64_final.push_back(make_pair(p, make_pair(make_pair((begin >> 8) & 0xffffffff, begin & 0xff), r)));
		  
		  begin = last;
		  p = (pattern_ids64[i].first & MIS_PATTERN[0]) >> BS_BITS;
		}
		
		uint64_t loc = (((uint64_t)pattern_ids64[i].second.first << 8) | pattern_ids64[i].second.second) + String_set<_val>::PERIMETER_PADDING;
		msa_4bytes[last] = (loc >> 8) & 0xffffffff;
		msa_1byte[last] = loc & 0xff;
		last++;
	      }
	      
	      cout << "# of unique patterns(40 bits): " << pattern_ids64_final.size() << endl;
	      
	      timer.go("Saving the index to disk");
	      store_to_file_long<uint64_t>(0, msa_4bytes, msa_1byte, msa_len, pattern_ids64_final, file_prefix.c_str());
	      
	      delete []msa_4bytes;
	      delete []msa_1byte;
	      
	      delete []seqs;
	      pattern_ids64.clear();
	      pattern_ids64_final.clear();
	      
	      timer.finish();
	      timer.go("Saving the references to disk");
	      ref_seqs<_val>::data_->save(main);
	      ref_ids::get().save(main);
	      
	      timer.go("Deallocating sequences");
	      delete ref_seqs<_val>::data_;
	      delete ref_ids::data_;
	      delete ss;
	    }
	    else {
	      uint64_t len = ref_seqs<_val>::data_->letters() + n_seq;
	      char *seqs = new char[len];
	      memset(seqs, 'A', sizeof(len));
	      uint64_t i, j, k, part_id, part_len;
	      unsigned r;
	      part_id = 0;
	      uint64_t pattern_id64;
	      unsigned id;
	      uint8_t *pat = new uint8_t[8];
	      unsigned msa_len;
	      //unsigned *msa = new unsigned[msa_len];
	      unsigned last, begin;
	      uint64_t p;

	      vector<pair<uint64_t, unsigned> > pattern_ids64;
	      vector<pair<uint64_t, pair<unsigned, unsigned> > > pattern_ids64_final;
	      
	      for (unsigned part = 0; part < N_PART; part++){
		k = 0;
		memset(seqs, 'A', sizeof(len));
		timer.go("Converting sequences");
		const Sequence_set<_val> *ref = ref_seqs<_val>::data_;
		for (i = part_id; i < n_seq; i++){
		  const sequence<const _val> seq = (*ref)[i];

		  if (k > len / N_PART) break;

		  for (j = 0; j < seq.length(); j++){
		    seqs[k] = Reduction<_val>::reduction(seq[j]) + 'B';
		    k++;
		  }
		  seqs[k] = 'A';
		  k++;
		}
		part_id = i;
		part_len = k;

		timer.go("Finding patterns");
		pattern_ids64.clear();
		pattern_ids64.reserve(part_len);

		i = 0;
		while (i < part_len - MAX_SEED_LEN){
		  id = 0;

		  for (j = 0; j < MAX_SEED_LEN; j++){
		    if (seqs[i+j] == 'A' || seqs[i+j] == 'M'){ // alphabet11
		      //if (seqs[i+j] == 'A' || seqs[i+j] == 'O'){ // alphabet13
		      i = i + j + 1;
		      id = 1;
		      break;
		    }
		  }

		  if (id == 1) continue;

		  // the pattern is valid
		  memset(pat, 0, 8);
		  pattern_id64 = 0;
		  for (j = 0; j < MAX_SEED_LEN; j += 2){
		    pat[j/2] = ((seqs[i+j+1] - 'A') << 4) | (seqs[i+j] - 'A');
		  }
		  pattern_id64 = (uint64_t)*((uint64_t*)pat);
		  pattern_ids64.push_back(make_pair(pattern_id64, i));
		  i++;
		}
		cout << "pattern_len = 16, # of patterns: " << pattern_ids64.size() << endl;

		msa_len = pattern_ids64.size();
		unsigned *msa = new unsigned[msa_len];

		timer.go("Sorting patterns");
		std::sort(pattern_ids64.begin(), pattern_ids64.end(), Compare<unsigned>(MIS_PATTERN[0]));

		timer.go("Constructing the modified suffix array");
		pattern_ids64_final.clear();
		memset(msa, 0, msa_len * sizeof(unsigned));

		p = (pattern_ids64[0].first & MIS_PATTERN[0]) >> BS_BITS;
		begin = 0;
		last = 0;

		for (i = 0; i < pattern_ids64.size(); i++){
		  if (((pattern_ids64[i].first & MIS_PATTERN[0]) >> BS_BITS) != p){
		    r = last - begin;
		    pattern_ids64_final.push_back(make_pair(p, make_pair(begin, r)));
		    
		    begin = last;
		    p = (pattern_ids64[i].first & MIS_PATTERN[0]) >> BS_BITS;
		  }

		  msa[last] = pattern_ids64[i].second + String_set<_val>::PERIMETER_PADDING;
		  last++;
		}

		cout << "# of unique patterns(40 bits): " << pattern_ids64_final.size() << endl;
		
		timer.go("Saving the index to disk");
		store_to_file<unsigned>(part, msa, part_len, msa_len, pattern_ids64_final, file_prefix.c_str());

		delete []msa;
	      }
	      
	      delete []seqs;
	      pattern_ids64.clear();
	      pattern_ids64_final.clear();

	      timer.finish();
	      timer.go("Saving the references to disk");
	      ref_seqs<_val>::data_->save(main);
	      ref_ids::get().save(main);

	      timer.go("Deallocating sequences");
	      delete ref_seqs<_val>::data_;
	      delete ref_ids::data_;
	      delete ss;
	    }	  
	  }  
	}

	timer.finish();
	ref_header.n_blocks = chunk;
	log_stream << "db seek" << endl;
	main.seekp(0);
	log_stream << "db write" << endl;
	main.write(&ref_header, 1);
	log_stream << "db close" << endl;
	main.close();
	verbose_stream << "Total time = " << boost::timer::format(total.elapsed(), 1, "%ws\n");
}

#endif /* MAKE_DB_H_ */
