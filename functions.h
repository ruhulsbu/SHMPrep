#include "library.h"
#include "constant.h"
#include "structure.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

//align_reads.cpp
extern int read_sequence_file(string& read_name_first, string& readseq_first, string& quality_first,
                                string& read_name_second, string& readseq_second, string& quality_second);
extern string return_primer_using_kmer_match(string &sequence, unordered_map<string, vector<int> >& primer_map);
extern string correct_alignment_direction(string& readseq_first, string& quality_first, meta_data& meta_details, ofstream& fp_detail);
extern void aligner_thread(meta_data& meta_details, vector<reference_index>& refindex);
//extern void align_reads(int& threadno, meta_data& meta_details, vector<reference_index>& refindex);
extern void align_reads(meta_data& meta_details, vector<reference_index>& refindex, int threadno);

//extern void align_reads(string& first_read_file, string& second_read_file, string& sam_file, vector<reference_index>& refindex);
extern void read_vs_reference(string& read, string& read_name, int dir, vector<reference_index>& refindex, 
			vector<pair<int, vector<pair<int, int> > > >& primary_chain);
extern void refine_kmer_index(vector<pair<int, int> >& kmer_index, vector<pair<int, vector<pair<int, int> > > >& primary_chain,
                                string read, int dir, vector<reference_index>& refindex, int index);
extern void create_primary_chain(vector<pair<int, int> >& kmer_index, vector<pair<int, int> >& chain,
                int start_ind, int next_ind, int readlen);
//end

//print_alignment.cpp
extern void prepare_fasta_output_heuristic(ifstream& fp_input, ofstream& fp_fasta, ofstream& fp_detail, string& reference, meta_data& meta_details);
extern void prepare_fasta_output_distribution(ifstream& fp_input, ofstream& fp_fasta, ofstream& fp_detail, string& reference, meta_data& meta_details);
extern void print_vector_alignment(vector<pair<char, char> >& alignment);
extern void print_alignment_back(vector<pair<char, char> >& alignment, int ref_position, int read_position, int step);
extern void print_alignment(vector<pair<char, char> >& alignment, string& ref, string& read, int ref_position, 
		int read_position, int step, fragment_alignment &fragment_alignment_info, bool flag, ofstream& fp_detail);
//end

//edit_distance
extern int find_kband_similarity(string& str1, string& str2, vector<pair<char, char> >& alignment, int &str1_end, int &str2_end, cell **matrix);
extern int find_local_similarity(string& str1, string& str2, vector<pair<char, char> >& alignment, int& str1_start, int& str2_start, 
						int& str1_end, int& str2_end, cell **matrix);
extern int find_banded_similarity(string& str1, string& str2, vector<pair<char, char> >& alignment, int &str1_end, int &str2_end);
extern int find_levenshtein_distance(string& str1, string& str2);
extern int find_longest_common_subsequence(string& str1, string& str2);
//end

//optimize_kband.cpp
extern void merge_both_alignments(vector<reference_index>& refindex, string sam_output_name, fragment_alignment& alignment_first, 
				fragment_alignment& alignment_second, ofstream& fp_sam, ofstream& fp_detail, 
				ofstream& fp_fastq, unordered_map<string, int>& umap);
extern int referenceless_alignment(string& readseq_first, string& readseq_second, string& quality_first, string& quality_second, 
					string& readsequence, string& readquality, cell **matrix, ofstream& fp_detail);
extern int create_gap_alignment(vector<pair<char, char> >& alignment);
extern int optimize_path(vector<pair<char, char> >& alignment);
extern void print_path_matrix(long long **path, int row, int column);
extern void print_path_cell_back(cell **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment);
extern void print_local_alignment_back(cell **path, int row, int column, string& str1, string& str2, 
						int& str1_start, int& str2_start, vector<pair<char, char> >& alignment);
extern void print_path_back(long long **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment);
extern void print_path(long long **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment);
//end

//sam_format.cpp
extern void sam_format(fragment_alignment final_alignment_info, vector<reference_index> &refindex, string& read, 
			string& read_name, vector<string>& output, ofstream& fp_detail);
//end

//utility.cpp
extern char complement(char ch);
extern int map_value(char ch);
extern void init_matrix(cell **matrix);
extern void remove_matrix(cell **matrix);

extern int similarity(char x, char y);
extern int min(int x, int y);
extern int max(int x, int y);

extern void upper_case(string& str);
extern void reverse_str(string &str);
extern string reverse_complement(string& str);
extern char calculate_basequal(char first, char second, int flag);

extern bool validate_alignment(vector<reference_index>& refindex, int ref_ind, vector<pair<char, char> >& alignment, 
				int ref_start, int read_start, string& read);
extern bool validate_chain(vector<reference_index>& refindex, int ref_ind, fragment_alignment& alignment_info, string& read);

extern void preprocess_thread(meta_data& meta_details, vector<reference_index>& refindex);
//end

//find_lcs.cpp
extern void create_automata(string& reference, automata_state automata[]);
extern string find_primer(automata_state automata[], string& query, int& refindex, int& queryindex, ofstream &fp_detail);

extern int find_lcs(string& reference, string& query, int& refindex, int& queryindex, int& max_len, 
			automata_state automata[], int& sz, int& last, bool extend, ofstream& fp_detail);
//end

#endif
