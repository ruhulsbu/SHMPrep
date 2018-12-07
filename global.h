#include "library.h"
#include "constant.h"
#include "structure.h"

#ifndef GLOBAL_H
#define GLOBAL_H

extern int KMER;// = 11;		//Size of KMER. Choose all KMER from Chromosome and Target and hash them. 
extern int SEED;// = 40;		//For the Set of KMER in due order find the SEED with 55% match
extern int MINREAD;// = 0;		//Start from Read No MINREAD
extern int MAXREAD;// = 100;		//End before MAXREAD is reached
extern int LOCALFLAG;
extern int PRINTDETAIL;
extern int GLOBAL_KBAND;
extern int PAIRED_READ_COUNT;

extern float OVERLAP_QUALITY;
extern float GLOBAL_QUALITY;
extern float PRIMER_QUALITY;

extern int DUPCOUNT;
extern int FASTQOUTPUT;
extern int STRIPOFF;
extern int MAXTHREAD;
extern ifstream fp_read_first;
extern ifstream fp_read_second;

extern int BARCODED_SEQUENCE;
extern int FORWARD_ADAPTER_LENGTH;
extern int REVERSE_ADAPTER_LENGTH;
extern vector<string> primer_vector;
extern unordered_map<char, float> log_error;
extern unordered_map<char, float> log_quality;
extern unordered_map<string, int> map_primers;
extern unordered_map<string, string> primer_name;
//extern unordered_map<string, vector<pair<string, string> > > map_barcode;

extern int t_lookup;// = 0;
extern int t_seed;// = 0;
extern int t_extend;// = 0;
extern int t_chain;// = 0;
extern int t_lis;// = 0;
extern int t_kmer_count;// = 0;
extern long base_power_kmer[21][4];
extern long error_free_seg[1000];
extern long error_dist[10];

extern int matching_score[51][51];
extern int mismatch_score[51][51];

extern vector<reference_index> refindex;	//vector of reference

extern int GAP;
extern int BASIC_CORRECTION;
extern int PERBASE_QSCORE;
extern int GLOBAL_QSCORE;
extern int EXCLUDE_INDELS;
extern int REFERENCELESS_ALIGNMENT;
extern int FILTER_LOW_FREQ_READS;
extern int FILTER_LOW_CONS_READS;
#endif
