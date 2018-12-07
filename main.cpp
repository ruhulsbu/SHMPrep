#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

//04-21-2015
int KMER = 8;			//Size of KMER. Find KMERs from Chromosome and hash them. 
int SEED = 40;			//SEEDs are not used in this program
int MINREAD = 0;		//Start from Read No MINREAD (used for Debugging)
int MAXREAD = 0;		//End before MAXREAD is reached (used for Debugging)
int MAXTHREAD = 1;		//Number of Threads 
int DUPCOUNT = 0;
int FASTQOUTPUT = 0;
int LOCALFLAG = 0;
int PRINTDETAIL = 0;
int PAIRED_READ_COUNT = 0;
int STRIPOFF = STRIP_OFF;
int GLOBAL_KBAND = KBAND;

int GAP = GAPX;			//Penalty score for GAP
int BASIC_CORRECTION = 0;	//Flag for Basic Correction 
int PERBASE_QSCORE = 30;	//Per base quality score
int GLOBAL_QSCORE = 25;		//Global quality score
int EXCLUDE_INDELS = 0;		//Flag for excluding indels
int FILTER_LOW_FREQ_READS = 0;	//Flag for filtering low quality reads
int FILTER_LOW_CONS_READS = 0;
int REFERENCELESS_ALIGNMENT = 0;//Flag for Referenceless Alignment

long base_power_kmer[21][4];
long error_free_seg[1000];
long error_dist[10];

int matching_score[51][51];	//Precalculated probability for matching bases
int mismatch_score[51][51];	//Precalculated probability for mismatching bases

float OVERLAP_QUALITY = OVERLAP_PERCENT_MATCH;
float GLOBAL_QUALITY = GLOBAL_PERCENT_MATCH;
float PRIMER_QUALITY = DEFAULT_PRIMER_QUALITY;

bool JOBDONE = false;
bool ALTERNATIVE_RUN = false;

string meta_file = "";// "./meta.txt";//met file path
ifstream fp_meta;

int BARCODED_SEQUENCE = 0;	//Flag for aligning Barcoded sequences
int FORWARD_ADAPTER_LENGTH = 0;	//Forward Adapter length
int REVERSE_ADAPTER_LENGTH = 0;	//Reverse Adapter length
vector<string> primer_vector;		
unordered_map<char, float> log_error;	//Precalculated map for log of error
unordered_map<char, float> log_quality;	//Precalculated map for log of quality
unordered_map<string, int> map_primers;	//Map for storing KMER to Primer List
unordered_map<string, string> primer_name;
//unordered_map<string, vector<pair<string, string> > > map_barcode; //Map to store set of alignments for a Barcode

std::mutex meta_mutex;

void input_reference(vector<pair<string, string> >& reference, string& ref_file, string& sam_file);
void kmer_inverted_index(vector<pair<string, string> >& reference, vector<reference_index>& refindex);


void init_variables()
{
	/*
	KMER = 8;			//Size of KMER. Find KMERs from Chromosome and hash them. 
	SEED = 40;			//SEEDs are not used in this program
	MINREAD = 0;			//Start from Read No MINREAD (used for Debugging)
	MAXREAD = 0;			//End before MAXREAD is reached (used for Debugging)
	MAXTHREAD = 1;			//Number of Threads 
	FASTQOUTPUT = 0;

	GAP = GAPX;			//Penalty score for GAP
	BASIC_CORRECTION = 0;		//Flag for Basic Correction 
	PERBASE_QSCORE = 30;		//Per base quality score
	GLOBAL_QSCORE = 25;		//Global quality score
	EXCLUDE_INDELS = 0;		//Flag for excluding indels
	FILTER_LOW_FREQ_READS = 0;	//Flag for filtering low quality reads
	REFERENCELESS_ALIGNMENT = 0;	//Flag for Referenceless Alignment

	BARCODED_SEQUENCE = 0;		//Flag for aligning Barcoded sequences
	FORWARD_ADAPTER_LENGTH = 0;	//Forward Adapter length
	REVERSE_ADAPTER_LENGTH = 0;	//Reverse Adapter length
	*/
	primer_vector.clear();		
	map_primers.clear();		//Map for storing KMER to Primer List
	//map_barcode.clear(); 		//Map to store set of alignments for a Barcode

}


//Create an inverted index for KMER to index of KMER in the reference
void kmer_inverted_index(vector<pair<string, string> >& reference, vector<reference_index>& refindex)
{
	int input, cases = 0, overlap, index_time;

	string str1, str2;
	for(int i = 0; i < KMER; i++)
		for(int k = 0; k < BASE; k++ )
			base_power_kmer[i][k] = pow(BASE, i) * k;

	long kmer_choices = (long) (pow(4, KMER) + 9);	
	for(int i = 0; i < reference.size(); i++)
	{
		reference_index refinfo;
		refinfo.name = reference.at(i).first;
		refinfo.ref = reference.at(i).second;
		refinfo.rev = reverse_complement(reference.at(i).second);

		refinfo.index = new int [kmer_choices];
		memset(refinfo.index, -1, (sizeof(int) * kmer_choices));
		refinfo.revind = new int [kmer_choices];
		memset(refinfo.revind, -1, (sizeof(int) * kmer_choices));

		//cout << "ref: " << refinfo.name << " and name = " << refinfo.ref.size() << endl;

		long hash_key = 0;
		int map_val = 0;
		for(int k = 0, i = KMER - 2; k < KMER - 1; k++, i--)
		{
			map_val = map_value(refinfo.ref.at(k));
			hash_key += base_power_kmer[i][map_val];
		}

		for(int k = KMER - 1; k < refinfo.ref.length(); k++)
		{
			map_val = map_value(refinfo.ref.at(k));
			hash_key = (hash_key << 2) + map_val;
			
			if(refinfo.index[hash_key] == -1)
			{
				vector<int> pos;
				refinfo.position.push_back(pos);
				refinfo.index[hash_key] = refinfo.position.size() - 1;
			} 
					
			refinfo.position[refinfo.index[hash_key]].push_back(k - KMER + 1);		
			map_val = map_value(refinfo.ref.at(k - KMER + 1));
			hash_key -= base_power_kmer[KMER - 1][map_val];
		}

		/////////////////////////////////////////////////////////////////////////////////////

		hash_key = 0;
		map_val = 0;
		for(int k = 0, i = KMER - 2; k < KMER - 1; k++, i--)
		{
			map_val = map_value(refinfo.rev.at(k));
			hash_key += base_power_kmer[i][map_val];
		}

		for(int k = KMER - 1; k < refinfo.rev.length(); k++)
		{
			map_val = map_value(refinfo.rev.at(k));
			hash_key = (hash_key << 2) + map_val;
			
			if(refinfo.revind[hash_key] == -1)
			{
				vector<int> pos;
				refinfo.position.push_back(pos);
				refinfo.revind[hash_key] = refinfo.position.size() - 1;
			} 
					
			refinfo.position[refinfo.revind[hash_key]].push_back(k - KMER + 1);
			map_val = map_value(refinfo.rev.at(k - KMER + 1));
			hash_key -= base_power_kmer[KMER - 1][map_val];
		}
	
		//cout << "Total index size = " << refinfo.position.size() << endl;	
		refindex.push_back(refinfo);
		
	}

}


//Create a Map for KMER to Primer List
void create_primer_map(string& sequence, unordered_map<string, vector<int> >& primer_map)
{
	string primer_kmer;
	int primer_id;

	for(int i = 0; i < sequence.length() - PRIMER_KMER_SIZE + 1; i++)
	{
		primer_kmer = sequence.substr(i, PRIMER_KMER_SIZE);
		if(primer_map.find(primer_kmer) == primer_map.end())
		{
			vector<int> primer_list;
			primer_map[primer_kmer] = primer_list;
		}

		primer_id = map_primers[sequence] / PRIMER_ID_CONSTANT;
		primer_map[primer_kmer].push_back(primer_id);
	}

}


void generate_iupac_primers(int index, string& sequence, string primer, vector<string>& primer_list)
{
	if(primer.length() == sequence.length())
	{
		//cout << "Primer: " << primer << endl;
		primer_list.push_back(primer);
		return;
	}
	switch(sequence[index])
	{
		case 'A':
		case 'C':
		case 'G':
		case 'T':
			primer += sequence[index];
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;
		case 'M':
			primer += 'A';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'C';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;
		case 'R':
			primer += 'A';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'G';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;
		case 'W':
			primer += 'A';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'T';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;	
		case 'S':
			primer += 'C';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'G';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;	
		case 'Y':
			primer += 'C';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'T';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;	
		case 'K':
			primer += 'G';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'T';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;	
		case 'V':
			primer += 'A';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'C';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'G';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;	
		case 'H':
			primer += 'A';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'C';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'T';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;	
		case 'D':
			primer += 'A';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'G';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'T';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;	
		case 'B':
			primer += 'C';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'G';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'T';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;	
		case 'N':
			primer += 'A';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'C';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'G';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			primer[index] = 'T';
			generate_iupac_primers(index + 1, sequence, primer, primer_list);
			break;	
				
	}
}

//Parse the Primers from the given primer file
void concat_primers(string& primer_file, string& primer_string, int adapter,
			unordered_map<string, vector<int> >& forward_primer_map, 
			unordered_map<string, vector<int> >& complement_fprimer_map,
			int& max_len, int& min_barcode)
{
	adapter = 0;	//disable adapters: FORWARD_ADAPTER_LENGTH, REVERSE_ADAPTER_LENGTH
	//int max_len = -1, find_n;
	int forward_primer_id, complement_fprimer_id, find_n;
	primer_vector.push_back("");

	ifstream fp_primer;
	fp_primer.open(primer_file.c_str(), ifstream::in);
	//cout << "Reading Primer File = " << primer_file << endl;

	string input = "";
	string line = "";
	string nchars = "";
	string sequence, revseq;
	primer_string = "";

	max_len = -1;
	min_barcode = 999;
	
	while(getline(fp_primer, input))
	{
		if(input.length() <= 0)
			break;

		cout << "Primer Name: " << input << ", Check!" << endl;
		while(!input.empty() && isspace(input.back())) 
			input.pop_back() ;
		cout << "Removed Trailing: " << input << ", Spaces!" << endl;

		getline(fp_primer, line);
		upper_case(line);
		find_n = line.rfind("N") + 1;//Use # or X instead of N
		sequence = line.substr(find_n, line.length() - find_n);
		revseq = reverse_complement(sequence);
		//cout << "Sequence = " << sequence << endl;
		if(find_n >= 1)
			nchars = line.substr(adapter, find_n - adapter);
		//cout << "NChars = " << nchars << endl;
		assert(adapter + sequence.length() + nchars.length() == line.length());	

		//forward_primer_id += 1;
		//complement_fprimer_id += 1;

		string primer = "";
		vector<string> primer_list;
		generate_iupac_primers(0, sequence, primer, primer_list);
		for(int i = 0; i < primer_list.size(); i++)
		{
			sequence = primer_list[i];
			revseq = reverse_complement(sequence);
			primer_string += ('$' + sequence);

			if(map_primers.find(sequence) != map_primers.end() || map_primers.find(revseq) != map_primers.end())
			{
				cout << "Conflict in Primer Sequence Found!!!" << endl;
				exit(0);
			}

			forward_primer_id = primer_vector.size();
			primer_vector.push_back(sequence);
			primer_name[sequence] = input.substr(1, input.length() - 1);
			map_primers[sequence] = forward_primer_id * PRIMER_ID_CONSTANT + nchars.length();
			create_primer_map(sequence, forward_primer_map);
			//cout << "Primer Map Size = " << forward_primer_map.size() << endl;

			complement_fprimer_id = primer_vector.size();
			primer_vector.push_back(revseq);
			primer_name[revseq] = input.substr(1, input.length() - 1);
			map_primers[revseq] = complement_fprimer_id * PRIMER_ID_CONSTANT + nchars.length();
			create_primer_map(revseq, complement_fprimer_map);
			//cout << "Complement Primer Map Size = " << complement_fprimer_map.size() << endl;

			max_len = max(max_len, sequence.length() + nchars.length());
		}

		min_barcode = min(nchars.length(), min_barcode);
	}

	if(DEBUG == 99)
		cout << endl;
	primer_string += '$';
	//upper_case(primer_string);
	fp_primer.close();
	
	//exit(0);//for testing this function
	return;// max_len;
}


//Parse the reference from the given Consensus file
void input_reference(vector<pair<string, string> >& reference, string& ref_file, string& sam_file)
{
	ifstream fp_ref;
	ofstream fp_sam;
	
	char *ref = new char[ref_file.length() + 1];
	strcpy(ref, ref_file.c_str());

	char *sam = new char[sam_file.length() + 1];
	strcpy(sam, sam_file.c_str());
	

	fp_ref.open(ref, ifstream::in);
	fp_sam.open(sam, ofstream::out);

	string space = " ";
	string input, output;

	string ref_name;
	string sequence;

	getline(fp_ref, input);
	
	while(!fp_ref.eof())
	{
		//getline(fp_ref, input);
		size_t find = input.find(space);
		if(find != string::npos)
			ref_name = input.substr(1, find);
		else
			ref_name = input.substr(1);
		sequence = "";		
		if(DEBUG == 99)
			cout << ref_name << endl;

		while(getline(fp_ref, input))
		{
			if(input.length() < 1)
				break;

			if(input.at(0) == '>')
			{
				break;
			}
			sequence += input;
		}

		//if(ref_name.find("chr13") == string::npos)
		//	continue;//added on 03-13-15
		
		upper_case(sequence);
		if(DEBUG == 99)
			cout << sequence << endl;		
		reference.push_back(make_pair(ref_name, sequence));
		output = "@SQ\tSN:" + ref_name + "\tLN:";
		fp_sam << output << sequence.length() << endl;

		if(input.length() < 1)
			break; //added on 03-11-15		
	}

	fp_ref.close();
	fp_sam.close();

	delete [] ref;
	delete [] sam;
}


//Parse the Job information from the given Meta file
void read_meta_file(meta_data& meta_details, int threadno, int jobno)
{
	
	meta_mutex.lock();
	
	string line = "", data[7];
	int i = 0;
	
	getline(fp_meta, line);
	if(line.length() == 0)
	{
		//cout << "Returning from Thread " << (threadno + 1) << " with Job Done " << jobno << endl;
		cout << "Total Number of Jobs Completed: " << jobno << endl;
		meta_details.index = "";
		meta_mutex.unlock();
		return;
	}
	else if(line.length() > 0 && line[0] == '#')
	{
		cout << "Skipping Meta file input starting with #" << endl;
		meta_details.index = "#";
		meta_mutex.unlock();
		return;
	}
	else
		cout << "Processing A Meta File" << endl;

	//Anne's Code for Reading Meta File and Windows Directives
	int nb_lines_read = 1;
	data[0] = std::string(line);
	while (std::getline(fp_meta, line))
	{
		if (line.empty())
			continue; // skip empty lines

		data[nb_lines_read] = std::string(line);
		nb_lines_read++;

		if (nb_lines_read == 7) break; // done reading current sample
	}

#if defined (_WIN32) || defined (_WIN64)
	i = CreateDirectory(data[6].c_str(), NULL);
	if (i == ERROR_ALREADY_EXISTS)
	{
		cerr << "Directory " << data[6].c_str() << " already exists" << endl;
	}
	else if (i == ERROR_PATH_NOT_FOUND)
	{
		cerr << "Path not found" << endl;
	}
#else
	i = mkdir(data[6].c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(i == -1)
	{
		cout << "The directory already Exists..." << endl;
		//cout << "Writing Output in Current Directory..." << endl;
		//data[6] = "./";
	}
#endif


	/*
	//cout << line << endl << endl; 
	//Original Code**************************************************
	char *str = new char[line.length() + 1];
	strcpy(str, line.c_str());
	
	char *pch = strtok(str, "\"\t ");
  	while(pch != NULL)
	{
		//printf ("%s\n", pch);
		data[i++] = string(pch);
		pch = strtok (NULL, "\"\t ");
	}	

	i = mkdir(data[6].c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	if(i == -1)
	{
		cout << "The directory already Exists..." << endl;
		//cout << "Writing Output in Current Directory..." << endl;
		//data[6] = "./";
	}
	//Original Code****************************************************
	*/

	int primer_len, min_barcode;
	meta_details.index = data[0];
	
	//meta_details.fprimer = data[1];
	cout << "Processing Forward Primer File" << endl;
	//primer_len = 
	concat_primers(data[1], meta_details.fprimer, FORWARD_ADAPTER_LENGTH,
		meta_details.forward_primer_map, meta_details.complement_fprimer_map,
			primer_len, min_barcode);
	meta_details.forward_offset = primer_len;
	meta_details.forward_minbarcode = min_barcode; 	

	//meta_details.rprimer = data[2];
	cout << "Processing Reverse Primer File" << endl;
	//primer_len = 
	concat_primers(data[2], meta_details.rprimer, REVERSE_ADAPTER_LENGTH,
		meta_details.reverse_primer_map, meta_details.complement_rprimer_map,
			primer_len, min_barcode);
	meta_details.reverse_offset = primer_len;
	meta_details.reverse_minbarcode = min_barcode;

	meta_details.reference = data[3];
	meta_details.first_read = data[4];
	meta_details.second_read = data[5];
	meta_details.output_path = data[6];

	string complement_fprimer = reverse_complement(meta_details.fprimer);
	create_automata(meta_details.fprimer, meta_details.forward_primer);
	create_automata(complement_fprimer, meta_details.complement_fprimer);

	string complement_rprimer = reverse_complement(meta_details.rprimer);
	create_automata(meta_details.rprimer, meta_details.reverse_primer);
	create_automata(complement_rprimer, meta_details.complement_rprimer);

	//delete [] str;

	cout << "Reference File = " << meta_details.reference << endl;
	if(DEBUG == 99)
	{
		cout << "Forward Primers = " << meta_details.fprimer << endl;
		cout << "Complement FPrimer = " << complement_fprimer << endl;
		cout << "Forward Offset Length = " << meta_details.forward_offset << endl;
		cout << endl;
		cout << "Reverse Primers = " << meta_details.rprimer << endl;
		cout << "Complement RPrimer = " << complement_rprimer << endl;
		cout << "Reverse Offset Length = " << meta_details.reverse_offset << endl;
		cout << endl;
	}
	cout << "Read1 Input File = " << meta_details.first_read << endl;
	cout << "Read2 Input File = " << meta_details.second_read << endl;
	cout << "Output Directory = " << meta_details.output_path << endl;
	//cout << endl;

	meta_mutex.unlock();
	return;
}


//Create a Job for starting Alignment process 
void thread_function(int threadno)
{
	//string detail_file = meta_details.index + "_Details.txt";
	//string sam_output_file = meta_details.index + "_All.sam";
	//fp_detail.open(detail_file.c_str(), ofstream::out);

	//cout << endl << "Starting Thread " << (threadno + 1) << endl;
	int jobno = 0;

	while(true)
	{
		meta_data meta_details;
		read_meta_file(meta_details, threadno, jobno);
		if(meta_details.index.length() == 0)
		{
			//cout << endl << "Returning From Thread " << (threadno + 1) << endl;
			return;
		}
		if(meta_details.index.length() == 1 && meta_details.index[0] == '#')
		{
			continue;
		}
		//string blastn = "output-" + logstr.str() + ".txt";
		//fp_blastn.open(blastn.c_str(), ofstream::out);

		vector<pair<string, string> > reference;
		vector<reference_index> refindex;	//vector of reference
	
		string sam_output_file = meta_details.output_path + meta_details.index + "_All.sam";

		if(REFERENCELESS_ALIGNMENT == 0)
		{
			input_reference(reference, meta_details.reference, sam_output_file);

			kmer_inverted_index(reference, refindex);	
		}
		
		if(ALTERNATIVE_RUN == true)
		{
			cout << "*Going to Preprocess Thread*" << endl;
			preprocess_thread(meta_details, refindex);
		}
		else
		{
			cout << "*SHMPrep Analyzing Sequence*" << endl;
			aligner_thread(meta_details, refindex);
		}

		for(int i = 0; i < reference.size(); i++)
		{
			reference[i].first.clear();
			reference[i].second.clear();
		}

		for(int i = 0; i < refindex.size(); i++)
		{
			refindex[i].ref.clear();
			refindex[i].rev.clear();
			refindex[i].name.clear();
		
			delete [] refindex[i].index;
			delete [] refindex[i].revind;
		}

		jobno += 1;
		init_variables();
		//cout << endl << "Completed a Job by Thread " << (threadno + 1) << endl;
	}
}

void print_help()
{
        cout << "Use the following options to run SHMPrep:" << endl;
        cout << "-f: To specify the name of META file (Ex: -f meta.txt)" << endl;
        cout << "-n: To specify the number of Threads (Ex: -n2)" << endl;
	cout << "-r: To enable Referenceless Alignment (Ex: -r)" << endl;
	cout << "-b: To process BARCODED sequences for Referenceless Alignment (Ex: -b)" << endl;
	cout << "-l: To apply Local Alignment for overlapping reads (Ex: -l)" << endl;
	cout << "-e: To exclude the Alignment with INDELs (Ex: -e)" << endl;
	cout << "-q: To print the final output in FASTQ format (Ex: -q)" << endl;
	cout << "-p: To print BARCODES found in the output read name (Ex: -p)" << endl;
	cout << "-s: If enabled then Alignment will include BARCODE/PRIMERs (Ex: -s)" << endl;
	cout << "-d: If enabled then both CONSCOUNT and DUPCOUNT will be shown (Ex: -d)" << endl;
        cout << "-Q: Optional filter for BASIC Correction by per base quality (Ex: -Q20)" << endl;
        cout << "-M: Optional filter for Mean Quality of the sequence (Ex: -M25)" << endl;
	cout << "-X: To specify the INDEL Penalty for alignment (Ex: -X9)" << endl;
	cout << "-C: to filter out Alignments with CONSCOUNT less than given number (Ex: -C5)" << endl;
	cout << "-D: To filter out Alignments with DUPCOUNT less than given number (Ex: -D10)" << endl;
	cout << "-F: To set the length of Forward Adapter (Ex: -F24) as Prefix in Read1" << endl;
	cout << "-O: To set the length of Reverse Adapter (Ex: -O24) as Prefix in Read2" << endl;
	cout << "-K: To set the size of KBAND in Global Alignment (Ex: -K10)" << endl;
	cout << "-P: To specify the minimum similarity score for PRIMERs (Ex: -P70)" << endl;
	cout << "-L: To specify Alignment Identity to apply Smith Waterman to the Overlapped region (Ex: -L50)" << endl;
	cout << "-G: To specify Alignment Identity to apply Global Alignment to a Consensus (Ex: -G90)" << endl;
        cout << "-h: To print this help menu" << endl << endl;
	//cout << "Pending: SW Alignment for Overlap and Option for Stripping off Barcode/Primers" << endl;
	//cout << "What about base N? Should we consider it for error correction. Previously returned -1" << endl;
	//cout << "Pending: Reorganize functions in multiple files and documentation" << endl;
	//cout << "Should I keep the -s option in the debugging option?" << endl;
}

void prepare_input(int argc, char *argv[])
{
	char c;

	while((c = getopt(argc, argv, "?harbleqpsdf:n:Q:M:X:D:F:O:K:L:G:")) != -1)
        {
                switch(c)
                {
                        case 'f':
				meta_file = string(optarg);
                                break;
			case 'n':
				MAXTHREAD = (int) atol(optarg);
				break;
			case 'r':
				REFERENCELESS_ALIGNMENT = 1;
				break;
			case 'b':
				BARCODED_SEQUENCE = 1;
				break;
			case 'l':
				LOCALFLAG = LOCAL;
				break;
			case 'e':
				EXCLUDE_INDELS = 1;
				break;
			case 'q':
				FASTQOUTPUT = FASTQ_OUTPUT;
				break;
			case 'p':
				PRINTDETAIL = PRINT_DETAIL;
				break;
			case 's':
				STRIPOFF = (int) -1 * STRIP_OFF;
				break;
			case 'd': 
				DUPCOUNT = 1;
				break;
			case 'a':
				ALTERNATIVE_RUN = true;
				break;
			case 'Q':
				BASIC_CORRECTION = 1;
				PERBASE_QSCORE = (int) atol(optarg);
				break;
			case 'M':
				GLOBAL_QSCORE = (int) atol(optarg);
				break;
			case 'X':
				GAP = (int) atol(optarg);
				break;
			case 'C':
				FILTER_LOW_CONS_READS = (int) atol(optarg);
				break;
			case 'D': 
				FILTER_LOW_FREQ_READS = (int) atol(optarg);
				break;
			case 'F':
				FORWARD_ADAPTER_LENGTH = (int) atol(optarg);
				break;
			case 'O':
				REVERSE_ADAPTER_LENGTH = (int) atol(optarg);
				break;
			case 'K':
				GLOBAL_KBAND = (int) atol(optarg);
				break;
			case 'P':
				PRIMER_QUALITY = (float) atol(optarg);
				break;
			case 'L':
				OVERLAP_QUALITY = (float) atol(optarg);
				LOCALFLAG = LOCAL;
				break;
			case 'G':
				GLOBAL_QUALITY = (float) atol(optarg);
				break;
			case 'h':
                        case '?':
                        default:
                                print_help();
                                exit(0);
		}
	}

	if(meta_file.length() == 0)
	{
		print_help();
		exit(0);
	}
}

int main(int argc, char *argv[])
{
	/*
	string primer = "", sequence = "ACNGTSACWGT";
	vector<string> primer_list;
	int index = 0;
	cout << "Original: " << sequence << endl << endl;
	generate_iupac_primers(index, sequence, primer, primer_list);
	exit(0);
	*/
	prepare_input(argc, argv);
	string header;

	cout << "KMER = " << KMER << ", SEED = " << SEED << ", MAXREAD = " << MAXREAD << endl;
	fp_meta.open(meta_file.c_str(), ifstream::in);
	getline(fp_meta, header);

	int x, y, z, novalue;
	FILE* fstream1 = fopen("matching_score.txt", "r");
	for(int i = 0; i < 51; i++)
	{
		for(int k = 0; k < 51; k++)
		{
			novalue = fscanf(fstream1, "%d,%d,%d", &x, &y, &z);
			//cout << "x = " << x << ", y = " << y << ", z = " << z << ", quality = " << char('!' + z) << endl;
			matching_score[x][y] = min(41, z);
		}
	}
	fclose(fstream1);

	FILE* fstream2 = fopen("mismatch_score.txt", "r");
	for(int i = 0; i < 51; i++)
	{
		for(int k = 0; k < 51; k++)
		{
			novalue = fscanf(fstream2, "%d,%d,%d", &x, &y, &z);
			//cout << "x = " << x << ", y = " << y << ", z = " << z << ", quality = " << char('!' + z) <<  endl;
			mismatch_score[x][y] = min(41, z);
		}
	}
	fclose(fstream2);

	log_quality['!'] = 0;
	for(char ch = '!' + 1; ch <= 'J'; ch++)
	{
		log_quality[ch] = log10(1 - pow(10.0, ('!' - ch) / 10.0));
		log_error[ch] = log10(pow(10.0, ('!' - ch) / 10.0) / 3.0);
		//cout << ch << ":: " << (1 - pow(10.0, ('!' - ch) / 10.0)) << ":: " << log_quality[ch] << endl;
		//cout << ch << "!: " << pow(10.0, ('!' - ch) / 10.0) << ":: " << log_error[ch] << endl;
	}


	//return 0;//for testing the functions in this file

	thread_function(0);

	fp_meta.close();
	return 0;
}


