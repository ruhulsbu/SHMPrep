#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

char complement(char ch)
{
	if(ch == 'A')
		return 'T';
	if(ch == 'C')
		return 'G';
	if(ch == 'G')
		return 'C';
	if(ch == 'T')
		return 'A';

	return ch;
}


int map_value(char ch)
{
	if(ch == 'A')
		return 0;
	if(ch == 'C')
		return 1;
	if(ch == 'G')
		return 2;
	if(ch == 'T')
		return 3;

	return -1;
}

void init_matrix(cell **matrix)
{
        matrix = new cell *[2 * FRAGMENT_SIZE + 5];
        for(int i = 0; i < 2 * FRAGMENT_SIZE + 5; i++)
        {
                matrix[i] = new cell[FRAGMENT_SIZE + 5];
        }
}

void remove_matrix(cell **matrix)
{
	for(int i = 0; i < 2 * FRAGMENT_SIZE + 5; i++)
        {
                delete [] matrix[i];
                matrix[i] = NULL;
        }

        delete [] matrix;

        matrix = NULL;
}

int similarity(char x, char y)
{
	if(x == y)
		return 1 * WEIGHT;
	else 
		return MISMATCH;
}

int min(int x, int y)
{
	if(x < y)
		return x;
	else 
		return y;
}

int max(int x, int y) 
{
	if(x < y)
		return y;
	else
		return x;
}

void upper_case(string& str)
{
	char *array = new char [str.length() + 1];
	strcpy(array, str.c_str());
	for(int i = 0; i < str.length(); i++)
		array[i] = toupper(array[i]);
	
	str.clear();
	str = string(array);
	delete [] array;
}

void reverse_str(string &str)
{
	int i, k;
	char ch;
	
	//cout << "Forward = " << str.substr(0, 80) << endl;
	for(i = 0, k = str.length() - 1; i < k; i++, k--)
	{
		ch = str[i];
		str[i] = str[k];
		str[k] = ch;
	}
	//cout << "Reverse = " << str.substr(str.length() - 80, 80) << endl;

}

string reverse_complement(string& str)
{
	int i, j;
	char x, y;
	char *input = new char [str.length() + 1];
	char *output = new char [str.length() + 1];
	strcpy(input, str.c_str());
	if(DEBUG == 1) cout << str << endl;	
	for(i = str.length() - 1, j = 0; i >= 0; i--, j++)
	{
		x = input[i];
		if(x == 'A')
			y = 'T';
		else if(x == 'T')
			y = 'A';
		else if(x == 'C')
			y = 'G';
		else if(x == 'G')
			y = 'C';
		else 
			y = x;//might also be = x here

		output[j] = y;
		if(DEBUG == 1) cout << x << " " << y << endl; 	
	}

	output[j] = '\0';
	if(DEBUG == 1) printf("%s\n", output);
	string rc(output);
	if(DEBUG == 1)
		cout << "Reverse(" << str << ") = " << rc << endl;
	
	delete [] input;
	delete [] output;

	return rc;
}

char calculate_basequal(char first, char second, int flag)
{
	int x = first - '!';
	int y = second - '!';

	if(x == 0)
		return second;
	if(y == 0)
		return first;
	
	if(flag == 1)
	{
		return (matching_score[x][y] + '!');
	}
	else
	{
		return (mismatch_score[x][y] + '!');
	}
}

bool validate_alignment(vector<reference_index>& refindex, int ref_ind, vector<pair<char, char> >& alignment, 
			int ref_start, int read_start, string& read)
{
        for(int i = 0; i < alignment.size(); i++)
        {
                //cout << "Refs: " << ref_start << " = " << alignment[i].first << " VS " << refindex[ref_ind].ref.at(ref_start) << endl;
                //cout << "Read: " << read_start << " = " << alignment[i].second << " VS " << read.at(read_start) << endl;

                if(alignment[i].first != '-')
                {
			//cout << "Refs: " << ref_start << " = " << alignment[i].first << " VS " << refindex[ref_ind].ref.at(ref_start) << endl;
                        assert(alignment[i].first == refindex[ref_ind].ref.at(ref_start));
                        ref_start += 1;
                }
                if(alignment[i].second != '-')
                {
			//cout << "Read: " << read_start << " = " << alignment[i].second << " VS " << read.at(read_start) << endl;
                        assert(alignment[i].second == read.at(read_start));
                        read_start += 1;
                }
        }

	cout << "alignment validate/validation is done properly" << endl;
}

bool validate_chain(vector<reference_index>& refindex, int ref_ind, fragment_alignment& fragment_alignment_info, string& read)
{
	int ref_start = fragment_alignment_info.ref_start;
	int read_start = fragment_alignment_info.read_start;
	for(int i = 0; i < fragment_alignment_info.alignment.size(); i++)
	{
		//cout << ref_start << " = " << fragment_alignment_info.alignment[i].first << " VS " << refindex[ref_ind].ref.at(ref_start) << endl;
		//cout << read_start << " = " << fragment_alignment_info.alignment[i].second << " VS " << read.at(read_start) << endl;

		if(fragment_alignment_info.alignment[i].first != '-')
		{
			assert(fragment_alignment_info.alignment[i].first == refindex[ref_ind].ref.at(ref_start));
			ref_start += 1;
		}
		if(fragment_alignment_info.alignment[i].second != '-')
		{
			assert(fragment_alignment_info.alignment[i].second == read.at(read_start));
			read_start += 1;
		}
	}
	cout << "chain validate/validation is done properly" << endl;
}

//Reorder the sequence based on the forward primer
//can this function be optimized by merging forward and complement suffix automata?
//Also, replace the erroneous substring with the correct one in the alignment sequence for higher dupcount
string correct_read_direction(string& readseq_first, string& quality_first, meta_data& meta_details, ofstream& fp_detail)
{
	int forward_primer_found = -1;
	int refindex, queryindex;
	string secondary_search_string, confirm_primer_string;
	string forward_offset_string, complement_offset_string;
	string forward_primer_string, complement_fprimer_string;
	string primer_string;

	//Find the maximum length substring between offset_string and primer automata
	//generate the consecutive 5-mers and find the primer where these exists
	//return the intersection of all the primer index set of the second step

	if(readseq_first.length() < meta_details.reverse_offset)
		return "";
	

	{
		forward_offset_string = readseq_first.substr(0, meta_details.reverse_offset);
		fp_detail << "Forward Offset String = " << forward_offset_string << endl;
		forward_primer_string = find_primer(meta_details.reverse_primer, 
					forward_offset_string, refindex, queryindex, fp_detail);
		fp_detail << "Forward Primer in First sequence Found from Automata = " << forward_primer_string << endl << endl;

		if(map_primers.find(forward_primer_string) != map_primers.end())
		{
			//cout << "Return Forward Primer Mapping Using Kmer Match = " << 
			//	return_primer_using_kmer_match(forward_primer_string, meta_details.forward_primer_map) << endl;


			return return_primer_using_kmer_match(forward_primer_string, meta_details.reverse_primer_map);
			//return forward_primer_string;
		}
	}
	/*
	string primer_string = return_primer_using_kmer_match(forward_primer_string, 
							meta_details.forward_primer_map);
	if(primer_string.length() > 0)
		return primer_string;//Testing Forward Reads in First Read Set Only
	*/
	{
		complement_offset_string = readseq_first.substr(readseq_first.length() - meta_details.reverse_offset, 
					meta_details.reverse_offset);
		fp_detail << "Complement Offset String = " << complement_offset_string << endl;
		complement_fprimer_string = find_primer(meta_details.complement_rprimer, 
							complement_offset_string, refindex, queryindex, fp_detail);
		fp_detail << "Complement FPrimer in First sequence  Found from Automata = " << complement_fprimer_string << endl << endl;

		if(map_primers.find(complement_fprimer_string) != map_primers.end())
		{
			readseq_first = reverse_complement(readseq_first);
			reverse_str(quality_first);

			//cout << "Return Forward Complement Primer Mapping Using Kmer Match = " << 
			//	return_primer_using_kmer_match(complement_fprimer_string, meta_details.complement_fprimer_map) << endl;
			//exit(0);//for testing this function
			
			//cout << "Complement FPrimer Found!" << endl;
			return return_primer_using_kmer_match(complement_fprimer_string, meta_details.complement_rprimer_map);
			//return reverse_complement(complement_fprimer_string);
		}

	}
		
	if(forward_primer_string.length() >= complement_fprimer_string.length())
	{
		if(forward_primer_string.length() > PRIMER_KMER_SIZE)
			return return_primer_using_kmer_match(forward_primer_string, meta_details.reverse_primer_map);
	}
	else
	{
		if(complement_fprimer_string.length() > PRIMER_KMER_SIZE)
			primer_string = return_primer_using_kmer_match(complement_fprimer_string, 
							meta_details.complement_rprimer_map);
		if(primer_string.length() > 0)
		{
			readseq_first = reverse_complement(readseq_first);
			reverse_str(quality_first);
			return primer_string;
		}
	}
	//exit(0);//for testing this function
	return "";
}



void preprocess_reads(meta_data& meta_details, vector<reference_index>& refindex, int threadno)
{
	/*
	fp_detail: It prints output of all the major computation to a Details.txt file for each thread per job
	fp_fastq: prints read_name, sam CIGAR information, alignment and quality string for each thread per job
	*/
	
	ofstream fp_detail, fp_fastq;
		
	string detail_file = meta_details.output_path + meta_details.index + "_" + to_string(threadno) + "_Details.txt";
	string fastq_output_file = meta_details.output_path + meta_details.index + "_" + to_string(threadno) + "_All.detail";

	if(DEBUG != 99)
		detail_file = "/dev/null";

	fp_detail.open(detail_file.c_str(), ofstream::out);
	fp_fastq.open(fastq_output_file.c_str(), ofstream::out);

	//variables for first and second read
	string read_name_first, readseq_first;
	string read_name_second, readseq_second;
	string quality_first, quality_second;

	int map = 0, count = 0, cant_map = 0, name_length, return_val;
	string return_fprimer_string, return_rprimer_string;
	string forward_unique_barcode, reverse_unique_barcode;
	//MAXREAD = 100000;

	while(read_sequence_file(read_name_first, readseq_first, quality_first, 
				read_name_second, readseq_second, quality_second) > 0)
	{
		if(readseq_first.length() < meta_details.forward_offset ||
			readseq_first.length() < meta_details.reverse_offset)
			continue;
		if(readseq_second.length() < meta_details.forward_offset ||
			readseq_second.length() < meta_details.reverse_offset)
			continue;

		if(count >= MAXREAD && MAXREAD != 0) break;
		count += 1;
		if(count <= MINREAD) continue;
		if(count % 100000 == 0)
		{
			cout << endl << "Total " << count << " Reads Aligned from " << meta_details.index << 
				"by Thread " << threadno << endl << endl;
		}

		fp_detail << "####################################################################################################################" << endl;
		fp_detail << "####################################################################################################################" << endl;
		fp_detail << endl;

		
		if(DEBUG == 99)
		{
			fp_detail << "First ) " << read_name_first << endl;
			fp_detail << "Read String = " << readseq_first << endl;//.substr(0, 80) << endl;
			fp_detail << "Quality =  :: " << quality_first << endl << endl;

			fp_detail << "Second) " << read_name_second << endl;
			fp_detail << "Read String = " << readseq_second << endl;//.substr(0, 80) << endl;
			fp_detail << "Quality =  :: " << quality_second << endl << endl;

			//referenceless_alignment(readseq_first, readseq_second, quality_first, 
			//			quality_second, matrix, fp_detail);
		}

		/*	
		transform each base of the two reads to upper case. find the referenceless alignment of 
		read pair and read quality, and keep the information in readsequence and readquality variables
		*/
		upper_case(readseq_first);		
		upper_case(readseq_second);
	
		fp_detail << "Zero Step: Find Forward Primer" << endl;
		fp_detail << "-------------------------------" << endl;

		return_fprimer_string = correct_alignment_direction(readseq_first, quality_first, meta_details, fp_detail);

		if(return_fprimer_string.length() > 0)
		{
			return_rprimer_string = correct_read_direction(readseq_second, quality_second, meta_details, fp_detail);	
		}
		else
		//	continue;			
		{
			fp_detail << "Analyze Second Read for FPrimer:" << endl;
			return_fprimer_string = correct_alignment_direction(readseq_second, quality_second, meta_details, fp_detail);
			if(return_fprimer_string.length() > 0)
			{
				//exit(0);//for test
				return_rprimer_string = correct_read_direction(readseq_first, quality_first, meta_details, fp_detail);
				//exit(0);//for test
				if(return_rprimer_string.length() > 0)
					continue;	
				
				string swap_string;
				swap_string = readseq_first;
				readseq_first = readseq_second;
				readseq_second = swap_string;

				swap_string = quality_first;
				quality_first = quality_second;
				quality_second = swap_string;
			}
			else
				continue;
		}
		
		forward_unique_barcode = readseq_first.substr(0, map_primers[return_fprimer_string] % PRIMER_ID_CONSTANT);
		reverse_unique_barcode = readseq_second.substr(0, map_primers[return_rprimer_string] % PRIMER_ID_CONSTANT);

		fp_detail << "Forward Primers: " << return_fprimer_string << endl;
		fp_detail << "Forward Barcode: " << forward_unique_barcode << endl;
		fp_detail << "Reverse Primers: " << return_rprimer_string << endl;
		fp_detail << "Reverse Barcode: " << reverse_unique_barcode << endl;

		fp_fastq << read_name_first << endl;
		fp_fastq << return_fprimer_string << endl;
		fp_fastq << forward_unique_barcode << endl;
		fp_fastq << readseq_first << endl;
		fp_fastq << quality_first << endl;

		fp_fastq << read_name_second << endl;
		fp_fastq << return_rprimer_string << endl;
		fp_fastq << reverse_unique_barcode << endl;
		fp_fastq << readseq_second << endl;
		fp_fastq << quality_second << endl;

	}

	fp_detail.close();
	fp_fastq.close();

	return;
}


int analyze_sequences(vector<pair<string, string> >& sequence_list, string& final_string, 
				string& final_quality, ofstream& fp_detail)
{
	unordered_map<int, int> length_frequency;
	int seq_length, max_count = 0;
	int TESTQUALITY = DEBUG;


	for(int l = 0; l < sequence_list.size(); l++)
	{
		if(TESTQUALITY == 99)
		{
			fp_detail << (l + 1) << ": " << sequence_list[l].first << endl;
			fp_detail << (l + 1) << ": " << sequence_list[l].second << endl;
		}
		seq_length = sequence_list[l].first.length();

		if(length_frequency.find(seq_length) == length_frequency.end())
			length_frequency[seq_length] = 0;
		length_frequency[seq_length] += 1;		
	}

	for(auto length: length_frequency)
	{
		if(max_count < length.second)
		{
			seq_length = length.first;
			max_count = length.second;
		}
	}

	vector<pair<string, string> > seq_information;
	for(int l = 0; l < sequence_list.size(); l++)
	{
		if(seq_length == sequence_list[l].first.size())
			seq_information.push_back(make_pair(sequence_list[l].first, sequence_list[l].second));
	}

	if(TESTQUALITY == 99)
		fp_detail << "Max Frequency = " << max_count << ", For Length = " 
				<< seq_length << ", Check = " << seq_information.size() << endl;

	char seq_base[] = {'A', 'T', 'C', 'G'}; //What about character N?

	for(int k = 0; k < seq_length; k++)
	{
		float base_quality[] = {0.0, 0.0, 0.0, 0.0};
		int count_base[] = {0, 0, 0, 0};

		//cout << endl << "Base Number: " << k << endl;
				
		for(int l = 0; l < seq_information.size(); l++)
		{
			char base = seq_information[l].first[k];//base k at string l
			char quality = seq_information[l].second[k];
			//cout << "Base: " << base << ", Qual: " << quality << endl;
			//cout << "Quality and Error = " << log_quality[quality] << "::: "
                        //                                        << log_error[quality] << endl;

			if(quality < '!' or quality > 'J')
				assert(false);

			for(int p = 0; p < 4; p++)
			{
				if(base == seq_base[p])
				{
					base_quality[p] += log_quality[quality];
					count_base[p] += 1;
				}
				else
					base_quality[p] += log_error[quality];

				//cout << base_quality[p] << "::: " << log_quality[quality] << "::: " 
				//				<< log_error[quality] << endl;
			}
		}

		float corrected_quality = 0.0, corrected_error = 0.0; //[] = {0.0, 0.0, 0.0, 0.0};
		float total_quality = 0.0, max_quality = 0.0;
			
		char corrected_base, corrected_qchar;

		for(int p = 0; p < 4; p++)
		{
			total_quality += pow(10.0, base_quality[p]);					
		}

		for(int p = 0; p < 4; p++)
		{
			corrected_quality = pow(10.0, base_quality[p]) / total_quality;
			//cout << "For Base = " << seq_base[p] << ", Quality = " << corrected_quality << endl; 
		
			if(max_quality < corrected_quality && count_base[p] > 0)
			{
				corrected_base = seq_base[p];
				max_quality = corrected_quality;
			}
		}

		//cout << "Corrected Base = " << corrected_base << ", Max Quality = " << max_quality << endl;

		corrected_error = -10.0 * log10(max(0.000079, 1.0 - max_quality));
		corrected_qchar = '!' + (int)(corrected_error + 0.5);

		if(corrected_qchar < '!' or corrected_qchar > 'J')
			assert(false);

		//cout << "Corrected Base: " << corrected_base << ", Corrected Qual: " << corrected_qchar << endl;

		final_string += corrected_base;
		final_quality += corrected_qchar;

	}

	return seq_information.size();
}


void preprocess_thread(meta_data& meta_details, vector<reference_index>& refindex)
{
	std::thread thread_object[MAX_THREAD];
	MAXTHREAD = min(MAXTHREAD, MAX_THREAD);

	fp_read_first.open(meta_details.first_read.c_str(), ifstream::in);
	fp_read_second.open(meta_details.second_read.c_str(), ifstream::in);


        //Launch a group of threads
        for(int i = 0; i < MAXTHREAD; ++i)
        {
                thread_object[i] = std::thread(preprocess_reads, std::ref(meta_details), std::ref(refindex), i);
        }

        std::cout << "Launched from the Thread Function\n";

        //Join the threads with the main thread
        for (int i = 0; i < MAXTHREAD; ++i) {
                thread_object[i].join();
        }

	fp_read_first.close();
	fp_read_second.close();

	//create file streams to merge the output from all the threads of a single job
	ofstream fp_fasta, fp_detail, fp_duplicate;

	string detail_file = meta_details.output_path + meta_details.index + "_Details.txt";
	if(DEBUG != 99)
		detail_file = "/dev/null";
	fp_detail.open(detail_file.c_str(), ofstream::out);

	string duplicate_file = meta_details.output_path + meta_details.index + "_All.detail";
	fp_duplicate.open(duplicate_file.c_str(), ofstream::out);

        string fasta_file = meta_details.output_path + meta_details.index + "_All.fasta";
        fp_fasta.open(fasta_file.c_str(), ofstream::out);


	//merge the output of all threads to _All.detail and _All.sam file
	for(int i = 0; i < MAXTHREAD; i++)
	{
		string input;
		ifstream fp_input;
		string fastq_output_file = meta_details.output_path + meta_details.index + "_" + to_string(i) + "_All.detail";
		fp_input.open(fastq_output_file.c_str(), ifstream::in);

		while(getline(fp_input, input))
		{
			fp_duplicate << input << endl;
		}
		fp_input.close();
		remove(fastq_output_file.c_str());
	}

	fp_duplicate.close();

	ifstream fp_input;
	fp_input.open(duplicate_file.c_str(), ifstream::in);

	struct read_pair {
		string read_name_first;
		string readseq_first;
		string quality_first;
		string forward_barcode;
		string forward_primer;

		string read_name_second;
		string readseq_second;
		string quality_second;
		string reverse_barcode;
		string reverse_primer;
	};

	string key, input;
	unordered_map<string, vector<read_pair>> hash_map;

	while(getline(fp_input, input))
        {
		read_pair reads;
		//Read first four lines: read_name, cigar string, alignment_string and quality_score
                //getline(fp_input, reads.read_name_first);
		reads.read_name_first = input;
 		getline(fp_input, reads.forward_primer);
      		getline(fp_input, reads.forward_barcode);
	 	getline(fp_input, reads.readseq_first);
		getline(fp_input, reads.quality_first);
	
		getline(fp_input, reads.read_name_second);
	 	getline(fp_input, reads.reverse_primer);
        	getline(fp_input, reads.reverse_barcode);
      		getline(fp_input, reads.readseq_second);
		getline(fp_input, reads.quality_second);
		
		key = reads.forward_barcode + "_" + reads.reverse_barcode;
		if(hash_map.find(key) == hash_map.end())
		{
			vector<read_pair> read_vector;
			hash_map[key] = read_vector;
		}
		hash_map[key].push_back(reads);	
		
		//fp_detail << "BARCODE = " << key << endl;
		//fp_detail << reads.read_name_first << endl;
		//fp_detail << reads.read_name_second << endl;
		 
	
	}

	cout << "Total Hash Map Size = " << hash_map.size() << endl;
	fp_detail << "Correction:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl << endl;

	
	int TESTQUALITY = DEBUG;
        //create a matrix for k-band semi global alignment program for this thread only
        cell **matrix = NULL;
        //init_matrix(matrix);//does not work
        matrix = new cell *[2 * FRAGMENT_SIZE + 5];
        for(int i = 0; i < 2 * FRAGMENT_SIZE + 5; i++)
        {
                matrix[i] = new cell[FRAGMENT_SIZE + 5];
        }

	for(auto& elem: hash_map)
	{
		string 	barcode_seq = elem.first;
		vector<read_pair> vector_read = elem.second;

		vector<pair<string, string> > sequence_list;
		for(int i = 0; i < vector_read.size(); i++)
			sequence_list.push_back(make_pair(vector_read[i].readseq_first, vector_read[i].quality_first));

		if(TESTQUALITY == 99){
			fp_detail << "#############################################################################" << endl;
			fp_detail << "Collaps Both the reads for BARCODE" << endl;
			fp_detail << "----------------------------------" << endl;
			fp_detail << ">\nFirst Read: BARCODE=" << barcode_seq 
					<< "|FPrimer=" << vector_read[0].forward_primer 
					<< "|RPrimer=" << vector_read[0].reverse_primer 
					<< "|TOTAL_READS=" << sequence_list.size() <<  endl;
		}
		string first_string = "", first_quality = "";
		int conscount = analyze_sequences(sequence_list, first_string, first_quality, fp_detail);
		if(TESTQUALITY == 99)
		{
			fp_detail << ">BARCODE=" << barcode_seq << "|CONSCOUNT=" << conscount <<  endl;
			fp_detail << first_string << endl;

			fp_detail << "+" << endl;	
			fp_detail << first_quality << endl;
		}


	
		sequence_list.clear();
		for(int i = 0; i < vector_read.size(); i++)
			sequence_list.push_back(make_pair(vector_read[i].readseq_second, vector_read[i].quality_second));

		if(TESTQUALITY == 99)
			fp_detail << ">\nScond Read: BARCODE=" << barcode_seq << "|TOTAL_READS=" << sequence_list.size() <<  endl;
	
		string second_string = "", second_quality = "";
		conscount = analyze_sequences(sequence_list, second_string, second_quality, fp_detail);
		if(TESTQUALITY == 99)
		{
			fp_detail << ">BARCODE=" << barcode_seq << "|CONSCOUNT=" << sequence_list.size() <<  endl;
			fp_detail << second_string << endl;

			fp_detail << "+" << endl;	
			fp_detail << second_quality << endl;
		}
		
		fp_detail << endl << endl;

		string readsequence = "", readquality = "";
		int return_val;
		return_val = referenceless_alignment(first_string, second_string, first_quality, second_quality,
                                                readsequence, readquality, matrix, fp_detail);
                if(return_val == -1)
			continue;

		string sam_output_name = "";
		string read_name_first = vector_read[0].read_name_first;
		if(sam_output_name.length() == 0)
                {
                        for(int i = 1; i < read_name_first.length() && read_name_first[i] != ' '; i++)
                                sam_output_name += read_name_first[i];
                }

		fp_fasta << ">" << sam_output_name << "|BARCODE=" << barcode_seq << "|CONSCOUNT=" << sequence_list.size() << endl;
		fp_fasta << readsequence << endl;
		fp_fasta << "+" << endl;
		fp_fasta << readquality << endl;
		//cout << read_name << endl << final_string << endl;
		//assert(false);

	}	

	fp_fasta.close();
	fp_input.close();
	fp_detail.close();

}


