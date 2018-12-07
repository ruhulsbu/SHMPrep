#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"


ifstream fp_read_first;
ifstream fp_read_second;
std::mutex read_mutex;

/*
We overlapped a read pair and created a new overlapped read which we simply refere as 'read' here. 
This function find the alignment direction between the read and reference. It also print the detail
in fp_detail stream for debugging. To find the direction, we created two KMER lookups table for the 
reference in ints both the 'forward' and 'reverse' direction. We count the KMER of the read for the 
first SEED number of bases using the two tables. Depending on the count of KMERs in either direction 
(forward_map_count >= reverse_map_count) we return the direction of alignment. 
*/

bool find_alignment_direction(string& read, int offset, reference_index& reference, ofstream& fp_detail)
{
	time_t start, end;

	time(&start);

	bool flag = true;
	long hash_key = 0;
	int map_val;
	int forward_map_count = 0;
	int reverse_map_count = 0; 

	for(int k = offset; k < read.length() && k < SEED + offset; k++)
	{
		if(flag == true && k + KMER > read.length())
			break;
		for(int l = k, end = k, i = KMER - 2; l < end + KMER - 1 && flag == true; l++, k++, i--)
		{
			map_val = map_value(read.at(l));
			if(map_val == -1)
				break;

			hash_key += base_power_kmer[i][map_val];
			
			//cout << "For character " << read.at(k) << ", at position = " << k <<
			//	", the hash_key = " << hash_key << endl;
		}

		map_val = map_value(read.at(k));
		if(map_val == -1)
		{
			//cout << "Encountered invalid character N ########################" << endl;
			flag = true;
			hash_key = 0;
			continue;
		}
		else
			flag = false;

		hash_key = hash_key * BASE + map_val;

		if(reference.index[hash_key] != -1)
		{
			forward_map_count += 1;
		}
		if(reference.revind[hash_key] != -1)
		{
			reverse_map_count += 1;
		}
		
		//cout << "For character " << read.at(k) << ", at position = " << k <<
		//		", the hash_key = " << hash_key << endl;
	
		map_val = map_value(read.at(k - KMER + 1));
		hash_key -= base_power_kmer[KMER - 1][map_val];
	}

	if(forward_map_count >= reverse_map_count)
		return true;
	else 
		return false;

}

/*
This function takes the following input:
readseq: read to be aligned
read_name: name of the read sequence
fragmnt_alignment_info: data structure to hold alignment information
quality: quality of the read sequence
final_result: vector to hold raead name, sam CIGAR information of an alignment, alignment and quality
meta_data: data structure that holds the meta information of an alignment task
**matrix: Two dimensional matrix used to calculate k-band semi-golobal alignment program
fp_sam: output stream to write the sam_format
fp_detail: output stream to write the detail inforation of the alignment for debugging
*/

int align_read_to_reference(string& readseq, string& read_name, fragment_alignment& final_alignment_info, string& quality,
				vector<string>& final_result, vector<reference_index>& refindex, meta_data& meta_details,  
				cell **matrix, string fprimer, string rprimer, ofstream& fp_sam, ofstream &fp_detail)
{
	int match_info, global_match = -1, indpos, find;
	int match = 0, gaps = 0, mismatch = 0, match_index;
	int ref_start, read_start, ref_end, read_end;
	int str1_end, str2_end, direction;

	int primer_length, except_primer;
	int ref_pointer, read_pointer;
	vector<pair<char, char> > alignment;
	
	bool flag;	
	string refgenome;
	char first, second;
	string final_read;
		
	for(int i = 0; i < refindex.size(); i++)
	{
		//find the forward primer length and reference length except primer, except_primer
		primer_length = fprimer.length();
		except_primer = fprimer.length() + refindex[i].ref.length();
		refgenome = fprimer + refindex[i].ref + rprimer;

		/*
		return if the reference length is ALIGNMENT_COVERAGE times shorter or larger than readseq length
		*/
		if(1.0 * refgenome.length() / readseq.length() > ALIGNMENT_COVERAGE)
			return -1;
		if(1.0 * readseq.length() / refgenome.length() > ALIGNMENT_COVERAGE)
			return -1;



		//initialize variable and find the direction of the alignment
		match = mismatch = gaps = str1_end = str2_end = 0;
		flag = find_alignment_direction(readseq, fprimer.length(), refindex[i], fp_detail);

		//for reverse direction, find the reverse complement of read sequence
		if(flag == false)
		{
			readseq = reverse_complement(readseq);
			reverse_str(quality);
		}
		
		//take equal length of substring from reference and read to find the alignment
		if(refgenome.length() > readseq.length())
			refgenome = refgenome.substr(0, readseq.length());
		else
		{
			readseq = readseq.substr(0, refgenome.length());
			quality = quality.substr(0, refgenome.length());
		}

		if(DEBUG == 99)
		{
			fp_detail << "Reference length = " << refgenome.length() << endl;
			fp_detail << "Reference = " << refgenome << endl;//.substr(0, 80) << endl;
			fp_detail << "Readseq   length = " << readseq.length() << endl;
			fp_detail << "Readseq   = " << readseq << endl;//.substr(0, 80) << endl;
		}

		find_kband_similarity(refgenome, readseq, alignment, str1_end, str2_end, matrix);

		if(DEBUG == 99)
		{
			fp_detail << "Reference alignment length = " << str1_end << endl;
			fp_detail << "Read      alignment length = " << str2_end << endl;
			fp_detail << "Quality      string length = " << quality.length() << endl;
		}

		final_alignment_info.alignment.clear();
		final_alignment_info.quality = "";
		final_read = "";

		ref_pointer = read_pointer = 0;
	
		/*	
		confirm that the alignment coordinates agree  to the actual read and reference
		following the forward primer as well as the start of the reverse primer
		*/
		for(int k = alignment.size() - 1; k >= 0 && ref_pointer < except_primer; k--)
		{
			assert(ref_pointer <= str1_end);
			assert(read_pointer <= str2_end);

			if(ref_pointer < primer_length)
			{
				if(alignment[k].first != '-')
					ref_pointer += 1;
				//else
				//	return -1;//primer should have perfect match
			
				if(alignment[k].second != '-')
					read_pointer += 1;

				continue;
			}
			
			final_alignment_info.alignment.push_back(alignment[k]);
			if(alignment[k].first == alignment[k].second && alignment[k].first != '-')
				match += 1;
			if(alignment[k].first != alignment[k].second && alignment[k].first != '-')
				mismatch += 1;
			if(alignment[k].first == '-' || alignment[k].second == '-')
				gaps += 1;

			if(alignment[k].first != '-')
				ref_pointer += 1;

			if(alignment[k].second != '-')
			{
				final_read += readseq.at(read_pointer);
				final_alignment_info.quality += quality.at(read_pointer);
				read_pointer += 1;
			}
			
		}

		/*
		if reference pointer does not match the start of reverse primer then the alignment is 
		not valid. Also, if ORIGINAL_PERCENT_MATCH number of read bases out of 100 bases does 
		not match then we return -1, meaning invalid alignment 
		*/
		if(ref_pointer != except_primer)
		{
			fp_detail << "Ref Pointer did not Match Except Pointer" << endl;
			return -1;
		}
		if(EXCLUDE_INDELS == 1 && gaps > 0)
		{
			fp_detail << "Returning due to INDELS in the alignment" << endl;
			return -1;
		}
		if(100.0 * match / read_pointer < GLOBAL_PERCENT_MATCH)
		{
			fp_detail << "Returnging for Less Percent Match: " << (100.0 * match / read_pointer) << endl;
			return -1;
		}

		read_start = 0;
		read_end = final_read.length();

		direction = 1;
		ref_start = 0;
		ref_end = except_primer - primer_length;
		//ref_start = primer_length;
		//ref_end = ref_pointer;

		//store the alignment information in to the final_alignment_information data structure
		final_alignment_info.total_len = final_alignment_info.alignment.size();
                final_alignment_info.identity_match = match;
                final_alignment_info.ref_start = ref_start;
                final_alignment_info.read_start = read_start;
                final_alignment_info.ref_end = ref_end;//total_ref_ind
                final_alignment_info.read_end = read_end;//total_read_ind
                final_alignment_info.ref_ind = i;
                final_alignment_info.read_dir = direction;
                final_alignment_info.gaps = gaps;
                final_alignment_info.mismatches = mismatch;

		//Create the sam_format for the alingment
		sam_format(final_alignment_info, refindex, final_read, read_name, final_result, fp_detail);

		fp_sam << final_result[0];
		for(int k = 1; k < final_result.size(); k++)
               	{
                       	fp_sam << "\t" << final_result[k];
                       	//cout << i << ": " << output[k] << endl;
               	}		
		fp_sam << endl;

		if(DEBUG == 99)
		{
			fp_detail << "Alignment Found For: " << read_name << endl;
			fp_detail << read_name << "," << refindex[i].name << "," << KMER << ",";
	       		fp_detail << final_alignment_info.read_start << "," << final_alignment_info.read_end << ",";
		       	fp_detail << final_alignment_info.ref_start << "," << final_alignment_info.ref_end << ",";
		       	fp_detail << final_alignment_info.gaps << "," << final_alignment_info.mismatches << ",";
		       	fp_detail << final_alignment_info.identity_match << "," << final_alignment_info.total_len << ",";
		       	fp_detail << (100.00 * final_alignment_info.identity_match / final_alignment_info.total_len) << ",";
		       	fp_detail << final_alignment_info.read_dir << "," << readseq.length() << ",";
		       	fp_detail << (100.00 * final_alignment_info.total_len / readseq.length()) << endl << endl << endl;
		}
		
		/*
	  	for(int k = 0; k < 80; k++)
			cout << complement(alignment[k].first);
		cout << endl;
		for(int k = 0; k < 80; k++)
			cout << complement(alignment[k].second);
		cout << endl << endl;

		cout << "Reference = " << refindex[i].ref << endl;
		cout << "Input ===== " << readseq << endl << endl;
		
		for(int k = 0; k < 80; k++)
			cout << final_alignment_info.alignment[k].first;
		cout << endl;
		for(int k = 0; k < 80; k++)
			cout << final_alignment_info.alignment[k].second;
		cout << endl << endl;

		cout << "Reference = " << refindex[i].ref << endl;
		cout << "Input ===== " << readseq <<endl << endl;
		*/

		//generate base pair <reference_i, read_i> from the alignment using print_alignment function
		if(direction == 1)
		{
			print_alignment(final_alignment_info.alignment, refindex[i].ref, final_read, ref_start, 
				read_start, 1, final_alignment_info, true, fp_detail);
				
			if(DEBUG == 99)
			{
				fp_detail << "reference = " << refindex[i].ref.substr(final_alignment_info.ref_start, 80) << endl;
				fp_detail << "for_read  = " << readseq.substr(0, 80) << endl << endl;
			}		
		}
		else
			assert(false);

	
		alignment.clear();
		return 0;
	}

	
	return 0;
}

//Scan the read pairs along with with their qaulity strings and read names
int read_sequence_file(string& read_name_first, string& readseq_first, string& quality_first,
                                string& read_name_second, string& readseq_second, string& quality_second)
{
	string input;
	read_mutex.lock();
	if(getline(fp_read_first, read_name_first) && getline(fp_read_second, read_name_second))
        {
                getline(fp_read_first, readseq_first);
		if(FORWARD_ADAPTER_LENGTH >= 1)
			readseq_first = readseq_first.substr(FORWARD_ADAPTER_LENGTH, readseq_first.length() - FORWARD_ADAPTER_LENGTH);

                getline(fp_read_first, input);
                getline(fp_read_first, quality_first);

                getline(fp_read_second, readseq_second);
		if(REVERSE_ADAPTER_LENGTH >= 1)
			readseq_second = readseq_second.substr(REVERSE_ADAPTER_LENGTH, readseq_second.length() - REVERSE_ADAPTER_LENGTH);

                getline(fp_read_second, input);
                getline(fp_read_second, quality_second);

		PAIRED_READ_COUNT += 1;		
		read_mutex.unlock();
		return 1;
	}
	else
	{
		cout << "Reached The End of File !" << endl;
		read_mutex.unlock();
		return -1;
	}
}

string return_primer_using_kmer_match(string &sequence, unordered_map<string, vector<int> >& primer_map)
{
	int primer_id_array[primer_vector.size()];
	memset(primer_id_array, 0, sizeof(int) * primer_vector.size());
	
	int primer_id, max_count = 0;
	string primer_kmer;

	int length = sequence.length() - PRIMER_KMER_SIZE + 1;
	//cout << "Primer Map Size = " << primer_map.size() << ", Sequence = " << sequence << endl;

	for(int i = 0; i < length; i++)
        {
                primer_kmer = sequence.substr(i, PRIMER_KMER_SIZE);
                if(primer_map.find(primer_kmer) == primer_map.end())
			continue;
		
                vector<int> primer_list = primer_map[primer_kmer];
               
		//cout << primer_kmer << " = " << primer_list.size() << endl;
 
		for(int k = 0; k < primer_list.size(); k++)
		{
			//cout << "\tPrimer ID From Vector = " << primer_list[k] << endl;
			primer_id = primer_list[k];
			primer_id_array[primer_id] += 1;
		}
        }

	primer_id = 0;
	for(int i = 0; i < primer_vector.size(); i++)
	{
		if(primer_id_array[i] > max_count)
		{
			max_count = primer_id_array[i];
			primer_id = i;
		}
	}

	//cout << "Max Count = " << max_count << ", Primer ID = " << primer_id << endl;
	return primer_vector[primer_id];
}


//Find the existence of reverse primer in the alignment
string identify_reverse_primer(string& alignment, string& quality, meta_data& meta_details, ofstream& fp_detail)
{
	string reverse_offset_string, reverse_primer_string;
	int refindex, queryindex;

	if(alignment.length() < meta_details.reverse_offset)
		return "";

	{
		//It seems the rprimer is the prefix
		//reverse_offset_string = alignment.substr(alignment.length() - meta_details.reverse_offset, 
		//			meta_details.reverse_offset);
		reverse_offset_string = alignment.substr(alignment.length() - meta_details.reverse_offset,
					meta_details.reverse_offset - meta_details.reverse_minbarcode);
		fp_detail << "Reverse Offset String = " << reverse_offset_string << endl;
		reverse_primer_string = find_primer(meta_details.complement_rprimer, 
					reverse_offset_string, refindex, queryindex, fp_detail);
		fp_detail << "Reverse Primer in First sequence  Found from Automata = " << reverse_primer_string << endl << endl;

		if(map_primers.find(reverse_primer_string) != map_primers.end())
		{
			//cout << "Return Reverse Primer Mapping Using Kmer Match = " << 
			//	return_primer_using_kmer_match(reverse_primer_string, meta_details.complement_rprimer_map) << endl;

			//exit(0);//for testing this function

			//return return_primer_using_kmer_match(reverse_primer_string, meta_details.complement_rprimer_map);
			return reverse_primer_string;
		}
	}

	if(reverse_primer_string.length() <= PRIMER_KMER_SIZE)
		return "";
	string reverse_primer_candidate = return_primer_using_kmer_match(reverse_primer_string, meta_details.complement_rprimer_map);
	float similarity = (100.0 * reverse_primer_string.length()) / reverse_primer_candidate.length();
	fp_detail << "Reverse Primer Candidate = " << reverse_primer_candidate << " with Score = " << similarity << endl;

	if(similarity > PRIMER_QUALITY)
		return reverse_primer_candidate;	
		
	similarity = (100.0 * find_longest_common_subsequence(reverse_primer_candidate, reverse_offset_string)) / reverse_primer_candidate.length();
	fp_detail << "From LCS = " << similarity << endl;
	fp_detail << "Reverse Primer Candidate = " << reverse_primer_candidate << " Offset String = " << reverse_offset_string << endl;


	if(similarity > LCS_PRIMER_SIMILARITY)
		return reverse_primer_candidate;	
	
	/*
	{
		//It seems the rprimer is the prefix
		reverse_offset_string = alignment.substr(0, meta_details.reverse_offset);
		fp_detail << "Reverse Offset String = " << reverse_offset_string << endl;
		reverse_primer_string = find_primer(meta_details.reverse_primer, 
					reverse_offset_string, refindex, queryindex, fp_detail);
		fp_detail << "Reverse Primer in First sequence  Found from Automata = " << reverse_primer_string << endl << endl;

		if(map_primers.find(reverse_primer_string) != map_primers.end())
		{
			alignment = reverse_complement(alignment);
			reverse_str(quality);
			
			//cout << "Return Reverse Primer Mapping Using Kmer Match = " << 
			//	return_primer_using_kmer_match(reverse_primer_string, meta_details.complement_rprimer_map) << endl;

			//exit(0);//for testing this function
			return reverse_primer_string;
		}
	}
	*/


	//exit(0);//for testing this function
	return "";
}

//Reorder the sequence based on the forward primer
//can this function be optimized by merging forward and complement suffix automata?
//Also, replace the erroneous substring with the correct one in the alignment sequence for higher dupcount
string correct_alignment_direction(string& readseq_first, string& quality_first, meta_data& meta_details, ofstream& fp_detail)
{
	int forward_primer_found = -1;
	int refindex, queryindex;
	
	string secondary_search_string, confirm_primer_string;
	string forward_offset_string, complement_offset_string;
	string forward_primer_string, complement_fprimer_string;
	string primer_string, offset_string;

	//Find the maximum length substring between offset_string and primer automata
	//generate the consecutive 5-mers and find the primer where these exists
	//return the intersection of all the primer index set of the second step

	if(readseq_first.length() < meta_details.forward_offset)
		return "";
	

	{
		//forward_offset_string = readseq_first.substr(0, meta_details.forward_offset);
		forward_offset_string = readseq_first.substr(meta_details.forward_minbarcode, 
						meta_details.forward_offset - meta_details.forward_minbarcode);
		fp_detail << "Forward Offset String = " << forward_offset_string << endl;
		forward_primer_string = find_primer(meta_details.forward_primer, 
					forward_offset_string, refindex, queryindex, fp_detail);
		fp_detail << "Forward Primer in First sequence Found from Automata = " << forward_primer_string << endl << endl;

		if(map_primers.find(forward_primer_string) != map_primers.end())
		{
			//cout << "Return Forward Primer Mapping Using Kmer Match = " << 
			//	return_primer_using_kmer_match(forward_primer_string, meta_details.forward_primer_map) << endl;


			//return return_primer_using_kmer_match(forward_primer_string, meta_details.forward_primer_map);
			return forward_primer_string;
		}
	}

	/*
	string primer_string = return_primer_using_kmer_match(forward_primer_string, 
							meta_details.forward_primer_map);
	if(primer_string.length() > 0)
		return primer_string;//Testing Forward Reads in First Read Set Only
	*/
	{
		//complement_offset_string = readseq_first.substr(readseq_first.length() - meta_details.forward_offset, 
		//			meta_details.forward_offset);
		complement_offset_string = readseq_first.substr(readseq_first.length() - meta_details.forward_offset,
					meta_details.forward_offset - meta_details.forward_minbarcode);
		fp_detail << "Complement Offset String = " << complement_offset_string << endl;
		complement_fprimer_string = find_primer(meta_details.complement_fprimer, 
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
			//return reverse_complement(return_primer_using_kmer_match(complement_fprimer_string, meta_details.complement_fprimer_map));
			return reverse_complement(complement_fprimer_string);
		}

	}


	string forward_primer_candidate = "";		
	float similarity = 0.0;

	if(forward_primer_string.length() >= complement_fprimer_string.length())
	{
		if(forward_primer_string.length() > PRIMER_KMER_SIZE)
		{
			forward_primer_candidate = return_primer_using_kmer_match(forward_primer_string, meta_details.forward_primer_map);
			similarity = (100.0 * forward_primer_string.length()) / forward_primer_candidate.length();
			fp_detail << "Forward Primer Candidate = " << forward_primer_candidate << " with Score = " << similarity << endl;

			if(similarity > PRIMER_QUALITY)
				return forward_primer_candidate;

			offset_string = forward_offset_string;	
		}
		else 
			return "";
	}
	else
	{
		if(complement_fprimer_string.length() > PRIMER_KMER_SIZE)
			primer_string = return_primer_using_kmer_match(complement_fprimer_string, 
							meta_details.complement_fprimer_map);
		if(primer_string.length() > 0)
		{
			readseq_first = reverse_complement(readseq_first);
			reverse_str(quality_first);
			forward_primer_candidate = reverse_complement(primer_string);
			
			similarity = (100.0 * complement_fprimer_string.length()) / forward_primer_candidate.length();
			fp_detail << "Forward Primer Candidate = " << forward_primer_candidate << " with Score = " << similarity << endl;

			if(similarity > PRIMER_QUALITY)
				return forward_primer_candidate;	

			offset_string = complement_offset_string;

		}
		else
			return "";
	}

	similarity = (100.0 * find_longest_common_subsequence(forward_primer_candidate, offset_string)) / forward_primer_candidate.length();
	fp_detail << "From LCS = " << similarity << endl;
	fp_detail << "Forward Primer Candidate = " << forward_primer_candidate << ", Offset String = " << offset_string << endl;

	if(similarity > LCS_PRIMER_SIMILARITY)
		return forward_primer_candidate;	

	//exit(0);//for testing this function
	return "";
}


/*
Each thread call this function with the meta details, reference and the threadno
*/
//How to ensure the direction of the referenceless alignment?
void align_reads(meta_data& meta_details, vector<reference_index>& refindex, int threadno)
{
	/*
	fp_detail: It prints output of all the major computation to a Details.txt file for each thread per job
	fp_fastq: prints read_name, sam CIGAR information, alignment and quality string for each thread per job
	fp_sam: prints the alignment information for each read in sam file for each thread per job
	*/
	
	ofstream fp_detail, fp_sam, fp_fastq;
		
	string detail_file = meta_details.output_path + meta_details.index + "_" + to_string(threadno) + "_Details.txt";
	string sam_output_file = meta_details.output_path + meta_details.index + "_" + to_string(threadno) + "_All.sam";
	string fastq_output_file = meta_details.output_path + meta_details.index + "_" + to_string(threadno) + "_All.detail";

	if(DEBUG != 99)
		detail_file = "/dev/null";

	fp_detail.open(detail_file.c_str(), ofstream::out);
	fp_sam.open(sam_output_file.c_str(), ofstream::out);
	fp_fastq.open(fastq_output_file.c_str(), ofstream::out);

	//variables for first and second read
	string input_first, read_name_first;
	string readseq_first, refgenome_first;
	
	string input_second, read_name_second;
	string readseq_second, refgenome_second;

	string quality_first, quality_second;
	string original_read_name, sam_output_name;

	string readsequence, readquality;

	fragment_alignment fragment_alignment_first;
	fragment_alignment fragment_alignment_second;
	vector<string> sam_output_first;
	vector<string> sam_output_second;

	int map = 0, count = 0, cant_map = 0;
	int name_length, return_val, read_pair_analyzed = 0;
	string return_fprimer_string, return_rprimer_string;
	string forward_unique_barcode, reverse_unique_barcode;

	//create a matrix for k-band semi global alignment program for this thread only
	cell **matrix = NULL;
	//init_matrix(matrix);//does not work
	matrix = new cell *[2 * FRAGMENT_SIZE + 5];
        for(int i = 0; i < 2 * FRAGMENT_SIZE + 5; i++)
        {
                matrix[i] = new cell[FRAGMENT_SIZE + 5];
        }

	//variables to hold forward and reverse primer information
	//string reverse_fprimer, reverse_rprimer;
	//reverse_fprimer = reverse_complement(meta_details.fprimer);
	//reverse_rprimer = reverse_complement(meta_details.rprimer);
	//cout << "Reverse Fprimer: " << reverse_fprimer << endl;
	//cout << "Reverse Rprimer: " << reverse_rprimer << endl;
	
	//if(DEBUG == 99)
	//	fp_detail << "size of matrix = " << sizeof(matrix) << endl;

	//use MINREAD and MAXREAD to debug a set of particular reads by using a single thread
	//MINREAD = 50000;
	//MAXREAD = 10000;
	int strtok_count;

	//scan all information for a pair of reads
        cout << "What the Hack!" << endl;
	while(read_sequence_file(read_name_first, readseq_first, quality_first, 
				read_name_second, readseq_second, quality_second) > 0)
	{
                //cout << read_name_first << endl;
                //cout << readseq_first << endl;
                //cout << read_name_second << endl;
                //cout << readseq_second << endl;

		if(readseq_first.length() < meta_details.forward_offset ||
			readseq_first.length() < meta_details.reverse_offset)
			continue;
		if(readseq_second.length() < meta_details.forward_offset ||
			readseq_second.length() < meta_details.reverse_offset)
			continue;
		
		strtok_count = 0;
		original_read_name = "";
		sam_output_name = "";
		/*
	        char *pch = strtok(str, ": ");
	        while(pch != NULL)
        	{
                	//printf ("%s\n", pch);
			if(strtok_count >= 4 && strtok_count <= 6)			
                		sam_output_name += string(pch) + "_";

                	pch = strtok (NULL, "\"\t ");
			strtok_count += 1;
        	}
		*/

		for(int i = 0; i < read_name_first.length() && read_name_first[i] != ' '; i++)
		{
			if(strtok_count >= 4)
			{
				if(read_name_first[i] >= '0' && read_name_first[i] <= '9')
					sam_output_name += read_name_first[i];
				else
					sam_output_name += '_';
			}

			if(read_name_first[i] == ':')
				strtok_count += 1;
		}

		if(sam_output_name.length() == 0)
		{
			for(int i = 1; i < read_name_first.length() && read_name_first[i] != ' '; i++)
				sam_output_name += read_name_first[i];
		}	
		original_read_name = sam_output_name;

		//else
		//	fp_detail << "File Name Chosen = " << sam_output_name << endl;

		//use MINRAD and MAXREAD for debugging purpose and print the information of 100,000th alignment
		if(count >= MAXREAD && MAXREAD != 0) break;
		count += 1;
		if(count <= MINREAD) continue;
		if(count % 100000 == 0)
		{
			//cout << endl << "Total " << count << " Reads Aligned from " << meta_details.index << 
			//	" by Thread " << (threadno + 1) << endl;
			cout << "Thread " << (threadno + 1) << " Aligned " << count << 
				" Reads, Contributing to Total: " << PAIRED_READ_COUNT << endl;
		}

		fp_detail << "####################################################################################################################" << endl;
		fp_detail << "####################################################################################################################" << endl;
		fp_detail << count << ") Given Dataset:" << original_read_name << endl << endl;

		
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

		return_val = referenceless_alignment(readseq_first, readseq_second, quality_first, quality_second,
                	readsequence, readquality, matrix, fp_detail);
		if(return_val == -1)
			continue;

	
		fp_detail << "Fourth Step: Find Forward Primer" << endl;
		fp_detail << "--------------------------------" << endl;

		return_fprimer_string = correct_alignment_direction(readsequence, readquality, meta_details, fp_detail);
		if(return_fprimer_string.length() == 0)
			continue;

	
		fp_detail << "Fifth Step: Find Reverse Primer" << endl;
		fp_detail << "-------------------------------" << endl;

		return_rprimer_string = identify_reverse_primer(readsequence, readquality, meta_details, fp_detail);
		if(return_rprimer_string.length() == 0)
		{
			fp_detail << "Conclusion: Valid Alignment is Not Found" << endl << endl;
			continue; 
			//exit(0);//for testing this function
		}
		else
			fp_detail << "Conclusion: A Valid Alignment is Found" << endl << endl;
	
		forward_unique_barcode = readsequence.substr(0, map_primers[return_fprimer_string] % PRIMER_ID_CONSTANT);
		reverse_unique_barcode = readsequence.substr(readsequence.length() - map_primers[return_rprimer_string] 
							% PRIMER_ID_CONSTANT, readsequence.length());

		fp_detail << "Forward Primers: " << return_fprimer_string << endl;
		fp_detail << "Forward Barcode: " << forward_unique_barcode << endl;
		fp_detail << "Reverse Primers: " << return_rprimer_string << endl;
		fp_detail << "Reverse Barcode: " << reverse_unique_barcode << endl;
		fp_detail << "Readsequence Length: " << readsequence.length() << endl;
		fp_detail << "Read quality Length: " << readquality.length() << endl;	

		string primer_barcode;
		sam_output_name = "";

		{
			if(return_fprimer_string.length() > 0)
			{
				sam_output_name += "|VPRIMER=" + primer_name[return_fprimer_string];
				primer_barcode += primer_name[return_fprimer_string] + "&";
			}
			else
				primer_barcode += "&";

			if(return_rprimer_string.length() > 0)
			{
				sam_output_name += "|CPRIMER=" + primer_name[return_rprimer_string];
				primer_barcode += primer_name[return_rprimer_string] + "#";
			}
			else
				primer_barcode += "#";
			
			if(forward_unique_barcode.length() > 0 &&
				reverse_unique_barcode.length() > 0)
			{
				//sam_output_name += "|VBARCODE=" + forward_unique_barcode;
				//sam_output_name += "|CBARCODE=" + reverse_unique_barcode;
				if(PRINTDETAIL == PRINT_DETAIL)
					sam_output_name += "|BARCODE=" + forward_unique_barcode + "/" + reverse_unique_barcode;
				primer_barcode += forward_unique_barcode + "&" + reverse_unique_barcode;
			}
			else
			{
				if(forward_unique_barcode.length() > 0)
				{
					if(PRINTDETAIL == PRINT_DETAIL)
						sam_output_name += "|BARCODE=" + forward_unique_barcode;
					primer_barcode += forward_unique_barcode + "&";
				}
				else
					primer_barcode += "&";

				if(reverse_unique_barcode.length() > 0)
				{
					if(PRINTDETAIL == PRINT_DETAIL)
						sam_output_name += "|BARCODE=" + reverse_unique_barcode;
					primer_barcode += reverse_unique_barcode;
				}
			}
		}

		//If -S enabled then it will not be stripped off; By default strip off	
		if(STRIPOFF == STRIP_OFF)
		{	
			//stripping off barcode adn primers from the alignment and quality
			int prefix_strip, total_strip;
			prefix_strip = forward_unique_barcode.length() + return_fprimer_string.length();
			total_strip = prefix_strip + (reverse_unique_barcode.length() + return_rprimer_string.length());
			if(readsequence.length() < prefix_strip + total_strip)
				continue;

			readsequence = readsequence.substr(prefix_strip, readsequence.length() - total_strip);
			readquality = readquality.substr(prefix_strip, readquality.length() - total_strip);
			fp_detail << "Final Alignment: " << readsequence << endl;
			fp_detail << "Final Quality  : " << readquality << endl;
		}
		/*
		//Original****************************************************************************************
		return_val = referenceless_alignment(readseq_first, readseq_second, quality_first, quality_second, 
					readsequence, readquality, matrix, fp_detail);	

		//continue if no overlap is found
		if(return_val == -1)
			continue;
		//Original***************************************************************************************
		*/

		if(REFERENCELESS_ALIGNMENT == 1)// && CONSERVATIVE_REFERENCELESS_ALIGNMENT == 0)
		{		
			
			fp_fastq << original_read_name << sam_output_name << endl;
			//fp_fastq << unique_barcode << "_" << return_fprimer_string << endl;
			fp_fastq << primer_barcode << endl;
			//fp_fastq << return_fprimer_string << endl;

			fp_fastq << readsequence << endl;
			fp_fastq << readquality << endl;

			fp_detail << "Alignment Title: " << original_read_name <<  sam_output_name << endl;
			read_pair_analyzed += 1;
			continue;

		}

		//align overlapped read to reference 
		fragment_alignment final_alignment_info;
		vector<string> final_result;

		if(STRIPOFF == STRIP_OFF)
		{
			return_val = align_read_to_reference(readsequence, original_read_name, final_alignment_info, readquality,
					final_result, refindex, meta_details,  matrix, "", "", fp_sam, fp_detail);
		}
		else
		{
			//Strip off Barcode for alignment with the primers at both end
			if(BARCODED_SEQUENCE == 1)
			{
				readsequence = readsequence.substr(forward_unique_barcode.length(), readsequence.length() - 
							reverse_unique_barcode.length());
				readquality = readquality.substr(forward_unique_barcode.length(), readquality.length() - 
							reverse_unique_barcode.length());
				fp_detail << endl << "Stripped Off Barcode, Now Calling for Alignment" << endl;
			}

			return_val = align_read_to_reference(readsequence, original_read_name, final_alignment_info, 
						readquality, final_result, refindex, meta_details,  matrix, 
						return_fprimer_string, return_rprimer_string, fp_sam, fp_detail);

		}
		//continue if no alignment is found
		if(return_val == -1)
			continue;

		//report all the information in the fp_fastq file stream
		assert(final_result[10].length() == final_result[9].length());

		if(BARCODED_SEQUENCE == 1)
		{
			fp_fastq << final_result[0] << "|CIGAR=" << final_result[5] << sam_output_name << endl;
			fp_fastq << primer_barcode << endl;
			//fp_fastq << forward_unique_barcode << "_" << reverse_unique_barcode << endl;
		}
		else
		{
			fp_fastq << final_result[0] << endl;
			fp_fastq << primer_barcode << "=" << final_result[5] << endl;
		}
		fp_fastq << final_result[9] << endl;
		fp_fastq << final_result[10] << endl;

		final_result.clear();
	        final_alignment_info.alignment.clear();

		read_pair_analyzed += 1;
		//break;	
	}

	//cout << "Total Read Pairs Analyzed: " << count << endl;
	//cout << "Successfully Overlapped: " << read_pair_analyzed << endl;

	if(DEBUG == 99)
	{
		fp_detail << endl << "Overall Statistics - " << endl;
		fp_detail << "Total read = " << count << endl;
	}

	//remove_matrix(matrix);
	for(int i = 0; i < 2 * FRAGMENT_SIZE + 5; i++)
        {
                delete [] matrix[i];
                matrix[i] = NULL;
        }

        delete [] matrix;
        matrix = NULL;


	fp_detail.close();
	fp_sam.close();
	fp_fastq.close();

	return;
}

/*
Create a number of threads for each job as defined in the meta_details.
refindex holds the reference string and its KMER look up table in both directions
Ideally, the size of this vector is always one, i.e. only one reference
*/
void aligner_thread(meta_data& meta_details, vector<reference_index>& refindex)
{
	std::thread thread_object[MAX_THREAD];
	MAXTHREAD = min(MAXTHREAD, MAX_THREAD);

	fp_read_first.open(meta_details.first_read.c_str(), ifstream::in);
	fp_read_second.open(meta_details.second_read.c_str(), ifstream::in);


        //Launch a group of threads
        for(int i = 0; i < MAXTHREAD; ++i)
        {
                thread_object[i] = std::thread(align_reads, std::ref(meta_details), std::ref(refindex), i);
        }

        std::cout << "Launched from the Thread Function\n";

        //Join the threads with the main thread
        for (int i = 0; i < MAXTHREAD; ++i) {
                thread_object[i].join();
        }

	fp_read_first.close();
	fp_read_second.close();

	//create file streams to merge the output from all the threads of a single job
	ofstream fp_detail, fp_duplicate, fp_fasta, fp_sam;

	string detail_file = meta_details.output_path + meta_details.index + "_Details.txt";
	if(DEBUG != 99)
		detail_file = "/dev/null";
	fp_detail.open(detail_file.c_str(), ofstream::out);

	string duplicate_file = meta_details.output_path + meta_details.index + "_All.detail";
	fp_duplicate.open(duplicate_file.c_str(), ofstream::out);


	string fasta_file = meta_details.output_path + meta_details.index + "_All.fasta";
	if(FASTQOUTPUT == FASTQ_OUTPUT)
		fasta_file = meta_details.output_path + meta_details.index + "_All.fastq";
	fp_fasta.open(fasta_file.c_str(), ofstream::out);

	string sam_file = meta_details.output_path + meta_details.index + "_All.sam";
        fp_sam.open(sam_file.c_str(), ofstream::out | ofstream::app);

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
		//remove(fastq_output_file.c_str());

		ifstream sam_input;
		string sam_output_file = meta_details.output_path + meta_details.index + "_" + to_string(i) + "_All.sam";
		sam_input.open(sam_output_file.c_str(), ifstream::in);
		
		while(getline(sam_input, input))
		{
			fp_sam << input << endl;
		}
		sam_input.close();
		//remove(sam_output_file.c_str());
	}

	fp_duplicate.close();
	fp_sam.close();

	ifstream fp_input;
	fp_input.open(duplicate_file.c_str(), ifstream::in);

	//create fasta file with alignment correction , DUPCOUNT and CIGAR string
	for(int k = 0; k < refindex.size() && REFERENCELESS_ALIGNMENT == 0; k++)
	{
		fp_fasta << ">CONSENSUS" << endl;
		fp_fasta << refindex[k].ref << endl;
		if(FASTQOUTPUT == FASTQ_OUTPUT)
		{
			fp_fasta << "+" << endl;
			for(int i = 0; i < refindex[k].ref.length(); i++)
				fp_fasta << "J";
			fp_fasta << endl;
		}
	}	

	//cout << endl;
	for(int k = 0; k < refindex.size() && REFERENCELESS_ALIGNMENT == 0; k++)
	{
		if(BASIC_CORRECTION == 1)
		{
			cout << "Applying Basic Correction" << endl;
			prepare_fasta_output_heuristic(fp_input, fp_fasta, fp_detail, refindex[k].ref, meta_details);
		}
		else
		{
			if(BARCODED_SEQUENCE == 1)
			{
				cout << "Applying Basic Correction" << endl;
                		prepare_fasta_output_heuristic(fp_input, fp_fasta, fp_detail, refindex[k].ref, meta_details);
			}
			else
			{
				cout << "Apply Poisson Binomial Distribution" << endl;
				prepare_fasta_output_distribution(fp_input, fp_fasta, fp_detail, refindex[k].ref, meta_details);
			}
		}
	}

	if(REFERENCELESS_ALIGNMENT == 1)
	{
		cout << "Applying Basic Correction" << endl;
		prepare_fasta_output_heuristic(fp_input, fp_fasta, fp_detail, refindex[0].ref, meta_details);
	}


	fp_input.close();
	fp_detail.close();
	fp_fasta.close();

}


