#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

bool compare_string(pair<string, string> first_pair, pair<string, string> second_pair)
{
        return (strcmp(first_pair.second.c_str(), second_pair.second.c_str()) <= 0);
}


/*
This function apply few heuristics for error correction. The input are:
1. fp_input: file stream from which we can scan the paired read alignment
2. fp_fasta: file stream to print the aligned reads with DUPCOUNT and SAM format
3. fp_detail: to print the log information in details
4. reference: reference string
5. meta_details: configurations in meta_details file
*/
void prepare_fasta_output_heuristic(ifstream& fp_input, ofstream& fp_fasta, ofstream& fp_detail, 
					string& reference, meta_data& meta_details)
{
	int testcount = 0, hashcount = 0;
        string read_name, annotation, cigar, alignment, quality, len;
	string alignment_string, quality_score, detail;
	int total_score, average_score, opt;
	int i, k, count, refindex, readindex;

        unordered_map<string, long long> umap;
	unordered_map<string, vector<int> > fastq;
	vector< pair<string, string> > sequences; 
	unordered_map<string, vector<pair<string, string> > > map_barcode;

	string barcode_seq;

	fp_detail << "###############################################Sequence Anaysis Starts##############################################" << endl;
	
	while(getline(fp_input, read_name))
        {
		
		//Read first four lines: read_name, annotation string, alignment_string and quality_score
                getline(fp_input, annotation);
                getline(fp_input, alignment_string);
		getline(fp_input, quality_score);
		
		//Calculate average quality_score. if it is less than GLOBAL_QSCORE then continue
		total_score = 0;
		for(i = 0; i < quality_score.length(); i++)
		{
			total_score += quality_score.at(i) - '!';
		}

		average_score = total_score / quality_score.length();
		fp_detail << endl << read_name << ": Global Average Score = " << average_score << endl;

		if(average_score < GLOBAL_QSCORE)
			continue;
		testcount += 1;

		if(REFERENCELESS_ALIGNMENT == 1)
		{			
			/*
			If we are doing REFERENCELESS_ALIGNMENT then concatenate cigar and 
			alignment_string and use this concatenated string as a key with the
			number of occurrence as a value which we print as DUPCOUNT.
			*/
			//testcount += 1;
			//Original Implementation**********************************
			if(BARCODED_SEQUENCE == 0)// || BASIC_CORRECTION == 1)
			{
				detail = annotation + "=" + alignment_string;			//quality score cigar+\n
				fp_detail << sequences.size() << ": " << detail << endl;
	                	if(umap.find(detail) == umap.end())
        	        	{
                	        	//umap[alignment] = 1;
                        		sequences.push_back(make_pair(read_name, detail));
                        		umap[detail] = 1;

					if(FASTQOUTPUT == FASTQ_OUTPUT)		//quality score
					{
						vector<int> quality_number;
						for(int k = 0; k < quality_score.length(); k++)
							quality_number.push_back(quality_score[k] - '!');
						fastq[detail] = quality_number;
					}
                		}
                		else
                		{
                        		umap[detail] += 1;
					
					if(FASTQOUTPUT == FASTQ_OUTPUT)		//quality score
					{
						vector<int> quality_number = fastq[detail];
						for(int k = 0; k < quality_number.size(); k++)
						{
							quality_number[k] += (quality_score[k] - '!');
						}
						fastq[detail] = quality_number;

						assert(quality_number.size() == quality_score.length());
					}
                		}
			
			}		
			//Original Implementation**********************************
			else
			{
				barcode_seq = annotation;
				if(map_barcode.find(barcode_seq) == map_barcode.end())
				{
					vector<pair<string, string> > sequence_vector;
					map_barcode[barcode_seq] = sequence_vector;
					sequences.push_back(make_pair(barcode_seq, read_name)); 
				}
				map_barcode[barcode_seq].push_back(make_pair(alignment_string, quality_score));
			}
			continue;
		}

		if(REFERENCELESS_ALIGNMENT == 0 && BARCODED_SEQUENCE == 1)
		{
			//testcount += 1;
			barcode_seq = annotation;
                        if(map_barcode.find(barcode_seq) == map_barcode.end())
                        {
                                vector<pair<string, string> > sequence_vector;
                                map_barcode[barcode_seq] = sequence_vector;
                                sequences.push_back(make_pair(barcode_seq, read_name));
                        }
                        map_barcode[barcode_seq].push_back(make_pair(alignment_string, quality_score));
			continue;
		}
	
		/*
	        We compare each read base to a reference base after the fprimer_length
	        */
	
		//refindex = meta_details.fprimer.length();
		refindex = 0;
		i = readindex = 0;
		alignment = "";
		quality = "";
		len = "";
		
		/*
		Rest of the code deal with the paired read alignment to a reference. In this part, 
		we modify the alignment and the quality string based on heuristics. We read the cigar
		file to identify each base as mtch, mismatch, insertion or deletion and keep it in the
		alignment if it was sequenced with a minimum quality score otherwise we consider the
		base from reference, and vice versa. The details are given below:
		*/

		int split_index = annotation.find("=") + 1;
		cigar = annotation.substr(split_index, annotation.length() - split_index);
		fp_detail << "CIGAR = " << cigar << endl;

		while(i < cigar.length())
		{
			if(cigar.at(i) >= '0' && cigar.at(i) <= '9')
				len += cigar.at(i);
			else
			{
				//Parse the cigar string for mismatch, insertion and deletion
				opt = cigar.at(i);
				count = atoi(len.c_str());
				fp_detail << "Count = " << count << ", and Option = " << opt << endl;
				len = "";			

				fp_detail << "Reference = " << reference << endl;
				fp_detail << "Alignment = " << alignment_string << endl;

				//For the mathcing or mismatching bases
				if(opt == 'M' || opt == 'X')
				{
					//For the next k number of bases
					for(k = 0; k < count; k++)
					{
						//If the read has 'N' then take the base from reference
						if(alignment_string.at(readindex) == 'N')
						{
							alignment += reference.at(refindex);
							quality += quality_score.at(readindex);

						}
						//If the quality score is less than user defined PERBASE_QSCORE
						//then consider the base from reference and quality_score from the read
						//otherwise, consider both the base and quality_score from read
						else if(quality_score.at(readindex) - '!' < PERBASE_QSCORE)
						{
							alignment += reference.at(refindex);
							quality += quality_score.at(readindex);
						}
						else//Do not apply error correction. Might be a mutation!!!
						{
							alignment += alignment_string.at(readindex);
							quality += quality_score.at(readindex);
						}

						//fp_detail << k << "= (" << refindex << ", " << readindex << ")" << endl;

						//Increment the read and reference index
						refindex += 1;
						readindex += 1;
					}	
				}
				//For the inserted bases, consider both the base and quality_score from read, 
				//then increment readindex
				else if(opt == 'I')
				{
					for(k = 0; k < count; k++)
					{
						alignment += alignment_string.at(readindex);
						quality += quality_score.at(readindex);

						readindex += 1;
					}
				}
				//For the deleted bases increment refindex only 
				else if(opt == 'D')
				{
					refindex += count;
				}
				else
				{
					assert(false);
				}
			} 
			
			i += 1;
		}

		/*
		Concatenate ciger and alignment_string and use this concatenated string as a key with the
		number of occurrence as a value which we print as DUPCOUNT.
		*/
		//testcount += 1;
		detail = annotation + "\n" + alignment;
		if(umap.find(detail) == umap.end())
                {
                        //umap[alignment] = 1;
                	sequences.push_back(make_pair(read_name, detail));
			umap[detail] = 1;

			if(FASTQOUTPUT == FASTQ_OUTPUT)		//quality score
			{
				vector<int> quality_number;
				for(int k = 0; k < quality_score.length(); k++)
					quality_number.push_back(quality_score[k] - '!');
				fastq[detail] = quality_number;
			}
                }
                else
                {
                        umap[detail] += 1;

			if(FASTQOUTPUT == FASTQ_OUTPUT)		//quality score
			{
				vector<int> quality_number = fastq[detail];
				for(int k = 0; k < quality_number.size(); k++)
				{
					quality_number[k] += (quality_score[k] - '!');
				}
				fastq[detail] = quality_number;

				assert(quality_number.size() == quality_score.length());
			}
                }

		/*
                fp_fasta << ">" << read_name;

                fp_fasta << "|DUPCOUNT=1";
                fp_fasta << "|CIGAR=" << cigar << endl;
                fp_fasta << alignment << endl;
		*/
        }

	//sort(sequences.begin(), sequences.end(), compare_string);

	//Create the fasta file with read_name, alignment, DUPCOUNT and CIGAR string
	//Original Implementation***************************************************
	if(BARCODED_SEQUENCE == 0)// || BASIC_CORRECTION == 1)
	{
		int split_index;
		string sequence_str = "";

		for(i = 0; i < sequences.size(); i++)
		{
			hashcount += umap[sequences[i].second];
			if(umap[sequences[i].second] < FILTER_LOW_FREQ_READS)
				continue;
			
			split_index = sequences[i].second.find("=") + 1;
			sequence_str = sequences[i].second.substr(split_index, sequences[i].second.length() - split_index);

	        	fp_fasta << ">" << sequences[i].first;
                	fp_fasta << "|DUPCOUNT=" << umap[sequences[i].second];
			if(REFERENCELESS_ALIGNMENT == 0)
                		//fp_fasta << "|CIGAR=" << sequences[i].second << endl;
				fp_fasta << "|CIGAR=" << sequence_str << endl;
			else
				fp_fasta << endl << sequence_str << endl;

			if(FASTQOUTPUT == FASTQ_OUTPUT)
			{
				fp_fasta << "+" << endl;
				string quality_score = "";
				int quality_count = umap[sequences[i].second];
				vector<int> quality_number = fastq[sequences[i].second];
				
				for(int k = 0; k < quality_number.size(); k++)
				{
					int quality = (int) lrint((1.0 * quality_number[k]) / (1.0 * quality_count));
					quality_score += ('!' + quality);
				}
				fp_fasta << quality_score << endl;
			}

		}
	}
	//Original Implementation***************************************************

	//cout << "Reached Here" << endl;
	else
	{
		testcount = 0;
		int TESTQUALITY = -99;
		vector<pair<string, string> > readname_detail;
		unordered_map<string, string> map_dupcount_info;

		for(i = 0; i < sequences.size(); i++)
		{

			unordered_map<int, int> length_frequency;
			int seq_length, max_count = 0;

			barcode_seq = sequences[i].first;//annotation of primers and barcodes
			read_name = sequences[i].second;

			vector<pair<string, string> >& sequence_list = map_barcode[barcode_seq];

			if(TESTQUALITY == 99)
				fp_fasta << ">" << read_name << "|TOTAL=" << sequence_list.size() <<  endl;
		
			for(int l = 0; l < sequence_list.size(); l++)
			{
				if(TESTQUALITY == 99)
				{
					fp_fasta << (l + 1) << ": " << sequence_list[l].first << endl;
					fp_fasta << (l + 1) << ": " << sequence_list[l].second << endl;
				}
				seq_length = sequence_list[l].first.length();

				if(length_frequency.find(seq_length) == length_frequency.end())
					length_frequency[seq_length] = 0;
				length_frequency[seq_length] += 1;
				
			}

			for(auto elements: length_frequency)
			{
				if(max_count < elements.second)
				{
					seq_length = elements.first;
					max_count = elements.second;
				}
			}

			vector<pair<string, string> > seq_information;
			for(int l = 0; l < sequence_list.size(); l++)
			{
				if(seq_length == sequence_list[l].first.size())
				{
					seq_information.push_back(make_pair(sequence_list[l].first, sequence_list[l].second));
				}
			}

			if(seq_information.size() < FILTER_LOW_CONS_READS)
				continue;
			testcount += seq_information.size();

			if(TESTQUALITY == 99)
				fp_fasta << "Max Frequency = " << max_count << ", For Length = " 
					<< seq_length << ", Check = " << seq_information.size() << endl;

			string final_string = "", final_quality = "";
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
				float total_quality = 0.0, max_quality = 0.0, qual;
				char corrected_base = 'N', corrected_qchar;
				
				//implement max_base_quality for the base that occurs more than 95%
				for(int p = 0; p < 4; p++)
				{
					//fp_fasta << "Count Base = " << count_base[p] << ", Quality = " << base_quality[p] << endl;
					qual = 1.0 * count_base[p] / seq_information.size();
					if(qual > 0.95)
					{
						//max_quality = qual;
						max_quality = pow(10.0, base_quality[p] / (seq_information.size() * qual));
						corrected_base = seq_base[p];
					}
				}
				//fp_fasta << "Max Quality = " << max_quality << endl;
				
				if(corrected_base == 'N')
				{	
					max_quality = 0.0;
					for(int p = 0; p < 4; p++)
					{
						if(base_quality[p] > -307)
						{
							total_quality += pow(10.0, base_quality[p]);
							//fp_fasta << "BaseQual = " << base_quality[p] << 
							//	", TotalQual = " << total_quality << endl;
							//assert(total_quality > 0.0);
						}
					}

					for(int p = 0; p < 4; p++)
					{
						if(base_quality[p] > -307 && total_quality > 0.0)
							corrected_quality = pow(10.0, base_quality[p]) / total_quality;
						else
							continue;
						//fp_fasta << "For Base = " << seq_base[p] << ", Quality = " << corrected_quality << endl; 
				
						if(max_quality < corrected_quality && count_base[p] > 0)
						{
							corrected_base = seq_base[p];
							max_quality = corrected_quality;
						}
					}
				}

				corrected_error = -10.0 * log10(max(0.000079, 1.0 - max_quality));
				corrected_qchar = '!' + (int)(corrected_error + 0.5);

				if(corrected_qchar < '!' or corrected_qchar > 'J')
					assert(false); 
				if(corrected_base < 'A' or corrected_base > 'T')
				{
					cout << "Error: " << corrected_base << endl;
					assert(false);
				}


				if(TESTQUALITY == 99)
				{
					fp_fasta << "Corrected Base = " << corrected_base << ", Max Quality = " << max_quality << endl;
					fp_fasta << "Corrected Base: " << corrected_base << ", Corrected Qual: " << corrected_qchar << endl << endl;
				}

				final_string += corrected_base;
				final_quality += corrected_qchar;
 			}

			assert(final_string.length() == final_quality.length());
			if(DUPCOUNT == 0)
			{
				//fp_fasta << ">" << read_name << "|BARCODE=" << barcode_seq << 
				//	"|CONSCOUNT=" << seq_information.size() <<  endl;
				//else
				fp_fasta << ">" << read_name << "|CONSCOUNT=" << seq_information.size() << endl;
				fp_fasta << final_string << endl;
			
				fp_fasta.flush();	
				if(FASTQOUTPUT == FASTQ_OUTPUT)
				{
					fp_fasta << "+" << endl;	
					fp_fasta << final_quality << endl;
				}
				
				hashcount += seq_information.size();
			}
			else
			{
				//read_name = read_name + "|CONSCOUNT=" + to_string(seq_information.size());
				//detail = "|CONSCOUNT=" + to_string(seq_information.size()) + "\n" + final_string + "\n";
				
				int split_index = barcode_seq.find("#");// + 1
				string primer_string = barcode_seq.substr(0, split_index);//, barcode_seq.length() - split_index);
				detail = primer_string + "#" + final_string;
				//detail = final_string; //pRESTO Edited

				if(DEBUG == 99)
					fp_detail << readname_detail.size() << ": " << detail << endl;

                                if(umap.find(detail) == umap.end())
                                {
                                        umap[detail] = DUPCOUNT_CONSTANT + seq_information.size();
					map_dupcount_info[detail] = read_name;
					readname_detail.push_back(make_pair(read_name, detail));
					if(FASTQOUTPUT == FASTQ_OUTPUT)		//quality score
					{
						vector<int> quality_number;
						for(int k = 0; k < final_quality.length(); k++)
							quality_number.push_back(seq_information.size() * (final_quality[k] - '!'));
						fastq[detail] = quality_number;

						assert(final_string.length() == final_quality.length());
					}
				}
                                else
                                {
                                        umap[detail] += DUPCOUNT_CONSTANT + seq_information.size();
					k = read_name.rfind("=") + 1;
					if(PRINTDETAIL == PRINT_DETAIL)
						map_dupcount_info[detail] = map_dupcount_info[detail] + "," + 
							read_name.substr(k, read_name.length() - k); 

					if(FASTQOUTPUT == FASTQ_OUTPUT)		//quality score
					{
						vector<int> quality_number = fastq[detail];
						for(int k = 0; k < quality_number.size(); k++)
						{
							quality_number[k] += (seq_information.size() * (final_quality[k] - '!'));
						}
						fastq[detail] = quality_number;

						assert(quality_number.size() == final_quality.length());
					}

                                }
				
				hashcount += seq_information.size();	
			}

			//cout << read_name << endl << final_string << endl;
			//assert(false);
		}

		//CONSCOUNT = sum of exact sequences spanning across 1 or multiple BARCODE pairs
		//DUPCOUNT = total number of BARCODE groups discovered for a unique sequence 

		int split_index;
		string sequence_str = "";

		for(i = 0; i < readname_detail.size() && DUPCOUNT == 1; i++)
                {
                        detail = readname_detail[i].second;
			split_index = detail.find("#") + 1;
			sequence_str = detail.substr(split_index, detail.length()- split_index);
			//sequence_str = detail; //pRESTO Edited

			if(umap[detail] < FILTER_LOW_FREQ_READS)	
				continue;
			
			if(FASTQOUTPUT == FASTQ_OUTPUT)
				fp_fasta << "@";
			else
				fp_fasta << ">";

			fp_fasta << map_dupcount_info[readname_detail[i].second];
			fp_fasta << "|DUPCOUNT=" << umap[detail] / DUPCOUNT_CONSTANT;
			fp_fasta << "|CONSCOUNT=" << umap[detail] % DUPCOUNT_CONSTANT;
			fp_fasta << endl << sequence_str << endl;//detail;
			
			if(FASTQOUTPUT == FASTQ_OUTPUT)
			{

				fp_fasta << "+" << endl;
				string quality_score = "";
				int quality_count = (int) (umap[detail] % DUPCOUNT_CONSTANT);
				vector<int> quality_number = fastq[detail];
				
				for(int k = 0; k < quality_number.size(); k++)
				{
					int quality = (int) lrint((1.0 * quality_number[k]) / (1.0 * quality_count));
					quality_score += ('!' + quality);
				}
				fp_fasta << quality_score << endl;

			}
			
						 
                }

	}

	assert(hashcount == testcount);
        //fp_input.close();
        //fp_fasta.close();

	cout << "Total Data Reported = " << testcount << endl;
	if(DEBUG == 99)
		cout << "Total by Hash Count = " << hashcount << endl;
	return;
}


/*
This function apply poisson distribution for error correction and  to detect mutation. The input are:
1. fp_input: file stream from which we can scan the paired read alignment
2. fp_fasta: file stream to print the aligned reads with DUPCOUNT and SAM format
3. fp_detail: to print the log information in details
4. reference: reference string
5. meta_details: configurations in meta_details file
*/

void prepare_fasta_output_distribution(ifstream& fp_input, ofstream& fp_fasta, ofstream& fp_detail, 
					string& reference, meta_data& meta_details)
{
        string read_name, annotation, cigar, alignment, quality, len;
	string alignment_string, quality_score, detail;
	int total_score, average_score, opt, reflength;
	int i, k, count, refindex, readindex;

	//int fprimer_length, except_primer;
        unordered_map<string, int> umap;
	unordered_map<string, vector<int> > fastq;
	vector< pair<string, string>> sequences; 
	double exponent;
	fp_detail << "###############################################Sequence Analysis Starts##############################################" << endl;

	/*
	For each base at a position in reference, initialize the following variables of error_correction structure: 
	1. match to reference
	2. mismatch to reference
	3. expected_mismatch
	4. standard_deviation
	5. Frequency of each quality score 'qscore'
	6. frequency of each mismatching baes
	*/
	distribution error_correction[reference.length()];
	for(i = 0; i < reference.length(); i++)
	{
		for(k = 0; k < ILLUMINA_SCORE; k++)
			error_correction[i].qscore[k] = 0;

		error_correction[i].matching_base = 0;
		error_correction[i].mismatch = 0;
		error_correction[i].expected_mismatch = 0.0;
		error_correction[i].standard_deviation = 0.0;

		for(k = 0; k < 26; k++)
			error_correction[i].mismatch_base[k] = 0;
	}
	
	//fprimer_length = meta_details.fprimer.length();
	//except_primer = reference.length() - meta_details.fprimer.length();
	//fprimer_length = 0;
	//except_primer = reference.length();

	while(getline(fp_input, read_name))
        {
		//Read first four lines: read_name, annotation string, alignment_string and quality_score
                getline(fp_input, annotation);
                getline(fp_input, alignment_string);
		getline(fp_input, quality_score);
	
		//Calculate average quality_score. if it is less than GLOBAL_QSCORE then continue	
		total_score = 0;
		for(i = 0; i < quality_score.length(); i++)
		{
			total_score += quality_score.at(i) - '!';
		}

		average_score = total_score / quality_score.length();
		fp_detail << endl << "Global Average Score = " << average_score << endl;

		if(average_score < GLOBAL_QSCORE)
			continue;
	
		/*
		We assume that the reference string is given as fprimer + actual reference + rprimer
		So we compare each read base to a reference base after the fprimer_length
		*/
		//refindex = fprimer_length;
		refindex = 0;
		i = readindex = 0;
		alignment = "";
		quality = "";
		len = "";

		/*
                Rest of the code deal with the paired read alignment to a reference. In this part,
                we update the error_correction structure considering the cigar string of alignment.
                The details are given below:
                */
		int split_index = annotation.find("=") + 1;
		cigar = annotation.substr(split_index, annotation.length() - split_index);
		fp_detail << "CIGAR = " << cigar << endl;
		while(i < cigar.length())
		{
			if(cigar.at(i) >= '0' && cigar.at(i) <= '9')
				len += cigar.at(i);
			else
			{
				opt = cigar.at(i);
				count = atoi(len.c_str());
				fp_detail << "Count = " << count << ", and Option = " << opt << endl;
				len = "";			

				//fp_detail << "Reference = " << reference << endl;
				//fp_detail << "Alignment = " << alignment_string << endl;

				//For a matching or mismatching base update the variables of error_correction structure
				if(opt == 'M' || opt == 'X')
				{
					for(k = 0; k < count; k++)
					{
						if(alignment_string.at(readindex) == reference.at(refindex))
						{
							error_correction[refindex].matching_base += 1;
						}
						else 
						{
							error_correction[refindex].mismatch += 1;
							fp_detail << "BASE Mismatch is Found at = " << refindex << endl;
						}

						error_correction[refindex].qscore[quality_score.at(readindex) - '!'] += 1;
						error_correction[refindex].mismatch_base[alignment_string.at(readindex) - 'A'] += 1;//edited 11.1.15

						//fp_detail << k << "= (" << refindex << ", " << readindex << ") = " 
						//	"(" << reference.at(refindex) << ", " << alignment_string.at(readindex) << ")" << endl;
						
						//Increment refindex and readindex
						refindex += 1;
						readindex += 1;
					}	
				}
				//For insertion and deletion we just increment the readindex and refindex
				//as we cannot compare the read base to the reference base in this position
				else if(opt == 'I')
				{
					readindex += count;
				}
				else if(opt == 'D')
				{
					refindex += count;
				}
				else
				{
					assert(false);
				}
			} 
			
			i += 1;
		}

	}

	/*
	So far we have computed the step 1 as below:
	1. For a consensus base, find the frequency of A, G, C, T from all the reads aligned at that base

	Now we will continue to calculate z_score:
	2. Find the expected number of mismatch using the quality score of all bases: 
		[mu = sum ( prob ( using quality score) )]
	3. Find the standard deviation using the quality score of all bases:
		[sd = sqrt( sum ( p * (1 - p) ) )]
	4. Calculate z_score for individual mismatching base:
		[z = ( BaseFrequency + 0.5 -  mu/3 ) / sd]

	Finally we will apply the following condition in the next step:
	5. If z_score > 1.645, then do not correct; otherwise change the base to consesnus
	*/
	fp_detail << endl;
	for(i = 0; i < reference.length(); i++)
	{
		fp_detail << "Showing Analysis for Index = " << i << endl;
		if(error_correction[i].mismatch + error_correction[i].matching_base == 0)
			continue;

		//Please follow the pring string to understand following code
		for(k = 0; k < ILLUMINA_SCORE; k++)
		{
			if(error_correction[i].qscore[k] == 0)
				continue;

			exponent = pow(10, -1 * k *  0.1);
			fp_detail << "QSCORE = " << k << ", and COUNT = " << error_correction[i].qscore[k] << ", and Exponent = " << exponent << endl;
			error_correction[i].expected_mismatch += (error_correction[i].qscore[k] * exponent);
			error_correction[i].standard_deviation += (error_correction[i].qscore[k] * exponent * (1 - exponent));
			
			fp_detail << "Expected Number of Mismatch = " << error_correction[i].expected_mismatch;
	                fp_detail << ", and Standard Deviation = " << error_correction[i].standard_deviation << endl;

		}

		//Calculated overall z_score and standard_deviation
		error_correction[i].standard_deviation = sqrt(error_correction[i].standard_deviation);
		error_correction[i].zscore = (error_correction[i].mismatch + 0.5 - 
						error_correction[i].expected_mismatch) / 
						error_correction[i].standard_deviation;

		//Calculated z_score for individual bases and also expected mismatch and standard_deviation
		error_correction[i].corrected_mm = error_correction[i].expected_mismatch / 3;
		error_correction[i].corrected_sd = error_correction[i].standard_deviation / sqrt(3.0);
		error_correction[i].zscore_a = (error_correction[i].mismatch_base['A' - 'A'] + 0.5 - 
						error_correction[i].corrected_mm) / error_correction[i].corrected_sd;
		error_correction[i].zscore_c = (error_correction[i].mismatch_base['C' - 'A'] + 0.5 - 
						error_correction[i].corrected_mm) / error_correction[i].corrected_sd;
		error_correction[i].zscore_g = (error_correction[i].mismatch_base['G' - 'A'] + 0.5 - 
						error_correction[i].corrected_mm) / error_correction[i].corrected_sd;
		error_correction[i].zscore_t = (error_correction[i].mismatch_base['T' - 'A'] + 0.5 - 
						error_correction[i].corrected_mm) / error_correction[i].corrected_sd;
		//edited 11.1.15

		//fp_detail << "Error Correction At Position = " << i - fprimer_length << endl;
		fp_detail << "Error Correction At Position = " << i << endl;
		fp_detail << "Mismatch = " << error_correction[i].mismatch << ", and Matching BASE = " 
				<< error_correction[i].matching_base << endl;

		fp_detail << "Expected Number of Mismatch = " << error_correction[i].expected_mismatch;
		fp_detail << ", and Standard Deviation = " << error_correction[i].standard_deviation << endl;
		fp_detail << "Finally calculated ZSCORE = " << error_correction[i].zscore << endl;
		
		fp_detail << endl;

		//if(error_correction[i].mismatch > 0)
		//	assert(false);
	}

	int testcount = 0;

	//Apply z_score for error correction. To do this we have to read the input again to apply correction.
	fp_input.clear();
	fp_input.seekg(0, fp_input.beg);

	while(getline(fp_input, read_name))
        {
		//read four lines and apply global quality score as a filter
                getline(fp_input, annotation);
                getline(fp_input, alignment_string);
		getline(fp_input, quality_score);
		
		total_score = 0;
		for(i = 0; i < quality_score.length(); i++)
		{
			total_score += quality_score.at(i) - '!';
		}

		average_score = total_score / quality_score.length();
		fp_detail << endl << "Global Average Score = " << average_score << endl;

		if(average_score < GLOBAL_QSCORE)
			continue;
		
		//refindex = meta_details.fprimer.length();
		refindex = 0;
		i = readindex = 0;
		alignment = "";
		quality = "";
		len = "";

		//read the cigar string and apply correction for matching and mismatching bases
		int split_index = annotation.find("=") + 1;
		cigar = annotation.substr(split_index, annotation.length() - split_index);
		fp_detail << "CIGAR = " << cigar << endl;
		while(i < cigar.length())
		{
			if(cigar.at(i) >= '0' && cigar.at(i) <= '9')
				len += cigar.at(i);
			else
			{
				opt = cigar.at(i);
				count = atoi(len.c_str());
				fp_detail << "Count = " << count << ", and Option = " << opt << endl;
				len = "";			

				fp_detail << "Reference = " << reference << endl;
				fp_detail << "Alignment = " << alignment_string << endl;

				if(opt == 'M' || opt == 'X')
				{
					for(k = 0; k < count; k++)
					{
						bool zscore_flag = true;
						if(APPLY_CORRECTION == 0)// This flag is not in use anymore!!!
						{
							//In this flag we apply z_score regardless of individual base frequency at a read position
							//We apply overall z_score for error eccorection
							if(alignment_string.at(readindex) != reference.at(refindex) &&
								error_correction[refindex].zscore < ZSCORE)
							{
								//Please follow the pring string to understand the application of z_score 
								fp_detail << endl << "Error CORRECTION Point at " << readindex;
								fp_detail << "ZSCORE at this point = " << error_correction[refindex].zscore << endl;
								fp_detail << "Changing " << alignment_string.at(readindex) << " To " <<
									reference.at(readindex) << endl;
						
								alignment += reference.at(refindex);
								quality += quality_score.at(readindex);

							}
							//If z_score > 1.645, then do not correct. Might be a mutation!!!
							else
							{
								alignment += alignment_string.at(readindex);
								quality += quality_score.at(readindex);
							}
						}
						//fp_detail << k << "= (" << refindex << ", " << readindex << ")" << endl;
						else 
						{
							//In this flag we consider the individual base frequency at a read position
							//Apply individual base z_score for error correction
							if(alignment_string.at(readindex) != reference.at(refindex))
							{
								if(alignment_string.at(readindex) == 'A' && 
									error_correction[refindex].zscore_a < ZSCORE)
								{
									zscore_flag = false;
									alignment += reference.at(refindex);
									quality += quality_score.at(readindex);
								}
								else if(alignment_string.at(readindex) == 'C' && 
									error_correction[refindex].zscore_c < ZSCORE)
								{
									zscore_flag = false;
									alignment += reference.at(refindex);
									quality += quality_score.at(readindex);
								}
								else if(alignment_string.at(readindex) == 'G' && 
									error_correction[refindex].zscore_g < ZSCORE)
								{
									zscore_flag = false;
									alignment += reference.at(refindex);
									quality += quality_score.at(readindex);
								}
								else if(alignment_string.at(readindex) == 'T' && 
									error_correction[refindex].zscore_t < ZSCORE)
								{
									zscore_flag = false;
									alignment += reference.at(refindex);
									quality += quality_score.at(readindex);
								}

							}
							//No error correction is required; thus keep the base from read. Might be a mutation!!!
							if(zscore_flag == true)
							{
								alignment += alignment_string.at(readindex);
								quality += quality_score.at(readindex);
							}
						}
						refindex += 1;
						readindex += 1;
					}	
				}
				//For insertion take the base and quality from read and increment readindex
				else if(opt == 'I')
				{
					for(k = 0; k < count; k++)
					{
						alignment += alignment_string.at(readindex);
						quality += quality_score.at(readindex);

						readindex += 1;
					}
				}
				//For deletion increment the refindex
				else if(opt == 'D')
				{
					refindex += count;
				}
				else
				{
					assert(false);
				}
			} 
			
			i += 1;
		}

	
		/*
                Concatenate ciger and alignment_string and use this concatenated string as a key with the
                number of occurrence as a value which we print as DUPCOUNT.
                */
		testcount += 1;
		detail = annotation + "\n" + alignment;
		if(umap.find(detail) == umap.end())
                {
                        //umap[alignment] = 1;
                	sequences.push_back(make_pair(read_name, detail));
			umap[detail] = 1;

			if(FASTQOUTPUT == FASTQ_OUTPUT)         //quality score
			{
				vector<int> quality_number;
				for(int k = 0; k < alignment.length(); k++)
					quality_number.push_back(quality[k] - '!');
				fastq[detail] = quality_number;
			}

                }
                else
                {
                        umap[detail] += 1;

			if(FASTQOUTPUT == FASTQ_OUTPUT)         //quality score
			{
				vector<int> quality_number = fastq[detail];
				for(int k = 0; k < quality_number.size(); k++)
				{
					quality_number[k] += (quality[k] - '!');
				}
				fastq[detail] = quality_number;

				assert(quality_number.size() == quality.length());
			}

                }

		/*
                fp_fasta << ">" << read_name;

                fp_fasta << "|DUPCOUNT=1";
                fp_fasta << "|CIGAR=" << cigar << endl;
                fp_fasta << alignment << endl;
		*/
        }

	int hashcount = 0;
	int split_index;
	string sequence_str = "";

	//Create the fasta file with read_name, alignment, DUPCOUNT and CIGAR string
	for(i = 0; i < sequences.size(); i++)
	{
		hashcount += umap[sequences[i].second];
		if(umap[sequences[i].second] < FILTER_LOW_FREQ_READS)
			continue;
	
		split_index = sequences[i].second.find("=") + 1;
		sequence_str = sequences[i].second.substr(split_index, sequences[i].second.length() - split_index);

	        fp_fasta << ">" << sequences[i].first;
                fp_fasta << "|DUPCOUNT=" << umap[sequences[i].second];
                fp_fasta << "|CIGAR=" << sequence_str << endl;

                //fp_fasta << alignment << endl;

		if(FASTQOUTPUT == FASTQ_OUTPUT)
		{
			fp_fasta << "+" << endl;
			string quality_score = "";
			int quality_count = umap[sequences[i].second];
			vector<int> quality_number = fastq[sequences[i].second];
				
			for(int k = 0; k < quality_number.size(); k++)
			{
				int quality = (int) lrint((1.0 * quality_number[k]) / (1.0 * quality_count));
				quality_score += ('!' + quality);
			}
			fp_fasta << quality_score << endl;
		}

	}

        fp_input.close();
        fp_fasta.close();
	assert(hashcount == testcount);

	cout << "Total Data Reported = " << testcount << endl;
	if(DEBUG == 99)
		cout << "DUPCOUNT Total by Hash Count = " << hashcount << endl;

	return;
}




void print_vector_alignment(vector<pair<char, char> >& alignment)
{
	int i, k, index = 0;

	while(index < alignment.size())
	{
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			cout << alignment[i].first;
		}

		cout << endl;

		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
				cout << "|";
			else
				cout << " ";
		}

		cout << endl;

		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			cout << alignment[i].second;
		}

		cout << endl << endl;
		index += BREAKAT;
	}
}

void print_alignment_back(vector<pair<char, char> >& alignment_given, int ref_position, int read_position, int step)
{
	int i, k, index = 0;
	pair<char, char> swap;
	vector<pair<char, char> > alignment;	

	for(k = alignment_given.size() - 1; k >= 0; k--)
	{
		alignment.push_back(alignment_given[k]);
	}
	
	while(index < alignment.size())
	{
		cout << ref_position << "\t    ";
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			cout << alignment[i].first;
			if(alignment[i].first != '-')
				ref_position += step;
		}
		//cout << "\t" << ref_position;
		printf("%12d", ref_position - step);
		cout << endl;

		cout << "\t    ";
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
				cout << "|";
			else
				cout << " ";
		}

		cout << endl;

		cout << read_position << "\t    ";
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].second == 'o')
				cout << '-';
			else
				cout << alignment[i].second;
			if(alignment[i].second != '-')
				read_position += step;
		}
		//cout << "\t" << read_position;
		printf("%12d", read_position - step);
		cout << endl;
		cout << endl;
		index += BREAKAT;
	}
}
//comment it out from actual implementation when testing is done...
void print_alignment(vector<pair<char, char> >& alignment, string& ref, string& read, int ref_position, 
			int read_position, int step, fragment_alignment &fragment_alignment_info, bool print, ofstream& fp_detail)
{
	int i, k, index = 0;
	int updated_index;
	int fragment_count = 0;
	int start = 0, end = -1;
	int  score = 0;
	bool flag_ref = false;
	bool flag_read = false;
	bool flag_break = false;
	int ref_start, ref_end;
	int read_start, read_end;

	ref_start = ref_position;//04-01-15
	read_start = read_position;

	while(index < alignment.size())
	{
		flag_break = false;
		if(DEBUG == 99)
			fp_detail << ref_position << "\t    ";
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].second == '.')//03-26-15
			{
				updated_index = i;
				flag_break = true;
				break;
			}

			if(DEBUG == 99)	
				fp_detail << alignment[i].first;
			
			if(alignment[i].first != '-')
			{
				/*
				if(flag_ref == false)
				{
					if(alignment[i].first != '-')//03-26-15
					//if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
					{
						start = i;
						flag_ref = true;
						ref_start = ref_position;
						ref_end = ref_position;
					}
				}
				else
				{
					if(alignment[i].first != '-')//03-26-15
					//if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
                                        {
						end = i;
                                                ref_end = ref_position;
                                        }

				}
				*/
				assert(alignment[i].first == ref.at(ref_position));//03-25-15
				ref_position += step;
			}
			
		}
		//cout << "\t" << ref_position;
		
		if(DEBUG == 99)
		{
			fp_detail << std::right << std::setw(12) << (ref_position - step); 
			//printf("%12d", ref_position - step);
			fp_detail << endl;
			fp_detail << "\t    ";
		}

		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].second == '.')
				break;
			if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
			{
				if(DEBUG == 99)
					fp_detail << "|";
				score += 1;
			}
			else
				if(DEBUG == 99)
					fp_detail << " ";
		}

		if(DEBUG == 99)
		{
			fp_detail << endl;
			fp_detail << read_position << "\t    ";
		}
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].second == '.')//03-26-15
				break;
			if(alignment[i].second == 'o')
			{
				if(DEBUG == 99)
					fp_detail << '-';
			}
			else
			{
				if(DEBUG == 99)
					fp_detail << alignment[i].second;
			}
			if(alignment[i].second != '-')
			{
				/*
				if(flag_read == false)
                                {
					if(alignment[i].first != '-')//03-26-15
                                        //if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
                                        {
                                                flag_read = true;
                                                read_start = read_position;
                                                read_end = read_position;
                                        }
                                }
                                else
                                {
					if(alignment[i].first != '-')//03-26-15
                                        //if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
                                        {
                                                read_end = read_position;
                                        }

                                }
				*/
				//cout << endl << alignment[i].second << " VS " << read.at(read_position) << endl;
				assert(alignment[i].second == read.at(read_position));
				read_position += step;
			}
			
		}

		if(DEBUG == 99)
		{	
			fp_detail << std::right << std::setw(12) << (read_position - step);
			//printf("%12d", read_position - step);
			fp_detail << endl;
			fp_detail << endl;
		}
		if(flag_break == true)
		{
			break;
			index = updated_index;
			while(alignment[index].second == '.')
				index += 1;
			//while(ref_position < fragment_alignment_info.fragment_ind[fragment_count].first &&
			//	read_position < fragment_alignment_info.fragment_ind[fragment_count].second)
			fragment_count += 1;
			ref_position = fragment_alignment_info.fragment_ind[fragment_count].first;
			read_position = fragment_alignment_info.fragment_ind[fragment_count].second;
			cout << "Next ref_position = " << ref_position << ", And read_position = " << read_position << endl;
			cout << "BREAK##########################################################################DANCE" << endl << endl;
		}
		else
			index += BREAKAT;
	}

	ref_end = ref_position;//04-01-15
	read_end = read_position;
	start = 0;
	end = alignment.size() - 1;

	if(print == true)
		return;

	fragment_alignment_info.ref_start = ref_start;
        fragment_alignment_info.ref_end = ref_end - 1;//<= end

	fragment_alignment_info.read_start = read_start;
	fragment_alignment_info.read_end = read_end - 1;//<= end

	for(i = start; i <= end; i++)
	{
		fragment_alignment_info.alignment.push_back(alignment[i]);
	}
	fragment_alignment_info.end_to_end = alignment;
	fragment_alignment_info.identity_match = score;
	fragment_alignment_info.total_len = fragment_alignment_info.alignment.size();
}


