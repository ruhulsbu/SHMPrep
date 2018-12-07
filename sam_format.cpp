#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

void sam_format(fragment_alignment final_alignment_info, vector<reference_index>& refindex, 
					string& read, string& read_name, vector<string>& output, ofstream& fp_detail)
{
	if(DEBUG == 99)
	{
		fp_detail << "####################################################################################\n";
		fp_detail << "Final Result - " << ":\n";
	}

	int read_dir = 0, ref_ind = 0, ref_position = 0, maximum = 0, editdist = 0;
	int match = 0, insert = 0, delet = 0, subst = 0, ignore = 0;
	int previousop = 0, currentop = 0;
	int total_len, total_score, ref_length;
	string cigar = "", cigarseq = "";

	time_t begin, end;
	time(&begin);

	vector<pair<char, char> > alignment = final_alignment_info.alignment;

	if(!alignment.empty())
	{	
		ref_ind = final_alignment_info.ref_ind;
		read_dir = final_alignment_info.read_dir;
		total_score = final_alignment_info.identity_match;
		total_len = final_alignment_info.total_len;
		ref_position = final_alignment_info.ref_start; 
		//ref_ind = final_ref_info / MAXLEN;
		//read_dir = final_ref_info % MAXLEN;
		//total_score = final_match_info / MAXLEN;
		//total_len = final_read_info / MAXLEN;
		//ref_position = MAXLEN - final_match_info % MAXLEN;
		ref_length = (final_alignment_info.ref_end - final_alignment_info.ref_start);
	
		//print_alignment_back(alignment);

		if(DEBUG == 99)
		{
			fp_detail << "RefLen = " << refindex[ref_ind].ref.length() << ", ReadLen = " << read.length() << 
				" AlignmentLen = " << total_score << endl;//alignment.size() << endl;
		}
		memset(error_dist, 0, sizeof(long) * 10);
		error_dist[7] = total_len;
		error_dist[8] = read.length();
		error_dist[9] = ref_length;

		int total_match = 0;
		int total_insert = 0;
		int total_delete = 0;
		int total_subst = 0;
		int total_nchar = 0;
		int total_ignore = 0;
		int real_match = 0;
		//for(int i = alignment.size() - 1; i >= 0; i--)
		for(int i = 0; i < alignment.size(); i++)
		{
			if(alignment[i].first == '-' && alignment[i].second == '.')
				continue;
			//if(alignment[i].first == '-' && alignment[i].second == 'N')
			//	continue;
			if(alignment[i].first == '-' && alignment[i].second == '-')
				continue;

			if(alignment[i].second == '.')
                        {
                                ignore += 1;
                                currentop = IGNORE_;
                                total_ignore += 1;
                                //cigarseq += alignment[i].second;
                        }
			else if(alignment[i].first == alignment[i].second)
			{
				match += 1;
				real_match += 1;
				currentop = MATCH_;
				total_match += 1;
				cigarseq += alignment[i].first;
			}
			else if(alignment[i].first != alignment[i].second && alignment[i].second == '-')
			{
				//if(real_match != 0)
				//	error_free_seg[real_match] += 1;
				delet += 1;
				editdist += 1;
				currentop = DELETE_;
				total_delete += 1;
				//cigarseq += alignment[i].first;//Abscent in SAM
			}
			else if(alignment[i]. first != alignment[i].second && alignment[i].first == '-')
			{
				//if(real_match != 0)
				//	error_free_seg[real_match] += 1;
				insert += 1;
				editdist += 1;
				currentop = INSERT_;
				total_insert += 1;
				cigarseq += alignment[i].second;
			}
			else if(alignment[i].first != alignment[i].second && alignment[i].second != '.')// && alignment[i].second != '-')
			{
				
				match += 1;
				real_match += 1;
				currentop = MATCH_;
				total_match += 1;
				cigarseq += alignment[i].second;
					
			}
			/*
			else if(alignment[i].first == '.' || alignment[i].second == '.')
			{
				ignore += 1;
				currentop = IGNORE_;
				total_ignore += 1;
				cigarseq += alignment[i].second;
			}
			*/
			if(i == 0) previousop = currentop;

			if(currentop != previousop)
			{
				ostringstream numstr;
				if(previousop == MATCH_)
				{
					numstr << match;
					cigar += numstr.str() + "M";
					error_free_seg[real_match] += 1;
					//error_dist[read.length()][2] += (match - real_match);
					error_dist[0] = max(error_dist[0], real_match);
					match = 0;
					real_match = 0;
				}
				else if(previousop == DELETE_)
				{
					numstr << delet;
					cigar += numstr.str() + "D";
					//error_dist[read.length()][1] += delet;
					delet = 0;
				}
				else if(previousop == INSERT_)
				{
					numstr << insert;
					cigar += numstr.str() + "I";
					//error_dist[read.length()][0] += insert;
					insert = 0;
				}
				else if(previousop == SUBSTITUTE_)
				{
					numstr << subst;
					cigar += numstr.str() + "X";
					subst = 0; 
				}
				else if(previousop == IGNORE_)
				{
					numstr << ignore;
					cigar += numstr.str() + "N";
					ignore = 0;
				}

				//cout << "Updating previos op = " << previousop << " with CIGAR = " << cigar << endl;
				previousop = currentop;

			}

			//cout << alignment[i].first << " and " << alignment[i].second << " where op = " << currentop << endl;
		}
		/*
		if(alignment.size() != 0)
		{
			print_alignment_back(alignment);
		}
		*/
		ostringstream numstr;
		if(previousop == MATCH_)
		{
			numstr << match;
			cigar += numstr.str() + "M";
			error_free_seg[real_match] += 1;
			//error_dist[read.length()][2] += (match - real_match);
			error_dist[0] = max(error_dist[0], real_match);
			match = 0;
			real_match = 0;
		}
		else if(previousop == DELETE_)
		{
			numstr << delet;
			cigar += numstr.str() + "D";
			//error_dist[read.length()][1] += delet;
			delet = 0;
		}
		else if(previousop == INSERT_)
		{
			numstr << insert;
			cigar += numstr.str() + "I";
			//error_dist[read.length()][0] += insert;
			insert = 0;
		}
		else if(previousop == SUBSTITUTE_)
		{
			numstr << subst;
			cigar += numstr.str() + "X";
			subst = 0;
		}
		else if(previousop == IGNORE_)
		{
			numstr << ignore;
			cigar += numstr.str() + "N";
			ignore = 0;
		}


		error_dist[1] = total_match;
		error_dist[2] = total_insert;
		error_dist[3] = total_delete;
		error_dist[4] = total_subst;
		error_dist[5] = total_nchar;
		error_dist[6] = total_ignore;

		if(DEBUG == 99)
		{
			fp_detail << "total_match = " << total_match << endl;
			fp_detail << "total_insert = " << total_insert << endl;
			fp_detail << "total_delete = " << total_delete << endl;
			fp_detail << "total_subst = " << total_subst << endl;
			fp_detail << "total_nchar = " << total_nchar << endl;
			fp_detail << "total_ignore = " << total_ignore << endl;
			fp_detail << "total reference = " << ref_length << endl;
		}
		//assert(cigarseq.length() == (alignment.size() - total_delete));
		assert(total_match + total_insert + total_delete + total_subst + total_nchar + total_ignore == 
					total_insert + ref_length);

		if(DEBUG == 99)
		{
			fp_detail << "Chain length = " << alignment.size()  << ", and Edit Distance: " << editdist <<  ", In Direction: " << 
				read_dir  << " at ref_positin = " << ref_position << ", while CIGAR = " << cigar.length() << endl;
		}
	}
	else
	{
		total_score = alignment.size();//final_match_info / MAXLEN;
		double seconds = difftime(end, begin);
		if(DEBUG == 99)
		{
			fp_detail << "Total Time Taken inside sam_format = " << seconds << " seconds" << endl;
			fp_detail << "####################################################################################\n";

			fp_detail << endl << endl;
		}
		return;	
	}

	output.push_back(read_name);	//1. read_name
	if(total_score == 0)		//2. direction
		output.push_back("4");
	else if(read_dir == 1)
		output.push_back("0");
	else
		output.push_back("16");
	output.push_back(refindex[ref_ind].name);	//3. reference name
	ostringstream numstr;		//4. index
	numstr << (ref_position + 1);//03-09-2015
	output.push_back(numstr.str());
	output.push_back("255");	//5. default mapq previously *
	if(cigar.length() == 0)		//6. cigar
		output.push_back("*");	
	else
		output.push_back(cigar);
	output.push_back("*");		//7. rnext
	output.push_back("0");		//8. pnext
	output.push_back("0");		//9. tlen 
	if(cigarseq.length() != 0)	//10. seq as cigar
		output.push_back(cigarseq);
	else 
		output.push_back("*");
	//output.push_back("*");		//11. default quality
	output.push_back(final_alignment_info.quality);

	time(&end);
	double seconds = difftime(end, begin);
	
	if(DEBUG == 99)
	{
		fp_detail << "Total Time Taken inside sam_format = " << seconds << " seconds" << endl;
		fp_detail << "####################################################################################\n";
	
		fp_detail << endl << endl;
	}
	return;
}

