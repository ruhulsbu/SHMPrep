#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

/*
Apply Referenceless Alignment. The input are:
1. readseq_first: first read string
2. readseq_second: second read string
3. quality_first: first read quality string
4. quality_second: second read quality string

return the results in:
5. readsequence: overlapped read string
6. readquality: overlapped read quality string
7. matrix: to calculated the similarity matrix
8. fp_detail: log the output in details
*/
int referenceless_alignment(string& readseq_first, string& readseq_second, string& quality_first, string& quality_second, 
				string& readsequence, string& readquality, cell **matrix, ofstream& fp_detail)
{
	fp_detail << "First Step: Starting Overlapping the Read Pairs" << endl;
	fp_detail << "-----------------------------------------------" << endl;

	int sz = 1, last = 0, matching_base;
	automata_state automata[STATE_SIZE];
	int rind = -1, qind = -1, len;
	int refindex = -1, queryindex = -1, max_len;

	//comment it when testing is done
	assert(readseq_first.length() == quality_first.length());
	assert(readseq_second.length() == quality_second.length());

	if(readseq_first.length() != quality_first.length() ||
		readseq_second.length() != quality_second.length())
			return -1;

	//Apply LCS to choose the overlap direction of the two reads: forward-forward or forward-reverse_complement
	find_lcs(readseq_first, readseq_second, refindex, queryindex,  max_len, automata, sz, last, true, fp_detail);	
	if(refindex == -1 || queryindex == -1)
		return -1;


	string query = reverse_complement(readseq_second);
	find_lcs(readseq_first, query, rind, qind, len, automata, sz, last, false, fp_detail);
	if(rind == -1 || qind == -1)
		return -1;

	//According to the overlap direction choose the second read and quality either in forward or reverse_complement
	if(max_len <= len)
	{
		max_len = len;
		refindex = rind;
		queryindex = qind;
		reverse_str(quality_second);
	}
	else
	{
		query = readseq_second;
	}
	
	/*
	fp_detail << "First ) " << endl;
	fp_detail << "Read String = " << readseq_first << endl;//.substr(0, 80) << endl;
        fp_detail << "Quality =  :: " << quality_first << endl << endl;

        fp_detail << "Second) " << endl;
        fp_detail << "Read String = " << query << endl;//.substr(0, 80) << endl;
        fp_detail << "Quality =  :: " << quality_second << endl << endl;
	*/
	assert(query.length() == quality_second.length());

	fp_detail << "Reference Index = " << refindex << endl;
	fp_detail << "Query Index = " << queryindex << endl;
	fp_detail << "Substring = " << query.substr(queryindex, max_len) << endl;
	fp_detail << "Longest Common Substring Length: " << max_len << endl << endl;

	int subtract = min(refindex, queryindex);
	int addition = min(readseq_first.length() - refindex, query.length() - queryindex);
	
	//fp_detail << "subtract from refindex = " << subtract << endl;
	//fp_detail << "addition to refindex = " << addition << endl;

	int ref_start_index = refindex - subtract;
	int ref_end_index = refindex + addition;
	string reference = readseq_first.substr(ref_start_index, ref_end_index - ref_start_index);

	int read_start_index = queryindex - subtract;
	int read_end_index = queryindex + addition;
	string read = query.substr(read_start_index, read_end_index - read_start_index);

	fp_detail << "Second Step: Overlapped Region:" << endl;
	fp_detail << "-------------------------------" << endl;
	fp_detail << "Read1::: " << reference << endl;
	fp_detail << "Quality: " << quality_first.substr(ref_start_index, ref_end_index - ref_start_index) << endl;
	fp_detail << "Read2::: " << read << endl;
        fp_detail << "Quality: " << quality_second.substr(read_start_index, read_end_index - read_start_index) << endl;

	vector<pair<char, char> > alignment;
	bool alignment_flag = true;
	matching_base = 0;
	for(int k = 0; k < reference.length(); k++)
	{
		if(reference[k] == read[k])
		{
			matching_base += 1;
		}
		alignment.push_back(make_pair(reference[k], read[k]));

	}

	fp_detail << "Overlap Quality: " << (100.0 * matching_base / reference.length()) << endl << endl;
		
	if(LOCALFLAG == LOCAL && 100.0 * matching_base / reference.length() < OVERLAP_QUALITY) 
	//	|| reference[0] != read[0] || reference[reference.length() - 1] != read[read.length() -1])
	{
		int str1_start, str2_start, str1_end, str2_end;
		alignment_flag = false;
		alignment.clear();

		//fp_detail << "Current Percent of Identity: " << (100.0 * matching_base / reference.length()) << endl;
		find_local_similarity(reference, read, alignment, str1_start, str2_start, str1_end, str2_end, matrix);//return -1;
	
		if(DEBUG == 99)
		{
			fp_detail << "Local Alignment: " << endl;
			for(int k = 0; k < alignment.size(); k++)
				fp_detail << alignment[k].first;
			fp_detail << endl;

			matching_base = 0;
			for(int k = 0; k < alignment.size(); k++)
			{
				fp_detail << alignment[k].second;
				if(alignment[k].first == alignment[k].second)
					matching_base += 1;
			}
			fp_detail << endl << endl;

			fp_detail << "Identity Match After Local Alignment: " << (100.0 * matching_base / reference.length()) << endl; 
			fp_detail << "Str1 Start = " << str1_start << endl;
			fp_detail << "Str2 Start = " << str2_start << endl;
			fp_detail << "Str1 End   = " << str1_end << endl;
			fp_detail << "Str2 End   = " << str2_end << endl;

		}
	}
	
	///////////////////////////////////////////////////////////////////////////////////

	//Apply Bayes theorem for choosing the correct base and quality for overlapped read bases
        fp_detail << "Result:: ";

	string merged_string = "";
	string merged_quality = "";

	char base_quality;
	char x_base, y_base;
	char x_quality, y_quality;
	int x_index = 0, y_index = 0;
	//edit here to jeewoen's alignment rule
	//for(int i = ref_start_index, k = read_start_index; i < ref_end_index && k < read_end_index; i++, k++)
	for(int i = 0; i < alignment.size(); i++)
	{
		if(alignment_flag == true)
		{
			x_base = alignment[i].first;
			y_base = alignment[i].second;
			x_quality = quality_first[i + ref_start_index];
			y_quality = quality_second[i + read_start_index];
		}
		else
		{
			if(alignment[i].first == '-')
			{
				x_base = 'N';
				x_quality = '!';
			}
			else
			{
				x_base = alignment[i].first;
				x_quality = quality_first[x_index + ref_start_index];
				x_index += 1;
			}

			if(alignment[i].second == '-')
			{
				y_base = 'N';
				y_quality = '!';
			}
			else
			{
				y_base = alignment[i].second;
				y_quality = quality_second[y_index + read_start_index];
				y_index += 1;
			}
		}
		

		//check for base quality; should not be necessary
		if(x_quality < '!' || x_quality > 'J' ||
			y_quality < '!' || y_quality > 'J')
			return -1;

		//Choose the base that is not 'N'
		if(x_base == 'N' || y_base == 'N')
		{
			if(x_base != y_base)
			{
				if(x_base != 'N')
				{
					merged_string += x_base;
					merged_quality += x_quality;
				}
				else
				{
					merged_string += y_base;
					merged_quality += y_quality;
				}
			}
			else
			{
				//if(HARD_RULE == true)
				return -1;
			}
		}
		//If both the bases are same then apply pre-caculated base quality using the quality score from two reads
		else if(x_base == y_base)
		{
			merged_string += x_base;

			if(APPLY_BAYES == 1)
			{
                                fp_detail << "-";
				base_quality = calculate_basequal(x_quality, y_quality, 1);

				if(base_quality < '!' or base_quality > 'J')
					assert(false);

				merged_quality += base_quality;	
				continue;
			}

			//Basic rule is applied if APPLY_BAYES is not one. Choose the base with higher quality
			if(x_quality > y_quality)
				merged_quality += x_quality;
			else
				merged_quality += y_quality;
		}		
		else
		{
			//If the bases do not match then choose the base with maximum quality and 
			//apply pre-calculated base quality using the quality score from two reads
			if(x_quality > y_quality)
			{
                                fp_detail << ">";
				merged_string += x_base;
				if(APPLY_BAYES == 1)
				{
					base_quality = calculate_basequal(x_quality, y_quality, 0);
					merged_quality += base_quality;

					if(base_quality < '!' or base_quality > 'J')
						assert(false);
				}
				else	//Basic rule is applied if APPLY_BAYES is not one.
					merged_quality += x_quality;
			}
			else if(x_quality == y_quality)
			{
                                fp_detail << "-";
                                if(x_base > y_base)
					merged_string += x_base;
				else
					merged_string += y_base;
				merged_quality += x_quality;
				//If two bases do not match but have same quality score then return -1.
				//if(HARD_RULE == true)
				//cout << "Shameful Case!!!";
				//exit(0);//for testing this function
				//return -1;//original code
			}
			else//This does the opposite of the if condition above.
			{
                                fp_detail << "<";
				merged_string += y_base;
				if(APPLY_BAYES == 1)
				{
					base_quality = calculate_basequal(x_quality, y_quality, 0);
					merged_quality += base_quality;

					if(base_quality < '!' or base_quality > 'J')
						assert(false);
				}
				else	//Basic rule is applied if APPLY_BAYES is not one.				
					merged_quality += y_quality;
			}
		}
	}

        fp_detail << endl;
	assert(merged_string.length() == merged_quality.length());

	fp_detail << "Merged:: " << merged_string << endl;
	fp_detail << "Quality: " << merged_quality << endl << endl;



	//Create overlapped read and it's quality string
	string str_first_part = "", str_second_part = "";
	string quality_first_part = "", quality_second_part = "";
	
	if(ref_start_index >= read_start_index)
	{
		str_first_part = readseq_first.substr(0, ref_start_index);
        	quality_first_part = quality_first.substr(0, ref_start_index);

        	str_second_part = query.substr(read_end_index, query.length() - read_end_index);
        	quality_second_part = quality_second.substr(read_end_index, query.length() - read_end_index);
	}
        else
        {
                str_first_part = query.substr(0, read_start_index);
                quality_first_part = quality_second.substr(0, read_start_index);

                str_second_part = readseq_first.substr(ref_end_index, readseq_first.length() - ref_end_index);
                quality_second_part = quality_first.substr(ref_end_index, readseq_first.length() - ref_end_index);
        }

	/*
	fp_detail << "CheckQuality >" << endl;
	
	fp_detail << "First  string:: " << str_first_part << endl;
	fp_detail << "First  quality: " << quality_first_part << endl << endl;

	fp_detail << "Second string:: " << str_second_part << endl;
	fp_detail << "Second quality: " << quality_second_part << endl << endl;
	*/
	fp_detail << "Total Alignment Length (Considering Substitution Only) = " 
			<< (str_first_part.length() + merged_string.length() + str_second_part.length()) << endl << endl;

	fp_detail << "Third Step: Create the Overlapped Alignment" << endl;
	fp_detail << "-------------------------------------------" << endl;
	string string_alignment = str_first_part + merged_string + str_second_part;
	string string_quality = quality_first_part + merged_quality + quality_second_part;
	fp_detail << "Alignment: " << string_alignment << endl;
	fp_detail << "Quality::: " << string_quality << endl << endl;


	assert(str_first_part.length() == quality_first_part.length());
	assert(str_second_part.length() == quality_second_part.length());
	assert(string_alignment.length() == string_quality.length());

	readsequence = string_alignment;
	readquality = string_quality;
	/*
	if(referenceless_flag == true && REFERENCELESS_ALIGNMENT == 1)//quickfix
	{
		readsequence = reverse_complement(string_alignment);
		reverse_str(readquality);
	}
	*/
	//if(IGNOREN == 0)
	{
		fp_detail << "Safe Return from Optimize KBAND" << endl << endl;
		return 0;
	}

	for(int i = 0; i < str_first_part.length(); i++)
	{
		if(str_first_part.at(i) == 'N')
			return -1;

		if(quality_first_part.at(i)< '!' || quality_first_part.at(i) > 'J')
                        return -1;

	}


	for(int i = 0; i < str_second_part.length(); i++)
	{
		if(str_second_part.at(i) == 'N')
			return -1;

		if(quality_second_part.at(i)< '!' || quality_second_part.at(i) > 'J')
                        return -1;

	}

	if(DEBUG == 99)
		fp_detail << "Safe Return from Optimize KBAND" << endl;
	return 0;
}



int create_gap_alignment(vector<pair<char, char> >& alignment)
{
	int mismatch = 0;
	int score = 0;
	int SEEDS = 100;
	int i;

	for(i = 0; i < alignment.size(); i++)
        {
                //cout << alignment[i].first << ", " << alignment[i].second << endl;

                if(alignment[i].first == alignment[i].second)
                {
                        if(alignment[i].first != '-')
                                score += 1;
                }
	}

	cout << "unoptimized length = " << alignment.size() << ", and score = " << score << endl;

	if(100.0 * score / alignment.size() > 1.0 * GAP_PERCENT_MATCH)
		return alignment.size();

	for(i = 0; i < alignment.size(); i++)
	{
		if(alignment[i].first == alignment[i].second)
		{
			if(alignment[i].first != '-')
			{
				if(100.0 * score / (alignment.size() - i - 1) > 1.0 * GAP_PERCENT_MATCH)
					break;

				score = score - 1;
			}

		}
	}	

	if(i > 0 && i <= alignment.size())
	{
                //alignment.erase(alignment.begin(), alignment.begin() + i);
		for(int k = 0; k < i; k++)
		{
			if(alignment[k].second != '-')
				alignment[k] = make_pair(alignment[k].first, 'o');
			else
				alignment[k] = make_pair(alignment[k].first, '-');
		}
	}
	
	cout << "alignment length = " << alignment.size() << ", optimized length = " << (alignment.size() - i) << 
			", and score = " << score << endl;

	return (alignment.size() - i);
}

int optimize_path(vector<pair<char, char> >& alignment)
{
	int match = 0, i;
	for(i = 0; i < alignment.size(); i++)
	{
		//cout << alignment[i].first << ", " << alignment[i].second << endl;

		if(alignment[i].first == alignment[i].second) 
		{
			if(alignment[i].first != '-')
				break;
			else
				match = match + MISMATCH;
		}
		else
			match = match - GAP;

	}

	if(i > 0 && i < alignment.size())
		alignment.erase(alignment.begin(), alignment.begin() + i);
	
	for(i = alignment.size() - 1; i >= 0; i--)
	{
		//cout << alignment[i].first << ", " << alignment[i].second << endl;

		if(alignment[i].first == alignment[i].second) 
		{
			if(alignment[i].first != '-')
				break;
			else
				match = match + MISMATCH;
		}
		else
			match = match - GAP;

	}

	if(i + 1 > 0 && i + 1 < alignment.size())
		alignment.erase(alignment.begin() + i + 1, alignment.end());

	for(i = 0; i < alignment.size(); i++)
	{
		//cout << alignment[i].first << ", " << alignment[i].second << endl;

		if(alignment[i].first == alignment[i].second) 
		{
			if(alignment[i].first != '-')
				match = match + WEIGHT;
			else
				match = match + MISMATCH;
		}
		else
			match = match - GAP;

	}


	return match;
}

void print_path_matrix(long long **path, int row, int column)
{
	int i, j;
	cout << "ROW = " << row << " and COL = " << column << endl;
	if(path == NULL)
	{
		cout << "NULL Pointer Detected" << endl;
		return;
	}
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < column; j++)
		{
			cout << "(" << i << "," << j << "): " << path[i][j] << "\t";
		}
		cout << endl;
	}

}

void print_path_cell_back(cell **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment)
{
	int dir, k;
	char x, y;

	while(row > 0 || column > 0)
	{	
		dir = (int) path[row][column].dir;
		k = path[row][column].str2_index;
		//cout << "Row = " << row << ", Column = " << column << endl;
		//cout << "K = " << k << " VS dir = " << path[row][column].dir << endl;
		if(dir == DIAG)
		{
			//assert(k > 0 && row > 0);
			x = str1.at(row - 1);
			y = str2.at(k - 1);
			//cout << x << ", " << y << " at " << k << endl;
		}
		else if(dir == UP)
		{
			//assert(row > 0);
			x = str1.at(row - 1);
			y = '-';
		}
		else
		{
			//assert(k > 0);
			x = '-';
			y = str2.at(k - 1);
			//cout << x << ", " << y  << " at " << k << endl;
		}

		//if(!(x == '-' && y == 'N'))
		alignment.push_back(make_pair(x, y));


		column = path[row][column].matrix_col;

		if(dir == DIAG || dir == UP) 
			row = row - 1;	
	}

}


void print_local_alignment_back(cell **path, int row_index, int column_index, string& str1, string& str2, 
				int& str1_start, int& str2_start, vector<pair<char, char> >& final_alignment)
{
	int dir;
	char x, y;
	int row = row_index;
	int column = column_index;
	vector<pair<char, char> > alignment;

	//while(row > 0 || column > 0)
	while(path[row][column].cost > 0)
	{	
		dir = (int) path[row][column].dir;
		//cout << "Row = " << row << ", Column = " << column << endl;
		//cout << " VS dir = " << path[row][column].dir << endl;
		if(dir == DIAG)
		{
			//assert(k > 0 && row > 0);
			x = str1.at(row - 1);
			y = str2.at(column - 1);
			//cout << x << ", " << y << " at " << k << endl;
		}
		else if(dir == UP)
		{
			//assert(row > 0);
			x = str1.at(row - 1);
			y = '-';
		}
		else
		{
			//assert(k > 0);
			x = '-';
			y = str2.at(column - 1);
			//cout << x << ", " << y  << " at " << k << endl;
		}

		//if(!(x == '-' && y == 'N'))
		alignment.push_back(make_pair(x, y));
		str1_start = row - 1;
		str2_start = column - 1;

		column = path[row][column].matrix_col;

		if(dir == DIAG || dir == UP) 
			row = row - 1;	
	}

	int i, k, min_index;
	if(row > column)
	{
		min_index = row - column;
		for(i = 0; i < min_index; i++)
		{
			final_alignment.push_back(make_pair(str1[i], '-'));
		}
		k = 0;
	}
	else if(row == column)
	{
		i = k = min_index = 0;
	}
	else
	{
		min_index = column - row;
		for(k = 0; k < min_index; k++)
		{
			final_alignment.push_back(make_pair('-', str2[k]));
		}
		i = 0;
	}

	min_index = min(row, column);
	for(int j = 0; j < min_index; j++)
	{
		final_alignment.push_back(make_pair(str1[i + j], str2[k + j]));
	}

	for(int j = alignment.size() - 1; j >= 0; j--)
	{
		final_alignment.push_back(alignment[j]);
	}

	min_index = min(str1.length() - row_index, str2.length() - column_index);
	for(int j = 0; j < min_index; j++)
	{
		final_alignment.push_back(make_pair(str1[row_index + j], str2[column_index + j]));
	}

	if(row_index < column_index)
	{
		i = row_index + min_index;
		for( ; i < str1.length(); i++)
		{
			final_alignment.push_back(make_pair(str1[i], '-'));
		}
	}
	
	if(row_index > column_index)
	{
		k = column_index + min_index;
		for( ; k < str2.length(); k++)
		{
			final_alignment.push_back(make_pair('-', str2[k]));
		}
	}


	if(DEBUG == -1)
	{	
		cout << "Testing Local Alignment: " << endl;
		for(int k = alignment.size() - 1; k >= 0; k--)
        		cout << alignment[k].first;
	        cout << endl;
        	for(int k = alignment.size() - 1; k >= 0; k--)
	        	cout << alignment[k].second;
        	cout << endl << endl;
	}

}


void print_path_back(long long **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment)
{
	int dir, k;
	char x, y;

	while(row != 0 || column != 0)
	{	
		dir = (int) ((path[row][column] / MAXLEN) / MAXLEN);
		k = path[row][column] % MAXLEN;

		if(dir == DIAG)
		{
			x = str1.at(row - 1);
			y = str2.at(k - 1);
		}
		else if(dir == UP)
		{
			x = str1.at(row - 1);
			y = '-';
		}
		else
		{
			x = '-';
			y = str2.at(k - 1);
		}

		//if(!(x == '-' && y == 'N'))
		alignment.push_back(make_pair(x, y));


		column = (path[row][column] / MAXLEN) % MAXLEN;

		if(dir == DIAG || dir == UP) 
			row = row - 1;	
	}

}

void print_path(long long **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment)
{
	int i, j, k, dir;
	char x, y;
	
	i = row;
	k = path[row][column] % MAXLEN;
	j = (path[row][column] / MAXLEN) % MAXLEN;
	dir = (int) ((path[row][column] / MAXLEN) / MAXLEN);


	assert(i >= 0);
	assert(k >= 0);
	assert(j >= 0);
	assert(dir >= 0);
	assert(i <= str1.length());
	assert(k <= str2.length());
	assert(j <= 2 * KBAND + 2);
	
	if(DEBUG == 1) cout << row << "," << column << ": " << path[row][column] << endl;
	if(row == 0 && column == 0)
		return;
	
	if(DEBUG == 1)
		cout << "i = " << i << ", k = " << k << ", j = " << j << endl;

	//if(i == 0 || k == 0)
	//	return;

	if(dir == DIAG)
		print_path(path, i - 1, j, str1, str2, alignment);
	else if(dir == UP)
		print_path(path, i - 1, j, str1, str2, alignment);
	else
		print_path(path, i, j, str1, str2, alignment);
	
	if(dir == DIAG)
	{
		//if(str1.at(i - 1) == str2.at(k - 1))
		//	x = y = str1.at(i - 1);

		x = str1.at(i - 1);
		y = str2.at(k - 1);
	}
	else if(dir == UP)
	{
		x = str1.at(i - 1);
		y = '-';
	}
	else
	{
		x = '-';
		y = str2.at(k - 1);
	}

	alignment.push_back(make_pair(x, y));
}


