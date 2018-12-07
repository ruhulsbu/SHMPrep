#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

//cell **matrix = NULL;

int find_kband_similarity(string& str1, string& str2, vector<pair<char, char> >& alignment, int &str1_end, int &str2_end, cell **matrix)
{
        int kband, match, preoff = -1, shift;
        int i, j, h, k, row, column, offset;
        float len_ratio;
        string swap;
	/*
        if(str1.length() < str2.length())
        {
                swap = str1;
                str1 = str2;
                str2 = swap;
        }
	*/

	kband = GLOBAL_KBAND;
        row = str1.length() + 1;
        column = 2 * kband + 2;

        len_ratio = 1.0 * str2.length() / str1.length();
     	
        for(i = 0; i < row; i ++)
        {
                matrix[i][0].cost = i * -GAP * LOCAL;// * 0;
                matrix[i][0].dir = UP;
		matrix[i][0].matrix_col = 0;
		matrix[i][0].str2_index = 0;
		matrix[i][0].match = 0;
		matrix[i][0].length = 0;
        }

        for(j = 1; j < column; j++)
        {
                matrix[0][j].cost = j * -GAP * LOCAL;
                matrix[0][j].dir = BACK;
		matrix[0][j].matrix_col = j - 1;
		matrix[0][j].str2_index = j;
		matrix[0][j].match = 0;
		matrix[0][j].length = 0;
        }
       
        int max_score = 0.0;
        int max_row, max_col;
	int max_length, max_match;
        max_row = max_col = 1;
	int gap_cost = 0;
	str1_end = str2_end = 0;
	max_length = max_match = 0;
	int start_col, end_col;
	int affinity = 0;

        for(i = 1; i < row; i++)
        {
                offset = (int) (len_ratio * i);
		start_col = max(1 - offset, -kband);
		end_col = min(str2.length() - offset, kband);
		//printf("i = %d, k_start = %d, k_end = %d\n", i, (offset + start_col), (offset + end_col));
                
                for(h = start_col; h <= end_col; h++)
                {
                        k = offset + h;

                        if(k >= 1 && k <= str2.length())
                        {
				
				if(str1.at(i - 1) == str2.at(k - 1))
				{
					match = WEIGHT;
					
					if(i > 1 && k > 1)//03-26-15
					{
						if(str1.at(i - 2) == str2.at(k - 2))
						{
							match = WEIGHT + MAT_OPEN;
						}
					}
				}
				else
				{
					match = MISMATCH;
				
					if(i > 1 && k > 1)
					{
						if(str1.at(i - 2) == str2.at(k - 2))
							match = MISMATCH - MAT_OPEN;
					}
				} 
				
                                if(offset > kband + 1)
                                        j = k - (offset - kband) + 1;
                                else
                                        j = k;

                                if(preoff == offset || offset <= kband + 1)
                                        shift = -1;
                                else
                                        shift = offset - preoff -1;

                                //assert(i >= 1 && i < row);
                                //assert(j >= 1 && j < 2 * kband + 2);


                                matrix[i][j].dir = DIAG;
				matrix[i][j].matrix_col = (j + shift);
				matrix[i][j].str2_index = k;
				matrix[i][j].length = matrix[i - 1][j + shift].length + 1;

				if(match >= WEIGHT)
					matrix[i][j].match = matrix[i - 1][j + shift].match + 1;
				else
					matrix[i][j].match = matrix[i - 1][j + shift].match;

                                if(LOCAL == 0)
                                        matrix[i][j].cost = max(0, matrix[i - 1][j + shift].cost + match);
                                else
                                        matrix[i][j].cost = matrix[i - 1][j + shift].cost + match;

                                if(j + shift + 1 <= 2 * kband + 1)
                                {
					if(matrix[i - 1][j + shift + 1].dir == DIAG)
                                   		gap_cost = GAP + GAP_OPEN;
                                	else
                                        	gap_cost = GAP;

                                        if(matrix[i][j].cost < matrix[i - 1][j + shift + 1].cost - gap_cost) 
					{ // <= changed to <
				                matrix[i][j].dir = UP;
						matrix[i][j].matrix_col = (j + shift + 1);
						matrix[i][j].str2_index = k;
                                        	matrix[i][j].cost = matrix[i - 1][j + shift + 1].cost - gap_cost;
						matrix[i][j].length = matrix[i - 1][j + shift + 1].length + 1;
						matrix[i][j].match = matrix[i - 1][j + shift + 1].match;
					}
                                }

                                if(j - 1 >= 1)// && str2.at(k - 1) != 'N')
                                {
					if(matrix[i][j - 1].dir == DIAG)
                                        	gap_cost = GAP + GAP_OPEN;
                                	else
                                        	gap_cost = GAP;

                                        if(matrix[i][j].cost < matrix[i][j - 1].cost - gap_cost)
					{ // <= changed to <
                                                matrix[i][j].dir = BACK;
						matrix[i][j].matrix_col = (j - 1);
						matrix[i][j].str2_index = k;
                                        	matrix[i][j].cost = matrix[i][j - 1].cost - gap_cost;
						matrix[i][j].length = matrix[i][j - 1].length + 1;
						matrix[i][j].match = matrix[i][j - 1].match;
					}
                                }

                                if(DEBUG == 1)
                                {
					//cout << "@(" << i << "," << j << "," << k << "): " << matrix[i][j] << "\t";
                                        //cout << "@(" << str1.at(i - 1) << ", " << str2.at(k - 1) << "): " << matrix[i][j] << "\t";
                                        cout << "(" << str1.at(i - 1) << "," << str2.at(k - 1) << "):"
                                                << "(" << i << "," << j << "): R = " << k << ", C = " << matrix[i][j].matrix_col 
							<< " at Direction: " << matrix[i][j].dir << "\n";
                                }

                               	max_score = matrix[i][j].cost;
                              	max_row = i;
                              	max_col = j;
                               	str1_end = i;
                              	str2_end = k;
				max_length = matrix[i][j].length;
				max_match = matrix[i][j].match;
				//cout << "Updating str1_end = " << i << ", and str2_end = " << str2_end <<
				//      ", with ratio = " << (100.0 * k / i) << endl;

                        }

                }

		//if(max_match + row - i < 0.51 * (max_length + row - i))
		//	break;

                preoff = offset; //to prevent updating same row values
		
		if(j + 1 < column)
		{
                        matrix[i][j + 1].cost = -(i + j + 1) * GAP * LOCAL;
			matrix[i][j + 1].dir = BACK;
			matrix[i][j + 1].matrix_col = j - 1;
			matrix[i][j + 1].str2_index = j;
			matrix[i][j + 1].match = 0;
			matrix[i][j + 1].length = 0;
		}
		
                if(DEBUG == 1) {
                        cout << "\n";
                        cout << "-------------------------------------------------------------------------------" << endl;
                }
		
		//printf("\n");

        }

        if(DEBUG == 1)
                cout << "$(" << str1.length() << "," << (kband + 1) << ") Score=" << matrix[str1.length()][kband + 1].cost << endl;
        
        //no use now//print_path(path, row - 1, column / 2, str1, str2, alignment);//delete
        {
        	//print_path_back(path, row - 1, column / 2, str1, str2, alignment);
        	//cout << "Total score = " << matrix[max_row][max_col] << ", with kband_score = " << max_score << endl;
        	//cout << "Max score found at the row = " << str1_end << ", and at the col = " << str2_end << endl;
                print_path_cell_back(matrix, max_row, max_col, str1, str2, alignment);
        }

        match = matrix[str1.length()][kband + 1].cost;

        if(DEBUG == 1)
	{
                cout << "Printing Path Completed" << endl;
		printf("The KBAND matching score = %d\n", match);
		printf("The line number is (end of kband: %d\n", __LINE__);
	}	

	return match;
}


int find_local_similarity(string& str1, string& str2, vector<pair<char, char> >& alignment, int& str1_start, int& str2_start, 
						int& str1_end, int& str2_end, cell **matrix)
{
        int i, k, row, column;

	if(str1.length() > FRAGMENT_SIZE || str2.length() > FRAGMENT_SIZE)
	{
		cout << "Change definition of FRAGMENT_SIZE in constant.h" << endl;
		exit(0);
	}
        
	row = str1.length() + 1;
        column = str2.length() + 1;
	
        for(i = 0; i < row; i ++)
        {
                matrix[i][0].cost = 0;
                matrix[i][0].dir = UP;
		matrix[i][0].matrix_col = 0;
		matrix[i][0].str2_index = 0;
		matrix[i][0].match = 0;
		matrix[i][0].length = 0;
        }

        for(k = 1; k < column; k++)
        {
                matrix[0][k].cost = 0;
                matrix[0][k].dir = BACK;
		matrix[0][k].matrix_col = 0;
		matrix[0][k].str2_index = 0;
		matrix[0][k].match = 0;
		matrix[0][k].length = 0;
        }

        int max_score = 0.0;
        int max_row, max_col;
	int max_length, max_match;
        max_row = max_col = 1;
	int match, gap_cost = 0;
	//int str1_end, star2_end;
	str1_start = str2_start = 0;
	str1_end = str2_end = 0;
	max_length = max_match = 0;

        for(i = 1; i < row; i++)
        {
                for(k = 1; k < column; k++)
                
		{
			/*
			if(str1.at(i - 1) == str2.at(k - 1))
			{
				match = WEIGHT;
					
				if(i > 1 && k > 1)
				{
					if(str1.at(i - 2) == str2.at(k - 2))
					{
						match = WEIGHT + MAT_OPEN;
					}
				}
			}
			else
			{
				match = MISMATCH;
			
				if(i > 1 && k > 1)
				{
					if(str1.at(i - 2) == str2.at(k - 2))
						match = MISMATCH - MAT_OPEN;
				}
			}
			*/

			if(str1[i - 1] == str2[k - 1])
				match = WEIGHT;
			else
				match = MISMATCH; 
				
     
                        matrix[i][k].dir = DIAG;
			matrix[i][k].matrix_col = (k - 1);
			matrix[i][k].str2_index = k;
			matrix[i][k].length = matrix[i - 1][k - 1].length + 1;

			if(match >= WEIGHT)
				matrix[i][k].match = matrix[i - 1][k - 1].match + 1;
			else
				matrix[i][k].match = matrix[i - 1][k - 1].match;

                        matrix[i][k].cost = max(0, matrix[i - 1][k - 1].cost + match);

			if(matrix[i - 1][k].dir == DIAG)
                        	gap_cost = GAP + GAP_OPEN;
                        else
                        	gap_cost = GAP;

                        if(matrix[i][k].cost < matrix[i - 1][k].cost - gap_cost) 
			{ // <= changed to <
				matrix[i][k].dir = UP;
				matrix[i][k].matrix_col = k;
				matrix[i][k].str2_index = k;
                                matrix[i][k].cost = matrix[i - 1][k].cost - gap_cost;
				matrix[i][k].length = matrix[i - 1][k].length + 1;
				matrix[i][k].match = matrix[i - 1][k].match;
			}
                        
			if(matrix[i][k - 1].dir == DIAG)
                        	gap_cost = GAP + GAP_OPEN;
                        else
                         	gap_cost = GAP;

                        if(matrix[i][k].cost < matrix[i][k - 1].cost - gap_cost)
			{ // <= changed to <
                        	matrix[i][k].dir = BACK;
				matrix[i][k].matrix_col = (k - 1);
				matrix[i][k].str2_index = k;
                                matrix[i][k].cost = matrix[i][k - 1].cost - gap_cost;
				matrix[i][k].length = matrix[i][k - 1].length + 1;
				matrix[i][k].match = matrix[i][k - 1].match;
			}
                        
			if(DEBUG == 1)
			{
				//cout << "@(" << i << "," << j << "," << k << "): " << matrix[i][j] << "\t";
                                //cout << "@(" << str1.at(i - 1) << ", " << str2.at(k - 1) << "): " << matrix[i][j] << "\t";
                                cout << "(" << str1.at(i - 1) << "," << str2.at(k - 1) << "):"
                	                << "(" << i << "," << k << "): R = " << k << ", C = " << matrix[i][k].matrix_col 
					<< " at Direction: " << matrix[i][k].dir << "\n";
			}

			if(max_score < matrix[i][k].cost)
			{
                      		max_score = matrix[i][k].cost;
                     		max_row = i;
                     		max_col = k;
                  		str1_end = i;
                    		str2_end = k;
				max_length = matrix[i][k].length;
				max_match = matrix[i][k].match;
				//cout << "Updating str1_end = " << i << ", and str2_end = " << str2_end <<
				//      ", with ratio = " << (100.0 * k / i) << endl;
			}

                }
	
                if(DEBUG == 1) {
                        cout << "\n";
                        cout << "-------------------------------------------------------------------------------" << endl;
                }
		
		//printf("\n");

        }

        if(DEBUG == 1)
                cout << "$(" << max_row << "," << max_col << ") Score=" << matrix[max_row][max_col].cost << endl;
        
        //no use now//print_path(path, row - 1, column / 2, str1, str2, alignment);//delete
        {
        	//print_path_back(path, row - 1, column / 2, str1, str2, alignment);
        	//cout << "Total score = " << matrix[max_row][max_col] << ", with kband_score = " << max_score << endl;
        	//cout << "Max score found at the row = " << str1_end << ", and at the col = " << str2_end << endl;
                print_local_alignment_back(matrix, max_row, max_col, str1, str2, str1_start, str2_start, alignment);
        }

	if(DEBUG == 1)
	{
                cout << "Printing Path Completed" << endl;
		printf("The KBAND matching score = %d\n", match);
		printf("The line number is (end of kband: %d\n", __LINE__);
	}

	return max_match;
	
}




int find_banded_similarity(string& str1, string& str2, vector<pair<char, char> >& alignment, int &end1, int &end2)
{
	//int FRAGMENT_SIZE = 400;

	float ratio;
	vector<pair<char, char> > sub_alignment;
	vector<pair<char, char> > current_alignment;

	int alignment_quality = 0;
	int str1_len, str2_len;
	int str1_start, str2_start;	
	int str1_end, str2_end;
	int rest_of_str1, rest_of_str2;

	string str1_substr, str2_substr;

	str1_len = str1.length();
	str2_len = str2.length();
	str1_start = str2_start = 0;
	cout << endl;

	int total_score = 0, total_len = 0;
	bool flag = false;

	while(str1_start < str1_len && str2_start < str2_len)
	{
		rest_of_str1 = str1.length() - str1_start;
		rest_of_str2 = str2.length() - str2_start;

		//cout << "rest_of_str1 = " << rest_of_str1 << ", and rest_of_str2 = " << rest_of_str2 << endl;

		ratio = 1.0 * rest_of_str2 / rest_of_str1;
		/*
		if(ratio < 0.65 || ratio > 1.35)
		{
			cout << "str1 is too bigger than str2 but continue" << endl;
		//	break;
		}
		*/
		if(rest_of_str1 > 2 * FRAGMENT_SIZE && rest_of_str2 > 2 * FRAGMENT_SIZE)
			rest_of_str1 = rest_of_str2 = FRAGMENT_SIZE;
		else
		{
			rest_of_str1 = rest_of_str2 = min(rest_of_str1, rest_of_str2);
			flag = true;
		}

		str1_substr = str1.substr(str1_start, rest_of_str1);
		str2_substr = str2.substr(str2_start, rest_of_str2);
		
		//cout << "str1 = " << str1_substr.substr(0, min(100, str1_substr.length())) << endl;
		//cout << "str2 = " << str2_substr.substr(0, min(100, str2_substr.length())) << endl;
		sub_alignment.clear();
		//find_similarity(str1_substr, str2_substr, sub_alignment, str1_end, str2_end);
		//find_kband_similarity(str1_substr, str2_substr, sub_alignment, str1_end, str2_end);
	
		int sub_alignment_score = 0;	
		for(int i = sub_alignment.size() - 1, preoccur = 0; i >= 0; i--)
		{
			if(sub_alignment[i].first == sub_alignment[i].second && sub_alignment[i].first != '-')
			{
				preoccur += 1;
				if(preoccur == 1)
					sub_alignment_score += WEIGHT;
				else
					sub_alignment_score += (WEIGHT + MAT_OPEN);
			}
			else
				preoccur = 0;
		}
		if(1.0 * sub_alignment_score / sub_alignment.size() < 3.40)
			break;
		/*
		print_alignment_back(sub_alignment, str1_substr.length(), str2_substr.length(), -1);
		*/
		if(str1_end < KMER || str2_end < KMER) // added_on 03-15-15
			break;

		for(int i = sub_alignment.size() - 1, preoccur = 0; i >= 0; i--)
		{
			current_alignment.push_back(sub_alignment[i]);
			if(sub_alignment[i].first == sub_alignment[i].second && sub_alignment[i].first != '-')
				total_score += 1;	
		}
		total_len += sub_alignment.size();
		alignment_quality += sub_alignment_score;

		str1_start += str1_end;
		str2_start += str2_end;
		/*
		cout << "str1_end = " << str1_end << ", and next str1_start = " << str1_start << endl;
		cout << "str2_end = " << str2_end << ", and next str2_start = " << str2_start << endl;
		cout << "total_score = " << total_score << ", and total_len = " << total_len << endl;
		*/
		if(100.00 * total_score / total_len < 1.0 * KBAND_PERCENT_MATCH)//fix the percentage here
		{
			cout << "breaking at the matching ratio = " << (100.00 * total_score / total_len) << endl;
			break;
		}
		if(str1_end * 4 < FRAGMENT_SIZE && str2_end * 4 < FRAGMENT_SIZE || flag == true)//Should be 2 for 400
		{
			cout << "breaking at the fragment_range problem" << endl;
			break;
		}

		ratio = 1.0 * str2_start / str1_start;
		if(ratio < 0.7 || ratio > 1.3)
		{
 			cout << "str1 and str2 ratio does not maintain the expected measurement" << endl;
 			break;
 		}

		//if(1.0 * alignment_quality / total_len < 3.40)
		//	exit(1);
	}

	int matching_score = 0;
	int ref_start = 0;
	int read_start = 0;
	for(int i = current_alignment.size() - 1; i >= 0; i--)
	{
		if(current_alignment[i].first != '-')
			ref_start += 1;
		if(current_alignment[i].second != '-')
			read_start += 1;
		if(current_alignment[i].first == current_alignment[i].second)
		{
			matching_score += 1;
		}
		alignment.push_back(current_alignment[i]);
	}
	current_alignment.clear();

	end1 = ref_start;
	end2 = read_start;
	cout << "alignment agrees to the index number (" << end1 << ", " << end2 << ")" << endl;

	return matching_score;

	alignment.clear();
	for(int i = current_alignment.size() - 1; i >= 0; i--)
		alignment.push_back(current_alignment[i]);

	ref_start = 0;
	read_start = 0;

	for(int i = 0; i < alignment.size(); i++) 
	{
		if(alignment[i].first != '-')
		{
			ref_start += 1;
		}
		if(alignment[i].second != '-')
		{
			read_start += 1;
		}
	}
	
	end1 = str1_start;
	end2 = str2_start;

	assert(str1_start == ref_start && str2_start == read_start);
	cout << "alignment agrees to the index number" << endl;
	
	return total_score;
}

int find_levenshtein_distance(string& s, string& t)
{
    // degenerate cases
    //if (s == t) return 0;
    //if (s.Length == 0) return t.Length;
    //if (t.Length == 0) return s.Length;

    // create two work vectors of integer distances
    int v0[t.length() + 1];
    int v1[t.length() + 1];
    int cost;

    // initialize v0 (the previous row of distances)
    // this row is A[0][i]: edit distance for an empty s
    // the distance is just the number of characters to delete from t
    for (int i = 0; i < t.length() + 1; i++)
        v0[i] = i;

    for (int i = 0; i < s.length(); i++)
    {
        // calculate v1 (current row distances) from the previous row v0

        // first element of v1 is A[i+1][0]
        //   edit distance is delete (i+1) chars from s to match empty t
        v1[0] = i + 1;

        // use formula to fill in the rest of the row
        for (int j = 0; j < t.length(); j++)
        {
            cost = (s[i] == t[j]) ? 0 : 1;
            v1[j + 1] = min(v1[j] + 1, min(v0[j + 1] + 1, v0[j] + cost));
        }

        // copy v1 (current row) to v0 (previous row) for next iteration
        for (int j = 0; j < t.length() + 1; j++)
            v0[j] = v1[j];
    }

    return v1[t.length()];
}


//http://www.geeksforgeeks.org/dynamic-programming-set-4-longest-common-subsequence/
int find_longest_common_subsequence(string& str1, string& str2)
{
	int m, n, i, j;
	m = str1.length();
	n = str2.length();

	int matrix[2][n + 1];
	//memset(matrix, 0, sizeof(matrix[0][0]) * 2 * (n + 1));
	
	/* Following steps build L[m+1][n+1] in bottom up fashion. Note 
	that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1] */
	
	for (i = 0; i <= m; i++)
	{
		for (j = 0; j <= n; j++)
		{
			if (i == 0 || j == 0)
				matrix[1 - i % 2][j] = 0;
  
			else if (str1[i-1] == str2[j-1])
				matrix[1 - i % 2][j] = matrix[i % 2][j-1] + 1;
  
			else
				matrix[1 - i % 2][j] = max(matrix[i % 2][j], matrix[1 - i % 2][j - 1]);
		}
	}
    
	/* L[m][n] contains length of LCS for X[0..n-1] and Y[0..m-1] */
	return matrix[1 - m % 2][n];
}

