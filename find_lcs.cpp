#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

/*
Suffix Automaton of a string “S" in simple terms is a directed acyclic graph where 
vertices or nodes are called “states” and the arcs or the edges between these nodes 
is called the “transition” between these states.

One of the states(nodes) is denoted by  “Initial State (t_o)”  from where we can reach to 
all other states in the suffix automaton.

One or more of this states are marked as  “Terminal States” such that if we  follow 
the link from any node to the next until reaching any of the terminal states 
and note down the labels of the edges then we get a suffix of the original string “S”.

The algorithm for construction of suffix automaton is online i.e we construct 
by adding a single character to our previous string, modifying the previous structure. 

Each state will consist of  len (length of  longest suffix), link (suffix link of 
the state), and transition edges to the next states. We also keep the string index of 
the character for which this node is created, and flag to indicate terminal state.

Initially our structure consist of a single state (t_o), which we shall assume to be 
zeroth state (all other states will get numbers 1,2,3, …). Assign this state len = 0 
and link = -1 ( meaning non-realistic link). Accordingly the whole task now is to add 
a single character “ch” to the end of the current line.

Let last  is the state(node) that corresponds to the entire current line before 
adding the symbol (initially last = 0) and after adding each new character, we will 
change the value to (1, 2, 3, … and so on). After adding a new character we will make 
a new node(state) ”cur“, we will assign len(cur) = len(last) + 1.

Till this time we  have created a new state, initialized it but we haven’t added 
it to our structure. For doing so we will run a loop, initially we are at the “last” 
node of tree, if there is no edge/transition with label “ch” we will add an edge 
between “last” and “cur” with the label “ch” and move to the node pointed by the suffix 
link of  “last” by changing last to suffix_link(last), we will keep on doing this until 
we reach (t_o) or we encounter the condition mentioned in next point.

If at some node the transition with label “ch” already exists, we will stop at that node, 
let that state be represented by “pnode”,  let the state which is connected to “pnode” 
via transition with label “ch” be represented by “qnode”. Now we have two cases depending 
upon if len(pnode) + 1 = len(qnode) or not –

if len(pnode) + 1 = len(qnode), we can simply assign link(“cur”)  = “qnode” and return. 
Otherwise we will have to create a clone of state qnode with everything remaining same as 
state qnode except for the value of len such that len(clone) = len(pnode) +  1 and assign 
link(“cur”)  = “clone”  and break. If it never occurred then we reach dummy state t_o, 
in this case we simply assigned link(“cur”) = 0 and return.
*/


//Initialize Suffix Automata of Size n defined 
void initialize_automata(int n, automata_state *automata)
{
	int i, j;

	for (i = 0; i < n; i++)
	{
		automata[i].index = -1;
		automata[i].link = -1;
		automata[i].len = 0;
		automata[i].terminal = false;
		automata[i].child.clear();	
	}
}

//Extend the Automata for a character in a sequence
//index is the position of ch in a sequence
//last is the index of the Automata node that will be extended
//sz is the size of the Suffix Automata created so far
void extend_automata(char ch, int index, int& last, int& sz, automata_state *automata)
{
	//cur is the current node in the Automata to hold ch
	//sz is increased by 1 after adding a new node with ch
	//automata[cur].index is updated with the index of ch in reference
	//automata[cur].len is updated with the length of the 
	//path from initial node to last node increased by 1
	int cur = sz++;
	automata[cur].index = index;
	automata[cur].len = automata[last].len + 1;

	//cout << "Processing index = " << cur << ", Character = " << ch << endl;
	
	int pnode;
	//Make a transition from all the possible pnodes to the cur node
	//break when reach the intiial state t_o or there exists a state
	//from which we created a transition to another state with ch label
	for(pnode = last; pnode != -1 && automata[pnode].child[ch] == 0; 
					pnode = automata[pnode].link)
	{
		automata[pnode].child[ch] = cur;
	}

	if(pnode == -1)
	{
		//We reached the initial state t_o, updating the cur link
		automata[cur].link = 0;
	}
	else
	{
		int qnode;
		qnode = automata[pnode].child[ch];

		//If the distance between pnode and qnode is one
		//then create a link from cur node to qnode
		if(automata[pnode].len + 1 == automata[qnode].len)
		{
			automata[cur].link = qnode;
		}
		else
		{
			//Otherwise create a clone node of the qnode
			//Create transions from the pnode and all it's link
			//for which there were a transition labeled with ch
			//Break if there exists no transion or t_o is reached 
			int clone = sz++;
			automata[clone].index = automata[qnode].index;
			automata[clone].link = automata[qnode].link;
			automata[clone].len = automata[pnode].len + 1;
			automata[clone].child = automata[qnode].child;

			for( ; pnode != -1 && automata[pnode].child[ch] == qnode; 
						pnode = automata[pnode].link)
			{
				automata[pnode].child[ch] = clone;
			}

			//update the link of qnode and cur node to the clone node
			automata[qnode].link = clone;
			automata[cur].link = clone;
		}
	}

	//update the last node to the cur node
	last = cur;
	
}


//Find the Longest Common Substring between query and reference using the Suffix Automata
//Return the reference index, query index, and length of LCS
int longest_common_substring(string& query, int& refindex, int& queryindex, int& max_len, 
				automata_state *automata, ofstream& fp_detail)
{
	int i, len, pnode, query_len;

	query_len = query.length();
	max_len = 0;
	pnode = 0;
	len = 0;

	//fp_detail << "Query = " << query << endl;
	//fp_detail << "Length = " << query.length() << endl;

	//For each character "ch" of the query, traverse the Automata
	for(i = 0; i < query_len; i++)
	{
		//printf("current pnode = %d\n", pnode);

		//If there exists no transition from a state of Automata, then move
		//to the next longest suffix following the link of last pnode and try to
		//find a match at different index of the reference string. Break if t_o
		//node is reached or if there exists a trnsition labeled with current ch
		for( ; pnode != -1 && automata[pnode].child[query.at(i)] == 0; )
						//pnode = automata[pnode].link)
		{
			//update the length as we follow the suffix link of pnode
			pnode = automata[pnode].link;
			//printf("go to parent = %d\n", pnode);
			len = automata[pnode].len;
			//printf("query[%d] = %c with len = %d\n", i, query[i], len);
		}

		//If t_o node is reached the length of current substring using current ch is 0
		if(pnode == -1)
		{
			pnode = 0;
			len = 0;
		}
		else
		{
			//Otherwise, update the pnode and length following the transitions with ch
			pnode = automata[pnode].child[query.at(i)];
			//printf("pnode updated = %c\n", query[i]);
			len++;
		}

		//Keep record of maximum length, reference index and query index
		if(max_len < len)
		{
			max_len = len;
			//if(max_len > 40)
			//	fp_detail << "max_len updated = " << max_len << endl;
			refindex = automata[pnode].index - max_len + 1;
			queryindex = i - max_len + 1;
		}
	}
	/*	
	if(refindex != -1 && queryindex != -1)
	{
		fp_detail << endl << "Reference Index = " << refindex << endl;
		fp_detail << "Query Index = " << queryindex << endl;
		fp_detail << "Substring = " << query.substr(queryindex, max_len) << endl << endl;
	}
	*/
	return max_len;
}

string find_primer(automata_state automata[], string& query, int& refindex, int& queryindex, ofstream &fp_detail)
{
	refindex = queryindex = -1;
	int max_len = -1;
	//find longes common substring
	longest_common_substring(query, refindex, queryindex, max_len, automata, fp_detail);

	fp_detail << "Max length of primer string in LCS = " << max_len << endl;
	
	/*
	if(query.length() == max_len)
		return true;
	else
		return false;
	*/

	if(max_len < 1)
		return "";
	else
		return query.substr(queryindex, max_len);
}

void create_automata(string& reference, automata_state automata[])
{
	int sz = 1;
	int last = 0;

	//initialize the automata
	initialize_automata(STATE_SIZE, automata);

	//extend the automata
	for(int i = 0; i < reference.length(); i++)
        {
                extend_automata(reference.at(i), i, last, sz, automata);
        }
		
}

int find_lcs(string& reference, string& query, int& refindex, int& queryindex, int& max_len, 
			automata_state automata[], int& sz, int& last, bool extend, ofstream& fp_detail)
{
	//If there exists a Suffix Automata then we do not create again
	if(extend == true)
		initialize_automata(STATE_SIZE, automata);

	//cout << "Reference = " << strref << endl;
	//cout << "Length = " << strref.length() << endl;

	//Assuming that there exists no Suffix Automata for a reference 
	//sequence and we have a Suffix Automata with an initial node t_o.
	//For each character of a reference sequence we extend that Automata
	//In the beginning the last node of Automata is the initial node,
	//hence for the first itertion last node is the 0th index of 
	//the automata array and size sz of the Automata is 1
	for(int i = 0; i < reference.length() && extend; i++)
	{
		extend_automata(reference.at(i), i, last, sz, automata);
	}

	//Find the Longest Common Substring between a reference and query
	//using Suffix Automata built on the reference sequence
	longest_common_substring(query, refindex, queryindex, max_len, automata, fp_detail);

	//fp_detail << endl << "Reference Index = " << refindex << endl;
	//fp_detail << "Query Index = " << queryindex << endl << endl;
	//fp_detail << "Substring = " << query.substr(queryindex, max_len) << endl;
	//fp_detail << "Longest Common Substrings: " << max_len << endl;
	
	return 0;
}
