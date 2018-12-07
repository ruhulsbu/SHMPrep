//constant variables

#ifndef CONSTANT_H
#define CONSTANT_H

//#define INFINITY (1 << 28)
#define SIZE 500
#define STATE_SIZE (4 * SIZE)

#define WINDOW 6	//size of window//should be 6
#define ANCHOR_WORD 3
#define KBAND 10	//size of band for string comparisons
#define HBAND 200	//currently unused

#define MAXLEN 10000000LL	//used to store x, y as (x * MAXLEN + y) as long long
#define MINLEN 100		//used to store x, y as (x * MINLEN + y) as long
#define MAXTRY 10		//number of trials to select reference for a target

#define GAPX 9		//delete or insert	//positive	6 default value 999 extreme value
#define MAT_OPEN 2	//start matching bases
#define GAP_OPEN 2	//start of gap
#define WEIGHT 5	//match bases
#define MISMATCH -4	//substitute		//negative	-4 default value
#define DEBUG 99	//debug code 99

#define UP 9		//symbol for up arrow
#define DIAG 6		//symbol for diagonal arrow
#define BACK 3		//symbol for back arrow
#define BREAKAT 60	//number for printing the string alignment

#define MATCH_ 1
#define INSERT_ 2
#define DELETE_ 3
#define SUBSTITUTE_ 4
#define NCHAR_ 5
#define IGNORE_ 6

#define FF 1		//forward forward alignment
#define FR 2		//forward reverse alignment
#define RF 3		//reverse forward alignment
#define RR 4		//reverse reverse alignment

#define HBAND_THRESHOLD 0.0	
#define KBAND_THRESHOLD 0.0

#define LOCAL 1
#define OPTIMIZE 0
#define FOREACHDIR 0 

#define HRATIO 0.95
#define ERROR_ 0.5
//#define KMER 11
//#define SEED 45 
#define FULLREAD 1
#define BASE 4

#define MAXREADLEN 50000
#define MINREADLEN 50

#define FRAGMENT_SIZE 400
#define GAP_PERCENT_MATCH 55
#define KBAND_PERCENT_MATCH 55
#define GLOBAL_PERCENT_MATCH 85.0f
#define OVERLAP_PERCENT_MATCH 50.0f
#define DEFAULT_PRIMER_QUALITY 90.0f
#define LCS_PRIMER_SIMILARITY 80.0f
#define CHAIN_PERCENT_MATCH 55
#define ALIGNMENT_COVERAGE 1.15f

#define PRIMER_KMER_SIZE 5
#define PRIMER_ID_CONSTANT 100000
#define DUPCOUNT_CONSTANT 10000000000LL

//#define FORWARD_ADAPTER_LENGTH 0
//#define REVERSE_ADAPTER_LENGTH 0
//#define FORWARD_BARCODE_LENGTH 15
//#define REVERSE_BARCODE_LENGTH 0

//#define GLOBAL_QSCORE 25
//#define PERBASE_QSCORE 30
#define ILLUMINA_SCORE 50
#define ZSCORE 1.645f
#define FASTQ_OUTPUT 1
#define STRIP_OFF 1
#define PRINT_DETAIL 1

//#define BASIC_CORRECTION 0
#define MAX_THREAD 16
#define APPLY_BAYES 1
#define APPLY_CORRECTION 1
#define TRUNCATE_REFERENCELESS_ALIGNMENT 1
#define CONSERVATIVE_REFERENCELESS_ALIGNMENT 0
#endif
