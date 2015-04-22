#include <iostream>
#define NREAD 1000000  //numbers of reads of the file, this splits the whole input file into N read chunks
//If we divide enough such that each segment can fit in computer cache, computating may speed up.
#define NBINS 1000 //number of bins, such as country index #

#define INPUT_INT
#define TYPE_BYTE

#ifdef INPUT_INT
	#define IN_BYT_PER_PIX 2
	#define in_t short
#endif

#ifdef TYPE_BYTE
	#define INDEX_BYT_PER_PIX 1
	#define index_t char
#endif


#include <fstream>
#include <sys/stat.h>


int main(int argc, char* argv[])
{
	if (argc != 4) {
		std::cout << "Syntax:  calc_stats <input file> <index file> <out file> " << argc << '\n';
		return 1;
	}

	//Test for existance of input files, and make sure output file does not exist
	struct stat in_file_stat, index_file_stat, out_file_stat;
	stat(argv[1], &in_file_stat);
	stat(argv[2], &index_file_stat);
	stat(argv[3], &out_file_stat);

	if (!S_ISREG(in_file_stat.st_mode)) {
		std::cout << "Error: opening input file: " << argv[1] << '\n';
	}
	if (!S_ISREG(in_index_stat.st_mode)) {
		std::cout << "Error: opening input file: " << argv[2] << '\n';
	}
	if (S_ISREG(out_File_stat.st_mode)) {
		std::cout << "Error: output file already exist " << argv[3] << '\n';
	}

	//check file size matching of input files
	if ((in_file_stat.st_size / IN_BYT_PER_PIX) != (index_file_stat.st_size / INDEX_BYT_PER_PIX)) {
		std::cout << "Error: input file sizes do not match for given data types" << '\n';
		return 1;
	}

	//Make arrays used for summation
	long long count[NBINS], temp_count[NBINS];
	double total[NBINS], temp_total[NBINS];

	


	return 0;
}
