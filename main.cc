#define NREAD 1000000  //numbers of reads of the file, this splits the whole input file into N read chunks
//If we divide enough such that each segment can fit in computer cache, computating may speed up.
#define NBINS 1000 //number of bins, such as country index #
#define PIXEL_SIZE 231.65635825  //pixel size in m  (assuming x_size = y_size)

#define INPUT_INT
#define TYPE_BYTE  //Make sure NBINS and type values play nice with each other

#ifdef INPUT_INT
	#define IN_BYT_PER_PIX 2
	#define in_t short
#endif

#ifdef TYPE_BYTE
	#define INDEX_BYT_PER_PIX 1
	#define index_t char
#endif

#ifdef TYPE_INT
	#define INDEX_BYT_PER_PIX 2
	#define index_t short
#endif


#include <iostream>
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
	if (!S_ISREG(index_file_stat.st_mode)) {
		std::cout << "Error: opening input file: " << argv[2] << '\n';
	}
	if (S_ISREG(out_file_stat.st_mode)) {
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

	unsigned long long npix = in_file_stat.st_size / IN_BYT_PER_PIX;
	unsigned long long cache_pix = npix / NREAD;
	unsigned long remainder_pix = npix % NREAD;
	in_t in_cache[cache_pix];
	index_t index_cache[cache_pix];

	//READ and process each block
	std::ifstream in_file(argv[1], std::ios::binary|std::ios::in);
	std::ifstream index_file(argv[2], std::ios::binary|std::ios::in);
	if (in_file.is_open() && index_file.is_open()) {	
		for(int i=0; i<NREAD; i++) {
			for(int j=0; j<NBINS; j++) {
				temp_count[j] = 0;
				temp_total[j] = 0;
			}
			in_file.read((char *) &in_cache, cache_pix*IN_BYT_PER_PIX);
			index_file.read((char *) &index_cache, cache_pix*INDEX_BYT_PER_PIX);
			for(int j=0; j<cache_pix; j++) {
				if((index_cache[j] >= 0) && (index_cache[j] < NBINS)) {   //ignore index that fall outside (i.e. -1 for water)
					temp_total[index_cache[j]] += in_cache[j];
					temp_count[index_cache[j]]++;
				}
			}
			//Dump temp values into final total variables
			for(int j=0; j<NBINS; j++) {
				if(temp_count[j]) {
					total[j] += temp_total[j];
					count[j] += temp_count[j];
				}
			}
		}

		//check for remainder
		if(remainder_pix) {
			in_file.read((char *) &in_cache, remainder_pix*IN_BYT_PER_PIX);
			index_file.read((char *) &index_cache, remainder_pix*INDEX_BYT_PER_PIX);
			for(int j=0; j<NBINS; j++) {
				temp_count[j] = 0;
				temp_total[j] = 0;
			}
			
			for(int j=0; j<remainder_pix; j++) {
				if((index_cache[j] >= 0) && (index_cache[j] < NBINS)) {   //ignore index that fall outside (i.e. -1 for water)
					temp_total[index_cache[j]] += in_cache[j];
					temp_count[index_cache[j]]++;
				}
			}
			//Dump temp values into final total variables
			for(int j=0; j<NBINS; j++) {
				if(temp_count[j]) {
					total[j] += temp_total[j];
					count[j] += temp_count[j];
				}
			}
			

		}


	} else {
		std::cout << "Error opening input files" << '\n';
		return 1;
	}

	//close input files
	in_file.close();
	index_file.close();

	//Write stats output
	std::ofstream out_file;
	out_file.open(argv[3]);
	out_file << "Country Code, Total, NPixels, Total AGB (Gg)\n";
	double pixel_area = PIXEL_SIZE * PIXEL_SIZE;  //pixel area in m^2
	for (int i=0; i<NBINS; i++) {
		if (count[i] > 0) {
			out_file << i << ',' << total[i]  << ',' << count[i] << ',' << total[i] / (pixel_area / 10000) / 1000 << '\n';
		}
	}

	out_file.close();

	return 0;
}
