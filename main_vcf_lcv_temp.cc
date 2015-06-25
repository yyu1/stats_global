//This version takes in a VCF file
//Also includes LandCover map to decide which carbon fraction to use

#define NREAD 10000  //numbers of reads of the file, this splits the whole input file into N read chunks
//If we divide enough such that each segment can fit in computer cache, computating may speed up.
#define NBINS 255 //number of bins, such as country index #
#define PIXEL_SIZE 231.65635825  //pixel size in m  (assuming x_size = y_size)

#define VCF_THRESH2 0
#define VCF_THRESH3 5
#define VCF_THRESH4 10

#define INPUT_INT
#define TYPE_BYTE  //Make sure NBINS and type values play nice with each other

#ifdef INPUT_INT
	#define IN_BYT_PER_PIX 2
	#define in_t short
#endif

#ifdef TYPE_BYTE
	#define INDEX_BYT_PER_PIX 1
	#define index_t unsigned char
#endif

#ifdef TYPE_INT
	#define INDEX_BYT_PER_PIX 2
	#define index_t short
#endif


#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include "sum_vcf_lcv_temp.c"

int main(int argc, char* argv[])
{
	if (argc != 6) {
		std::cout << "Syntax:  calc_stats_vcf <input file> <index file> <vcf file> <lcv file> <out file> " << argc << '\n';
		return 1;
	}

	//Test for existance of input files, and make sure output file does not exist
	struct stat in_file_stat, index_file_stat, vcf_file_stat, lcv_file_stat, out_file_stat;
	stat(argv[1], &in_file_stat);
	stat(argv[2], &index_file_stat);
	stat(argv[3], &vcf_file_stat);
	stat(argv[4], &lcv_file_stat);
	stat(argv[5], &out_file_stat);

	if (!S_ISREG(in_file_stat.st_mode)) {
		std::cout << "Error: opening input file: " << argv[1] << '\n';
		return 1;
	}
	if (!S_ISREG(index_file_stat.st_mode)) {
		std::cout << "Error: opening input file: " << argv[2] << '\n';
		return 1;
	}
	if (!S_ISREG(vcf_file_stat.st_mode)) {
		std::cout << "Error: opening input file: " << argv[3] << '\n';
		return 1;
	}
	if (!S_ISREG(lcv_file_stat.st_mode)) {
		std::cout << "Error: opening input file: " << argv[4] << '\n';
		return 1;
	}
	if (S_ISREG(out_file_stat.st_mode)) {
		std::cout << "Error: output file already exist " << argv[5] << '\n';
		return 1;
	}

	//check file size matching of input files
	if ((in_file_stat.st_size / IN_BYT_PER_PIX) != (index_file_stat.st_size / INDEX_BYT_PER_PIX) || ((index_file_stat.st_size / INDEX_BYT_PER_PIX) != vcf_file_stat.st_size) || ((index_file_stat.st_size / INDEX_BYT_PER_PIX) != lcv_file_stat.st_size) )
	{
		std::cout << "Error: input file sizes do not match for given data types" << '\n';
		return 1;
	}

	//Make arrays used for summation
	long long count[NBINS], temp_count[NBINS];
	long long count2[NBINS], temp_count2[NBINS];
	long long count3[NBINS], temp_count3[NBINS];
	long long count4[NBINS], temp_count4[NBINS];
	double total[NBINS], temp_total[NBINS];
	double total2[NBINS], temp_total2[NBINS];
	double total3[NBINS], temp_total3[NBINS];
	double total4[NBINS], temp_total4[NBINS];
	double global_total, temp_global_total;
	double global_total2, temp_global_total2;
	double global_total3, temp_global_total3;
	double global_total4, temp_global_total4;
	long long global_count, temp_global_count;
	long long global_count2, temp_global_count2;
	long long global_count3, temp_global_count3;
	long long global_count4, temp_global_count4;

	for (int i=0; i<NBINS; i++) {
		count[i] = 0;
		total[i] = 0;
		count2[i] = 0;
		total2[i] = 0;
		count3[i] = 0;
		total3[i] = 0;
		count4[i] = 0;
		total4[i] = 0;
	}

	global_total = 0;
	temp_global_total = 0;
	global_count = 0;
	temp_global_count = 0;
	global_total2 = 0;
	temp_global_total2 = 0;
	global_count2 = 0;
	temp_global_count2 = 0;
	global_total3 = 0;
	temp_global_total3 = 0;
	global_count3 = 0;
	temp_global_count3 = 0;
	global_total4 = 0;
	temp_global_total4 = 0;
	global_count4 = 0;
	temp_global_count4 = 0;

	unsigned long long npix = in_file_stat.st_size / IN_BYT_PER_PIX;
	unsigned long long cache_pix = npix / NREAD;
	unsigned long remainder_pix = npix % NREAD;
	in_t in_cache[cache_pix];
	index_t index_cache[cache_pix];
	unsigned char vcf_cache[cache_pix];
	unsigned char lcv_cache[cache_pix];
	std::cout << "Total pix: " << npix << " Cache pix: " << cache_pix << " Remainder: " << remainder_pix << '\n';
	std::cout.flush();

	//READ and process each block
	std::ifstream in_file(argv[1], std::ios::binary|std::ios::in);
	std::ifstream index_file(argv[2], std::ios::binary|std::ios::in);
	std::ifstream vcf_file(argv[3], std::ios::binary|std::ios::in);
	std::ifstream lcv_file(argv[4], std::ios::binary|std::ios::in);
	std::cout << "Processing input file...";
	std::cout.flush();
	global_count_block tmp_global_count;

	if (in_file.is_open() && index_file.is_open() && vcf_file.is_open()) {

		for(int i=0; i<NREAD; i++) {
			if ((i % (NREAD/100)) == 0) {
				std::cout << i / (NREAD/100) << ' ';
				std::cout.flush();
			}

			in_file.read((char *) &in_cache, cache_pix*IN_BYT_PER_PIX);
			index_file.read((char *) &index_cache, cache_pix*INDEX_BYT_PER_PIX);
			vcf_file.read((char *) &vcf_cache, cache_pix);
			lcv_file.read((char *) &lcv_cache, cache_pix);

			tmp_global_count = sum_vcf_lcv(in_cache, index_cache, vcf_cache, lcv_cache, cache_pix, temp_total, temp_count, temp_total2, temp_count2, temp_total3, temp_count3, temp_total4, temp_count4);
			
			//Dump temp values into final total variables
			for(int j=0; j<NBINS; j++) {
				if(temp_count[j]) {
					total[j] += temp_total[j];
					count[j] += temp_count[j];
				}
				if(temp_count2[j]) {
					total2[j] += temp_total2[j];
					count2[j] += temp_count2[j];
				}
				if(temp_count3[j]) {
					total3[j] += temp_total3[j];
					count3[j] += temp_count3[j];
				}
				if(temp_count4[j]) {
					total4[j] += temp_total4[j];
					count4[j] += temp_count4[j];
				}
			}
			global_total += tmp_global_count.global_total;
			global_count += tmp_global_count.global_count;
			global_total2 += tmp_global_count.global_total2;
			global_count2 += tmp_global_count.global_count2;
			global_total3 += tmp_global_count.global_total3;
			global_count3 += tmp_global_count.global_count3;
			global_total4 += tmp_global_count.global_total4;
			global_count4 += tmp_global_count.global_count4;
		}

		std::cout << "Processing remainder...\n";
		std::cout.flush();
		//check for remainder
		if(remainder_pix) {
			in_t *remain_in_cache = new in_t[remainder_pix];
			index_t *remain_index_cache = new index_t[remainder_pix];
			unsigned char *remain_vcf_cache = new unsigned char[remainder_pix];
			unsigned char *remain_lcv_cache = new unsigned char[remainder_pix];
			in_file.read((char *) remain_in_cache, remainder_pix*IN_BYT_PER_PIX);
			index_file.read((char *) remain_index_cache, remainder_pix*INDEX_BYT_PER_PIX);
			vcf_file.read((char *) &remain_vcf_cache, remainder_pix);
			lcv_file.read((char *) &remain_lcv_cache, remainder_pix);

			tmp_global_count = sum_vcf_lcv(remain_in_cache, remain_index_cache, remain_vcf_cache, remain_lcv_cache, remainder_pix, temp_total, temp_count, temp_total2, temp_count2, temp_total3, temp_count3, temp_total4, temp_count4);
	
			//Dump temp values into final total variables
			for(int j=0; j<NBINS; j++) {
				if(temp_count[j]) {
					total[j] += temp_total[j];
					count[j] += temp_count[j];
				}
				if(temp_count2[j]) {
					total2[j] += temp_total2[j];
					count2[j] += temp_count2[j];
				}
				if(temp_count3[j]) {
					total3[j] += temp_total3[j];
					count3[j] += temp_count3[j];
				}
				if(temp_count4[j]) {
					total4[j] += temp_total4[j];
					count4[j] += temp_count4[j];
				}
			}
			global_total += tmp_global_count.global_total;
			global_count += tmp_global_count.global_count;
			global_total2 += tmp_global_count.global_total2;
			global_count2 += tmp_global_count.global_count2;
			global_total3 += tmp_global_count.global_total3;
			global_count3 += tmp_global_count.global_count3;
			global_total4 += tmp_global_count.global_total4;
			global_count4 += tmp_global_count.global_count4;
			
			delete[] remain_in_cache;
			delete[] remain_index_cache;
			delete[] remain_vcf_cache;
			delete[] remain_lcv_cache;
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
	out_file.open(argv[5]);
	out_file << "Country Code, NPixels, Total_AGB_Tg, NPixels_vcf" << VCF_THRESH2 << ", Total_AGB_Tg_vcf" << VCF_THRESH2;
	out_file << ", NPixels_vcf" << VCF_THRESH3 << ", Total_AGB_Tg_vcf" << VCF_THRESH3;
	out_file << ", NPixels_vcf" << VCF_THRESH4 << ", Total_AGB_Tg_vcf" << VCF_THRESH4 << '\n';
	double pixel_area = PIXEL_SIZE * PIXEL_SIZE;  //pixel area in m^2
	for (int i=0; i<NBINS; i++) {
		if (count[i] > 0) {
			#ifdef INPUT_INT
				//integer input time is scaled by 10
				out_file << i << ',' << count[i] << ',' << total[i] / 10 * (pixel_area / 10000) / 1000000;
				out_file << ',' << count2[i] << ',' << total2[i] / 10 * (pixel_area / 10000) / 1000000;
				out_file << ',' << count3[i] << ',' << total3[i] / 10 * (pixel_area / 10000) / 1000000;
				out_file << ',' << count4[i] << ',' << total4[i] / 10 * (pixel_area / 10000) / 1000000 << '\n';
			#else
				out_file << i << ',' << count[i] << ',' << total[i] * (pixel_area / 10000) / 1000000;
				out_file << ',' << count2[i] << ',' << total2[i] * (pixel_area / 10000) / 1000000;
				out_file << ',' << count3[i] << ',' << total3[i] * (pixel_area / 10000) / 1000000;
				out_file << ',' << count4[i] << ',' << total4[i] * (pixel_area / 10000) / 1000000 << '\n';
			#endif
		}
	}

	#ifdef INPUT_INT
		//integer input time is scaled by 10
		out_file << NBINS << ',' << global_count << ',' << global_total/10 * (pixel_area / 10000) / 1000000;
		out_file << ',' << global_count2 << ',' << global_total2/10 * (pixel_area / 10000) / 1000000;
		out_file << ',' << global_count3 << ',' << global_total3/10 * (pixel_area / 10000) / 1000000;
		out_file << ',' << global_count4 << ',' << global_total4/10 * (pixel_area / 10000) / 1000000 << '\n';
	#else
		out_file << NBINS << ',' << global_count << ',' << global_total * (pixel_area / 10000) / 1000000;
		out_file << ',' << global_count2 << ',' << global_total2 * (pixel_area / 10000) / 1000000;
		out_file << ',' << global_count3 << ',' << global_total3 * (pixel_area / 10000) / 1000000;
		out_file << ',' << global_count4 << ',' << global_total4 * (pixel_area / 10000) / 1000000 << '\n';
	#endif

	out_file.close();

	return 0;
}

