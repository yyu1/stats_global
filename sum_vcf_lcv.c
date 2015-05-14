//Given pointers to caches that have been read from disk already. Calculate the sums


inline double get_carbon_frac(unsigned char lcv) {
  switch (lcv) {
    case 40:
      return (double) 0.48;
		case 160:
			return (double) 0.48;
    case 50:
      return (double) 0.48;
    case 60:
      return (double) 0.48;
    case 70:
      return (double) 0.51;
    case 90:
      return (double) 0.51;
		case 100:
			return (double) 0.495;
    default:
      return (double) 0.47;
  }
}

struct global_count_block {
	double global_total;
	long long global_count;
	double global_total2;
	long long global_count2;
	double global_total3;
	long long global_count3;
};


global_count_block sum_vcf_lcv(in_t *in_cache, index_t *index_cache, unsigned char *vcf_cache, unsigned char *lcv_cache, unsigned long long npix, double *tmp_total, long long *tmp_count, double *tmp_total2, long long *tmp_count2, double *tmp_total3, long long *tmp_count3) {

	global_count_block my_count;
	my_count.global_total = 0;
	my_count.global_count = 0;
	my_count.global_total2 = 0;
	my_count.global_count2 = 0;
	my_count.global_total3 = 0;
	my_count.global_count3 = 0;

	//Reset all the temp counters to 0 since we are just starting the local sum here
	for (int i=0; i<NBINS; i++) {
		tmp_total[i] = 0;
		tmp_count[i] = 0;
		tmp_total2[i] = 0;
		tmp_count2[i] = 0;
		tmp_total3[i] = 0;
		tmp_count3[i] = 0;
	}

	register double carbon_frac;

	for (unsigned long long i=0; i<npix; i++) {

		if (in_cache[i] > 0) {
			
			carbon_frac = get_carbon_frac(lcv_cache[i]);
			
			my_count.global_total += in_cache[i] * carbon_frac;
			++my_count.global_count;

			if ((index_cache[i] >= 0) && (index_cache[i] < NBINS)) {    // ignore index that fall outside (i.e. -1 for water)
				tmp_total[index_cache[i]] += in_cache[i] * carbon_frac;
				++tmp_count[index_cache[i]];
			}

			//VCF Thresh2
			if (vcf_cache[i] >= VCF_THRESH2) {
				my_count.global_total2 += in_cache[i] * carbon_frac;
				++my_count.global_count2;

				if ((index_cache[i] >= 0) && (index_cache[i] < NBINS)) {
					tmp_total2[index_cache[i]] += in_cache[i] * carbon_frac;
					++tmp_count2[index_cache[i]];
				}

				//VCF Thresh3

				if (vcf_cache[i] >= VCF_THRESH3) {
					my_count.global_total3 += in_cache[i] * carbon_frac;
					++my_count.global_count3;

					if ((index_cache[i] >= 0) && (index_cache[i] < NBINS)) {
						tmp_total3[index_cache[i]] += in_cache[i] * carbon_frac;
						++tmp_count3[index_cache[i]];
					}
				}
			}
		}
	}
	return my_count;
}
