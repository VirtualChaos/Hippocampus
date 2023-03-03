#ifndef _spikeBinner_h
#define _spikeBinner_h

void spikeBinner(const double* flat_spiketimes, const double* stc, const double* conditions, int repeat, double* consol_arr, double* binArr, int flat_spiketimes_size, int stc_size, int total_nBins);

#endif
