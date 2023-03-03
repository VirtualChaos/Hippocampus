#include <stdio.h>
#include <stdlib.h>
#include "spikeBinner.h"

void spikeBinner(const double* flat_spiketimes, const double* stc, const double* conditions, int repeat, double* consol_arr, double* binArr, int flat_spiketimes_size, int stc_size, int total_nBins)
{
    int interval = 0;
    int bins_hit = 0;
    int tbins_hit = 0;
    int sp = 0;
    int temp_col;
    
    for (sp=0; sp<flat_spiketimes_size; sp++) {
        while (interval < stc_size) {
            if ((flat_spiketimes[sp] >= stc[interval]) && (flat_spiketimes[sp] < stc[interval+1])) {
                break;
            }
            interval++;
        }
        bins_hit = stc[stc_size*(3-1)+(interval-1)]; //[interval][2]
        if (conditions[interval] && (bins_hit >= 1)) {
            temp_col = flat_spiketimes[flat_spiketimes_size+sp-1];
            consol_arr[100*(temp_col-1)+(bins_hit-1)]++; //[bins_hit][flat_spiketimes[sp][1]]
        }
        
        if (repeat == 1) {
            tbins_hit = stc[stc_size*(2-1)+(interval-1)]; //[interval][1]
            if (conditions[interval] && (tbins_hit >= 1)) {
            temp_col = flat_spiketimes[flat_spiketimes_size+sp-1];
            binArr[total_nBins*(temp_col-1)+(tbins_hit-1)]++; //[tbins_hit][flat_spiketimes[sp][1]]
            }
        }
    }
    
    // int row = 4;
    // int col = 3;
    // binArr[total_nBins*(col-1)+(row-1)] = 5; // row, col as in matlab
}
