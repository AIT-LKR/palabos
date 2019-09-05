#ifndef TOOLS_GLS_HH
#define TOOLS_GLS_HH

#include <vector>
#include <gsl/gsl_histogram.h>

#include "./tools_gls.h"


// https://www.gnu.org/software/gsl/doc/html/histogram.html
template<typename T1, typename T2>
void histFromVect(gsl_histogram* hist, T1 min, T1 max, std::vector<T2> data) {
    gsl_histogram_set_ranges_uniform(hist, min, max);
    for(unsigned int i = 0; i < data.size(); i++ )
        if( gsl_histogram_increment( hist, data.at(i) ) != 0 )
            exit(1);
}


template<typename T>
std::vector<T> vectFromHist(gsl_histogram * hist) {
    unsigned int length = gsl_histogram_bins(hist);
    std::vector<T> vect(length);
    for( unsigned int i=0; i<length; i++ ) vect[i] = gsl_histogram_get(hist, i);

    return vect;
}


#endif
