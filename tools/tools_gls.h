#ifndef TOOLS_GLS_H
#define TOOLS_GLS_H

#include <vector>
#include <gsl/gsl_histogram.h>


template<typename T1, typename T2>
void histFromVect(gsl_histogram* hist, T1 min, T1 max, std::vector<T2> data);
template<typename T>
std::vector<T> vectFromHist(gsl_histogram * hist);


#endif
