#ifndef TOOLS_HH
#define TOOLS_HH

#include <string>
#include <vector>
#include <unordered_map>
#include <iterator>

#include "../tools/tools.h"

using namespace std;


template<typename T1, typename T2>
void changeElementsInMap(std::unordered_map<T1,T2> &dataMap, std::vector<T1> sorting, std::vector<T2> newData) {

    for (unsigned int i=0; i<sorting.size(); i++) {
        dataMap[sorting[i]] = newData[i];
	}
}

template<typename T>
std::string createStringFromVector(std::vector<T> data, std::string separator) {
    std::stringstream ss;
    for (unsigned int i=0; i<data.size(); i++) {
        ss << data[i];
        if (i == data.size()-1) break;
        ss << separator;
	}
    ss << endl;
    return ss.str();
}

template<typename T1, typename T2>
std::string createStringFromMap(std::unordered_map<T1,T2> dataMap, std::vector<T1> sorting, std::string separator) {
    std::stringstream ss;
    for (unsigned int i=0; i<sorting.size(); i++) {
        ss << dataMap[sorting[i]];
        if (i == sorting.size()-1) break;
        ss << separator;
	}
    ss << endl;
    return ss.str();
}

std::vector<string> splitString(std::string input) {
    std::istringstream iss(input);
    std::vector<string> result{istream_iterator<string>{iss},
                          istream_iterator<string>{}};
    return result;
}

template<typename T, typename ofstream>
void writeVec2File(std::vector<T> vec, std::string fileName, ofstream ofile) {
    ofile.open(fileName, std::ios_base::out | std::ios_base::trunc );
    for (const auto &e : vec) ofile << e << "\n";
    return;
}

std::vector<double> convertStringToDouble(std::vector<string> stringVector) {
    vector<double> doubleVector(stringVector.size());
    transform(stringVector.begin(), stringVector.end(), doubleVector.begin(),
              [](string const& val) {return stod(val);});
    return doubleVector;
}

std::string readLastLineInFile(std::string fileName) {
    std::string result;
    std::ifstream read(fileName.c_str(), std::ios_base::ate );
    if( !read ) return "";

    int length = 0;
    char c = '\0';

    length = read.tellg();//Get file size

    // loop backward over the file
    for(int i = length-2; i > 0; i-- ) {
        read.seekg(i);
        c = read.get();
        if( c == '\r' || c == '\n' )//new line?
             break;
    }
    std::getline(read, result);//read last line
    return result;
}

string messageLost(string fileName) {
    std::stringstream message;
    message << "Cannot find:" << endl << " '" << fileName << "'" << endl;
    return message.str();
}

string messageSaving(std::string fileName) {
    std::stringstream message;
    message << "Saving:" << endl << " '" << fileName << "'" << endl;
    return message.str();
}

string messageLoading(std::string fileName) {
    std::stringstream message;
    message << "Loading: '" << fileName << "'" << endl;
    return message.str();
}

template<typename T1>
T1 randNumber<T1>::operator()(T1 randMin, T1 randMax) {
    float factor = ( randoom() /(m-1) );
    T1 value = randMin + static_cast<float>(randMax-randMin)*factor;
    return value; }

template<typename T1>
unsigned int randNumber<T1>::randoom() {
    return lastEntry = (a*lastEntry) %m; }

/// Time ramp for velocity
template<typename T1, typename T2>
T2 rampVelocity(T2 vel, T1 currentTime, T1 timeMax) {
    if (currentTime > timeMax) return vel;
    else return vel*currentTime/timeMax;
}

template<typename T>
T getGradient(T firstValue, T secondValue, T step) {
    return (secondValue-firstValue)/step;
}

belowThresholdCounter::belowThresholdCounter(double threshold_, int characteristicCount_)
    : threshold(threshold_),
      characteristicCount(characteristicCount_),
      counter(0)
{ }

bool belowThresholdCounter::repeatedlyBelowThrehold(double currentValue) {
    if (abs(currentValue) <= threshold) counter++;
    else counter=0;

    if (counter >= characteristicCount) return true;
    else return false;
}

void belowThresholdCounter::reset() {
    counter=0;
}

std::string ReadNthLine(const std::string& filename, int N)
{
   std::ifstream in(filename.c_str());

   std::string s;
   //for performance
   s.reserve(100);

   //skip N lines
   for(int i = 0; i < N; ++i)
       std::getline(in, s);

   std::getline(in,s);
   return s;
}

double numberInString(std::string text)
{
    std::string tmp1, tmp2;
    double value;
    std::stringstream ss(text);
    ss >> tmp1 >> tmp2 >> value;
    return value;
}

template<typename T>
inRange<T>::inRange(T min_, T max_)
    : min(min_), max(max_)
{}
template<typename T>
bool inRange<T>::operator()(T value) const {
    if (value < min || value > max)
        return false;
    else
        return true;
}


#endif
