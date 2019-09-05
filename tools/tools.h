#ifndef TOOLS_H
#define TOOLS_H

#include <string>
#include <vector>
#include <unordered_map>


using namespace std;


template<typename T1, typename T2>
void changeElementsInMap(std::unordered_map<T1,T2> &dataMap, std::vector<T1> sorting, std::vector<T2> newData);

template<typename T>
std::string createStringFromVector(std::vector<T> data, std::string separator);

template<typename T1, typename T2>
std::string createStringFromMap(std::unordered_map<T1,T2> dataMap, std::vector<T1> sorting, std::string separator);

std::vector<string> splitString(std::string input);

template<typename T>
void writeVec2File(std::vector<T> vec, std::string name);

std::vector<double> convertStringToDouble(std::vector<string> stringVector);

std::string readLastLineInFile(std::string fileName);

string messageLost(string fileName);

string messageSaving(std::string fileName);

string messageLoading(std::string fileName);

template<typename T1>
class randNumber {
    public:
        randNumber(unsigned int seed_) : lastEntry(seed_) { }
        T1 operator()(T1 randMin, T1 randMax);
        unsigned int randoom();
    private:
        unsigned int a = 16807;
        unsigned int m = 2147483647;
        unsigned int lastEntry;
};

/// Time ramp for velocity
template<typename T1, typename T2>
T2 rampVelocity(T2 vel, T1 currentTime, T1 timeMax);

template<typename T>
T getGradient(T firstValue, T secondValue, T step);

class belowThresholdCounter {
    public:
        belowThresholdCounter(double threshold_, int characteristicCount_);
        bool repeatedlyBelowThrehold(double currentValue);
        void reset();
    private:
        double threshold;
        int characteristicCount;
        int counter;
};

template<typename T1, typename T2>
void addToAverage(T1& average, T1 newValue, T2& n) {
    T1 nMinusOne = (T1) n - (T1) 1;
    T1 oneOverN = (T1) 1 / (T1) n;

    average = oneOverN * (nMinusOne * average + newValue);
}

template<typename T>
class inRange{
public:
    inRange(T min_, T max_);
    bool operator()(T value) const;
private:
    T min, max;
};


#endif
