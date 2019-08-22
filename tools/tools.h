#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <random>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <math.h>       /* sqrt */

#ifndef TOOLS_H
#define TOOLS_H

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

/// Velocity on the parabolic Poiseuille profile
template<typename T1, typename T2>
T2 poiseuilleVelocity(T1 r, T2 uMax, T1 powerPoiseuilleVel, T1 inletRadius);


template<typename T>
T getGradient(T firstValue, T secondValue, T step);

template<typename T>
bool isInCircle(T x,T y, T rSq, T cX,T cY);

template<typename T>
bool isInTube(T x,T y,T z, T rSq, T cX,T cY,T Z0,T Z1);

template<typename T>
bool isInSphere(T x,T y,T z, T rSq, T cX,T cY,T cZ);

template<typename T1, typename T2>
bool isInParabola(T1 x,T1 y, T2 steepness, T1 cX,T1 cY);

template<typename T1, typename T2>
bool isInParaboloid(T1 x,T1 y,T1 z, T2 steepness, T1 cX,T1 cY,T1 cZ);

template<typename T1>
double dotProduct(T1 aX,T1 aY,T1 aZ,
                  T1 bX,T1 bY,T1 bZ);

template<typename T1>
std::vector<double> normalizeVector(T1 x,T1 y,T1 z);

template<typename T1>
bool isBelowPlane(T1  x,T1  y,T1  z,
                  T1 nX,T1 nY,T1 nZ, double d);

template<typename T>
bool heaviside(T x, T cX, bool orientedPositive=true);

template<typename T>
int sign(T value);

template<typename T>
double distanceFromPoint(T x,T y, T PX,T PY);

bool correctArguments(const int argc, char const* const* argv);

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

//https://stackoverflow.com/questions/26312570/calculate-surface-area-of-a-3d-mesh
template<typename T>
T areaOfTriangle(T p1X, T p1Y, T p1Z,
                 T p2X, T p2Y, T p2Z,
                 T p3X, T p3Y, T p3Z);

template<typename T>
T meshSurface(T *X, T *Y, T *Z,
              int numT, int *V1, int *V2, int *V3);


#endif
