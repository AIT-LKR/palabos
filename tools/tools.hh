#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <random>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <math.h>       /* sqrt */

#include "../tools/tools.h"

#ifndef TOOLS_HH
#define TOOLS_HH

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

/// Velocity on the parabolic Poiseuille profile
template<typename T1, typename T2>
T2 poiseuilleVelocity(T1 r, T2 uMax, T1 powerPoiseuilleVel, T1 inletRadius) {
    // r_dash represents shift from zero to new inletCentre and scaling to
    // a parabola across inlet
    T2 r_dash = T2(r)/inletRadius;
    if (r_dash > 1.) return 0;
    T2 r_dash_sqr = pow( r_dash, powerPoiseuilleVel);
    T2 vel = -1*uMax*( r_dash_sqr - 1);
    return vel;
}

template<typename T>
T getGradient(T firstValue, T secondValue, T step) {
    return (secondValue-firstValue)/step;
}

template<typename T>
bool isInCircle(T x,T y, T rSq, T cX,T cY) {
    if (rSq <= 0) return false;
    T xr = x-cX; // relative to centre
    T yr = y-cY;
    if (xr*xr + yr*yr > rSq) return false;
    return true;
}

template<typename T>
bool isInTube(T x,T y,T z, T rSq, T cX,T cY,T Z0,T Z1) {
    if (Z0 > Z1) return false;
    if (z < Z0)  return false;
    if (z > Z1)  return false;
    return isInCircle(x,y,rSq,cX,cY);
}

template<typename T>
bool isInSphere(T x,T y,T z, T rSq, T cX,T cY,T cZ) {
    if (rSq <= 0) return false;
    T xr = x-cX; // relative to centre
    T yr = y-cY;
    T zr = z-cZ;
    if (xr*xr + yr*yr + zr*zr > rSq) return false;
    return true;
}

template<typename T1, typename T2>
bool isInParabola(T1 x,T1 y, T2 steepness, T1 cX,T1 cY) {
    if (steepness == 0) {
        return false;
    }
    // check if lies below parabola:
    if (sign(steepness) == 1) {
        if (y < cY) return false;
    } else {
        if (y > cY) return false;
    }
    // check if too far away:
    double horizontalDistance = abs(x - cX);
    double radiusAtHeight = sqrt( abs(double(y - cY)/double(steepness)) );
    if (horizontalDistance > radiusAtHeight) return false;

    return true;
}

template<typename T1, typename T2>
bool isInParaboloid(T1 x,T1 y,T1 z, T2 steepness, T1 cX,T1 cY,T1 cZ) {
    if (steepness == 0) {
        return false;
    }
    // check if lies below paraboloid:
    if (sign(steepness) == 1) {
        if (z < cZ) return false;
    } else {
        if (z > cZ) return false;
    }
    // check if too far away:
    double horizontalDistance = distanceFromPoint(x,y, cX,cY);
    double radiusAtHeight = sqrt( abs(double(z - cZ)/double(steepness)) );
    if (horizontalDistance > radiusAtHeight) return false;

    return true;
}


template<typename T1>
double dotProduct(T1 aX,T1 aY,T1 aZ,
                  T1 bX,T1 bY,T1 bZ) {
    double result = aX*bX + aY*bY + aZ*bZ;
    return result;
}

template<typename T1>
std::vector<double> normalizeVector(T1 x,T1 y,T1 z) {
    double length =  sqrt(dotProduct(x,y,z, x,y,z));
    double nX = x / length;
    double nY = y / length;
    double nZ = z / length;
    std::vector<double> normedVec = {nX,nY,nZ};
    return normedVec;
}

template<typename T1>
bool isBelowPlane(T1  x,T1  y,T1  z,
                  T1 nX,T1 nY,T1 nZ, double d) {
    double distance = dotProduct(x,y,z, nX,nY,nZ);
    if (distance < d)
        return true;
    else
        return false;
}

template<typename T>
bool heaviside(T x, T cX, bool orientedPositive) {
    if (x == cX) return true;
    if (x > cX) {
        if (orientedPositive) return true;
        else                  return false;
    } else {
        if (orientedPositive) return false;
        else                  return true;
    }
}

template<typename T>
int sign(T value) {
    return (T(0) < value) - (value < T(0));
}


template<typename T>
double distanceFromPoint(T x, T y, T PX,T PY) {
    T xr = x-PX; // relative to centre
    T yr = y-PY;
    return sqrt(xr*xr + yr*yr);
}

template<typename T>
double distanceFromPoint(T x,T y,T z, T PX,T PY,T PZ) {
    T xr = x-PX; // relative to centre
    T yr = y-PY;
    T zr = z-PZ;
    return sqrt(xr*xr + yr*yr, zr*zr);
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
T areaOfTriangle(T p1X, T p1Y, T p1Z,
                 T p2X, T p2Y, T p2Z,
                 T p3X, T p3Y, T p3Z)
{
    T ax = p2X - p1X;
    T ay = p2Y - p1Y;
    T az = p2Z - p1Z;
    T bx = p3X - p1X;
    T by = p3Y - p1Y;
    T bz = p3Z - p1Z;
    T cx = ay*bz - az*by;
    T cy = az*bx - ax*bz;
    T cz = ax*by - ay*bx;

    return 0.5 * sqrt(cx*cx + cy*cy + cz*cz);
}

template<typename T>
T meshSurface(T *X, T *Y, T *Z,
              int numT, int *V1, int *V2, int *V3)
{
    T area = 0.;
    for (int n=0; n<numT; n++) {
        area += areaOfTriangle(X[V1[n]], Y[V1[n]], Z[V1[n]],
                               X[V2[n]], Y[V2[n]], Z[V2[n]],
                               X[V3[n]], Y[V3[n]], Z[V3[n]]);
    }
    return area;
}


#endif
