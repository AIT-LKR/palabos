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

#ifndef TOOLS_GEOMETRIC_H
#define TOOLS_GEOMETRIC_H

using namespace std;

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

//https://stackoverflow.com/questions/26312570/calculate-surface-area-of-a-3d-mesh
template<typename T>
T areaOfTriangle(T p1X, T p1Y, T p1Z,
                 T p2X, T p2Y, T p2Z,
                 T p3X, T p3Y, T p3Z);

template<typename T>
T meshSurface(T *X, T *Y, T *Z,
              int numT, int *V1, int *V2, int *V3);

template<typename T1, typename T2>
std::vector<T2> calculateCOM(std::vector<std::vector<T1>> coordinates);

template<typename T1, typename T2>
T2 diameterFromSphereVolume(T1 volume);


#endif
