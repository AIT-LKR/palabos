#ifndef TOOLS_GEOMETRIC_HH
#define TOOLS_GEOMETRIC_HH

#include <vector>
#include <math.h>       /* sqrt */

#include "./tools_geometric.h"

using namespace std;


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

template<typename T1, typename T2>
std::vector<T2> calculateCOM(std::vector<std::vector<T1>> points) {
    std::vector<T2> center;
    T2 X = 0;
    T2 Y = 0;
    T2 Z = 0;
    int numberOfVoxels = points.size();
    for (int i=0; i<numberOfVoxels; i++) {
        std::vector<T1> currentVoxel = points.at(i);
        X += currentVoxel.at(0);
        Y += currentVoxel.at(1);
        Z += currentVoxel.at(2);
    }
    X /= numberOfVoxels;
    Y /= numberOfVoxels;
    Z /= numberOfVoxels;

    center.push_back( X );
    center.push_back( Y );
    center.push_back( Z );
    return center;
}

template<typename T1, typename T2>
T2 diameterFromSphereVolume(T1 volume) {
    T2 diameter = 0;
    diameter = 2*cbrt( 3.*(T2)volume /4./M_PI );
    return diameter;
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

std::vector<float> rotate(float x, float y, float angle) {
    float x_dash = cos(angle)*x + sin(angle)*y;
    float y_dash = -1*sin(angle)*x + cos(angle)*y;
    std::vector<float> rotatedCoords = {x_dash,y_dash};
    return rotatedCoords;
}


#endif
