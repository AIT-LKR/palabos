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

/* **************************************************************************** */

template<typename T1, typename T2>
PoiseuilleVelocity<T1,T2>::PoiseuilleVelocity(T2 uLB_, T2 uDev_, T1 powerPoiseuilleVel_, char dir_,
        T1 inletCentre_, T1 inletRadius_)
    : uLB(uLB_), uDev(uDev_), powerPoiseuilleVel(powerPoiseuilleVel_), dir(dir_),
     inletCentreA(inletCentre_), inletCentreB(-1), inletRadius(inletRadius_),
     seed(std::chrono::system_clock::now().time_since_epoch().count())
{
    // random number generator
    generator = new std::default_random_engine(seed);
    main_distribution = new std::normal_distribution<T2>(uLB-uDev, uDev);
    side_distribution = new std::normal_distribution<T2>(0., uDev);
}

template<typename T1, typename T2>
PoiseuilleVelocity<T1,T2>::PoiseuilleVelocity(T2 uLB_, T2 uDev_, T1 powerPoiseuilleVel_, char dir_,
        T1 inletCentreA_, T1 inletCentreB_, T1 inletRadius_)
    : uLB(uLB_), uDev(uDev_), powerPoiseuilleVel(powerPoiseuilleVel_), dir(dir_),
      inletCentreA(inletCentreA_), inletCentreB(inletCentreB_), inletRadius(inletRadius_),
      seed(std::chrono::system_clock::now().time_since_epoch().count())
{
    // random number generator
    generator = new std::default_random_engine(seed);
    main_distribution = new std::normal_distribution<T2>(uLB-uDev, uDev);
    side_distribution = new std::normal_distribution<T2>(0., uDev);
}

template<typename T1, typename T2>
void PoiseuilleVelocity<T1,T2>::operator()(T1 iX, T1 iY, T1 iZ, Array<T2,3>& u) const {
    T2 main_velocity = uLB;
    T2 side_velocity = 0.;
    T2 down_velocity = 0.;
    if (uDev) {
        main_velocity = (*main_distribution)(*generator);
        side_velocity = (*side_distribution)(*generator)/4;
        down_velocity = (*side_distribution)(*generator)/4;
    }
    T1 r;

    if (dir == 'x') {
        if (inletCentreB == -1)
            r = distanceFromPoint<T2>(iY,iZ, inletCentreA,iZ);
        else
            r = distanceFromPoint<T2>(iY,iZ, inletCentreA,inletCentreB);
        u[0] = poiseuilleVelocity(r, main_velocity, powerPoiseuilleVel, inletRadius);
        u[1] = side_velocity;
        u[2] = down_velocity;
    }
    if (dir == 'y') {
        if (inletCentreB == -1)
            r = distanceFromPoint<T2>(iX,iZ, inletCentreA,iZ);
        else
            r = distanceFromPoint<T2>(iX,iZ, inletCentreA,inletCentreB);
        u[0] = side_velocity;
        u[1] = poiseuilleVelocity(r, main_velocity, powerPoiseuilleVel, inletRadius);
        u[2] = down_velocity;
    }
}

template<typename T1, typename T2>
void PoiseuilleVelocity<T1,T2>::operator()(T1 iX, T1 iY, Array<T2,2>& u) const {
    T2 main_velocity = uLB;
    T2 side_velocity = 0.;
    if (uDev) {
        main_velocity = (*main_distribution)(*generator);
        side_velocity = (*side_distribution)(*generator)/2;
    }

    if( dir == 'x' ) {
        u[0] = poiseuilleVelocity(iY-inletCentreA, main_velocity, powerPoiseuilleVel, inletRadius);
        u[1] = side_velocity;
    }
    if( dir == 'y' ) {
        u[0] = side_velocity;
        u[1] = poiseuilleVelocity(iX-inletCentreA, main_velocity, powerPoiseuilleVel, inletRadius);
    }
}

#endif
