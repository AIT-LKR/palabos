#include <../src/core/globalDefs.h>
#include <../src/io/parallelIO.h>
#include <../src/parallelism/mpiManager.h>

#include <string>
#include <sstream>
#include <math.h>       /* sqrt */

typedef double T;

#ifndef PHYSICAL_FLOW_PARAM_H
#define PHYSICAL_FLOW_PARAM_H

using namespace std;
using namespace plb;

/// Numeric parameters for isothermal, incompressible flow.
template<typename T>
class PhysicalFlowParam {
public:
    /// Constructor
    PhysicalFlowParam(T lbU0_, T Re_, T phNu_, T phL0_,
            plint nodesPerMeter_, T lx_, T ly_, T lz_)
        : lbU0(lbU0_),
          Re(Re_),
          phNu(phNu_),
          phL0(phL0_),
          nodesPerMeter(nodesPerMeter_),
          lx(lx_), ly(ly_), lz(lz_)
    {}

    T getLatticeU() const             { return lbU0; }
    T getRe() const                   { return Re; }
    T getPhysicalNu() const           { return phNu; }
    T getPhysicalLength(T l=1.) const { return l * phL0; } // l=1. -> char. length (physical units)
    plint getNodesPerMeter() const    { return nodesPerMeter; }
    // dimensionless extend of system:
    T getlx() const                { return lx; }
    T getly() const                { return ly; }
    T getlz() const                { return lz; }

    plint getResolution() const    { return nodesPerMeter * getPhysicalLength(); }
    T getPhysicalU() const         { return Re * getPhysicalNu() / getPhysicalLength(); }
    T getPhysicalT(T t=1.) const   { return t * getPhysicalLength() / getPhysicalU(); } // l=1 -> char time (physical units)
    T getDeltaX() const            { return 1 / (T)getResolution(); }
    T getDeltaT() const            { return getDeltaX() * getLatticeU(); }
    T getLatticeNu() const         { return getLatticeU()*(T)getResolution()/Re; }
    //T getLatticeNu() const         { return getDeltaT() / getDeltaX() / getDeltaX() / Re; }
    T getTau() const               { return (T)3 * getLatticeNu() + (T)0.5; }
    T getOmega() const             { return 1 / getTau(); }
    T getMa() const                { return getLatticeU() / sqrt(1/3); }
    // physical extend of system:
    T getLx() const                { return lx * getPhysicalLength(); }
    T getLy() const                { return ly * getPhysicalLength(); }
    T getLz() const                { return lz * getPhysicalLength(); }
    plint getNx() const                { return getLx() * (T)nodesPerMeter + 1; }
    plint getNy() const                { return getLy() * (T)nodesPerMeter + 1; }
    plint getNz() const                { return getLz() * (T)nodesPerMeter + 1; }

    /// conversion from physical to dimensionless units for time coordinate
    T dimlessTime (T phTime) const   { return phTime / getPhysicalT(); }
    /// conversion from lattice units to dimensionless for time coordinate
    T stepToDimless(plint iT) const  { return iT*(T)getDeltaT(); }
    /// conversion from physical to lattice units for time coordinate
    plint phStep(T phTime) const   { return nStep(dimlessTime(phTime));        }
    /// conversion from dimensionless to lattice units for time coordinate
    plint nStep(T t) const         { return (int)(t / (T)getDeltaT() +(T)0.5); }
    /// calculation of mass flow
    T getMassFlow(T rho0, T area) const { return rho0 * area * getPhysicalU(); }
    /// physical pressure
    T getPhPressure0(T lbRho) const { return lbRho * getDeltaX()*getDeltaX()/getDeltaT()/getDeltaT() *(1/3); }
    T getPhPressure1(T lbRho) const { return lbRho * getPhysicalU()*getPhysicalU() /getMa()/getMa(); }

private:
    T lbU0;               // lattice velocity
    T Re;                 // Reynolds number =l0^2/(t0*nu) = u0*l0/nu
    T phNu;               // viscosity               [kg/m/s]
    T phL0;               // characteristic lenght   [m]
    plint nodesPerMeter;  // cells per meter
    T lx, ly, lz;         // extent of the system (dimensionless units)
};

/*
template<typename T>
T getPressure1(T rho_ph, T dx, T dt, T rho) {
return NSDESCRIPTOR<T>::cs2*(rho-1.)*dx*dx/(dt*dt)*rho_ph;
}
template<typename T>
T getPressure2(T rho_ph, T dx, T dt, T rho) {
return rho_ph * (dx*dx/dt/dt) * rho;
}
*/


#endif
