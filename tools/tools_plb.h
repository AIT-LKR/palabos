#include <../src/core/globalDefs.h>
#include <../src/io/parallelIO.h>
#include <../externalLibraries/tinyxml/tinyxml.h>

#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <random>

#include "./tools.h"
#include "./tools.hh"
#include "./plb_physicalFlowParam.h"

typedef double T;

#ifndef TOOLS_PLB_H
#define TOOLS_PLB_H

using namespace std;
using namespace plb;


void renameAllExisting(string outDir, string fileName);

bool copyFile(string sourceName, string destinationName);

void readFile(std::vector<std::vector<T>> vecs, std::string name);

void rewriteOutputFile(std::string data, std::string fileName);

void appendOutputFile(std::string data, std::string fileName);


template<typename T>
class PoiseuilleVelocity {
    public:
        PoiseuilleVelocity(T uLB_, T uDev_, plint powerPoiseuilleVel_, char dir_,
                plint inletCentre_, plint inletRadius_);
        PoiseuilleVelocity(T uLB_, T uDev_, plint powerPoiseuilleVel_, char dir_,
                plint inletCentreA_, plint inletCentreB_, plint inletRadius_);
        void operator()(plint iX, plint iY, Array<T,2>& u) const;
        void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const;
    private:
        T uLB;
        T uDev;
        plint powerPoiseuilleVel;
        char dir;
        plint inletCentreA;
        plint inletCentreB;
        plint inletRadius;
        unsigned seed;
        std::default_random_engine* generator;
        std::normal_distribution<>* main_distribution;
        std::normal_distribution<>* side_distribution;
};


Box3D layer(plint Nx, plint Ny, plint i) {
    // create slice with width=Nx, length=Ny, height=1
    return Box3D(0,Nx-1, 0,Ny-1, i,i);
}

// inspired by src/io/imageWriter.h
namespace plb {
template<typename T>
class TopoWriter {
public:
    TopoWriter(plint valueRange_);
    void writePpm(
        std::string const& fName, bool colored,
        MultiScalarField2D<T>& field) const;
    void writePpm (
        std::string const& fName, bool colored,
        MultiScalarField3D<T>& field, int dir=2) const;
    void writeBinary(
        std::string const& fName,
        MultiScalarField3D<T>& field) const;

private:
    void writePpmImplementation (
        std::string const& fName, bool colored,
        ScalarField2D<T>& localField) const;
    void writeBinaryImplementation (
        std::string const& fName,
        ScalarField3D<T>& localField) const;
private:
    plint valueRange;
};
}

template<typename T>
class inRange{
public:
    inRange(T min_, T max_);
    bool operator()(T value) const;
private:
    T min, max;
};

template<typename T>
bool isAnyNearestNeighbourValue2D(ScalarField2D<T> field,
        plint iX, plint iY, T value);
template<typename T>
bool isAnyNearestNeighbourValue3D(ScalarField3D<T> field,
        plint iX, plint iY, plint iZ, T value);
template<typename T>
bool isAllNearestNeighbourValue2D(ScalarField2D<T> field,
        plint iX, plint iY, T value);
template<typename T>
bool isAllNearestNeighbourValue3D(ScalarField3D<T> field,
        plint iX, plint iY, plint iZ, T value);

template<typename T>
bool isAnyDirectNeighbourValue2D(ScalarField2D<T> field,
        plint iX, plint iY, T value);
template<typename T>
bool isAnyDirectNeighbourValue3D(ScalarField3D<T> field,
        plint iX, plint iY, plint iZ, T value);
template<typename T>
bool isAllDirectNeighbourValue2D(ScalarField2D<T> field,
        plint iX, plint iY, T value);
template<typename T>
bool isAllDirectNeighbourValue3D(ScalarField3D<T> field,
        plint iX, plint iY, plint iZ, T value);

template<typename T>
bool isAnyNextNeighbourValue3D(ScalarField3D<T> field,
        plint iX, plint iY, plint iZ, T value);

template<typename T>
void lessEqual(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<int>& result, Box3D domain);

template<typename T>
class A_le_alpha_functional3D : public BoxProcessingFunctional3D_SS<T,int>
{
public:
    A_le_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<int>& result);
    virtual A_le_alpha_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
void greaterEqual(MultiScalarField3D<T>& A, MultiScalarField3D<T>& B, MultiScalarField3D<int>& result, Box3D domain);

template<typename T>
class A_ge_alpha_functional3D : public BoxProcessingFunctional3D_SS<T,int>
{
public:
    A_ge_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<int>& result);
    virtual A_ge_alpha_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};


#endif
