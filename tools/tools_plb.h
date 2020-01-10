#ifndef TOOLS_PLB_H
#define TOOLS_PLB_H

#include <../src/core/globalDefs.h>
#include <../src/io/parallelIO.h>
#include <../externalLibraries/tinyxml/tinyxml.h>

#include <string>
#include <vector>
#include <random>


using namespace std;
using namespace plb;


void renameAllExisting(string outDir, string fileName);

bool copyFile(string sourceName, string destinationName);

void rewriteOutputFile(std::string data, std::string fileName);

void appendOutputFile(std::string data, std::string fileName);


Box3D layer(plint Nx, plint Ny, plint i) {
    // create slice with width=Nx, length=Ny, height=1
    return Box3D(0,Nx-1, 0,Ny-1, i,i);
}

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

/* ******** MaskedReplaceAlphaFunctional3D ************************************** */
// based on: MaskedBoxScalarMaxFunctional3D : public ReductiveBoxProcessingFunctional3D_SS<T,int>
// in: dataAnalysisFunctional3D.h:512
template<typename T1, typename T2>
class MaskedReplaceAlphaFunctional3D : public plb::ReductiveBoxProcessingFunctional3D_SS<T1, T2>
{
public:
    MaskedReplaceAlphaFunctional3D(T1 alpha_, std::vector<int>& flags_);
    MaskedReplaceAlphaFunctional3D(T1 alpha_, int flag_);
    virtual void process( Box3D domain,
                          ScalarField3D<T1>& scalarField,
                          ScalarField3D<T2>& mask );
    virtual MaskedReplaceAlphaFunctional3D<T1,T2>* clone() const {
        return new MaskedReplaceAlphaFunctional3D<T1,T2>(*this); }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
private:
    T1 alpha;
    std::vector<int> flags;
};

/* ******** MaskedBoxScalarListFunctional3D ******* */
template<typename T, class BoolMask>
class MaskedBoxScalarListFunctional3D : public ReductiveBoxProcessingFunctional3D_SS<T,int> {
public:
    MaskedBoxScalarListFunctional3D(BoolMask condition_)
        : listScalarId(this->getStatistics().subscribeList()),
          condition(condition_) { }
    void process ( Box3D domain,
                   ScalarField3D<T>& scalarField,
                   ScalarField3D<int>& mask );
    MaskedBoxScalarListFunctional3D<T,BoolMask>* clone() const {
        return new MaskedBoxScalarListFunctional3D<T,BoolMask>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    std::vector<T> getListScalar() const;
private:
    plint listScalarId;
    BoolMask condition;
};

/* *************** PseudomaskedSmoothen3D ******************************************* */

template<typename T>
class PseudomaskedSmoothen3D : public BoxProcessingFunctional3D_SS<T,T> {
public:
    PseudomaskedSmoothen3D(T weightDirectNeighbour_)
        : weightDirectNeighbour(weightDirectNeighbour_),
          weightCentreCell(1.) { }
    virtual void process(Box3D domain, ScalarField3D<T>& data, ScalarField3D<T>& result);
    virtual PseudomaskedSmoothen3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T weightDirectNeighbour;
    T weightCentreCell;
};


/* *************** PseudomaskedSmoothen3D ******************************************* */
// This is a "bulk" version. It has no special boundary treatment.

template<typename T>
void pseudomaskedSmoothen(MultiScalarField3D<T>& data, MultiScalarField3D<T>& result, Box3D domain, T weightDirectNeighbour);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > pseudomaskedSmoothen(MultiScalarField3D<T>& data, Box3D domain, T weightDirectNeighbour);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > pseudomaskedSmoothen(MultiScalarField3D<T>& data, T weightDirectNeighbour);

template<typename T1, typename T2>
class PoiseuilleVelocity {
    public:
        PoiseuilleVelocity(T2 uLB_, T2 uDev_, char dir_,
                T1 inletCentre_, T1 InR);
        PoiseuilleVelocity(T2 uLB_, T2 uDev_, char dir_,
                T1 inletCentreA_, T1 inletCentreB_, T1 InR);
        PoiseuilleVelocity(T2 uLB_, T2 uDev_, char dir_,
                T1 inletCentreA_, T1 inletCentreB_, T1 InR, T1 InR2_);
        void operator()(T1 iX, T1 iY, Array<T2,2>& u) const;
        void operator()(T1 iX, T1 iY, T1 iZ, Array<T2,3>& u) const;
    private:
        T2 uLB;
        T2 uDev;
        char dir;
        T1 inletCentreA;
        T1 inletCentreB;
        T1 InRSq;
        T1 InR2Sq;
        unsigned seed;
        std::default_random_engine* generator;
        std::normal_distribution<>* main_distribution;
        std::normal_distribution<>* side_distribution;
};

 template<typename T>
 class OwnUpdateAveScalarTransientStatistics3D : public BoxProcessingFunctional3D_SS<T,T>
 {
 public:
     OwnUpdateAveScalarTransientStatistics3D(plint n_);
     virtual void process(Box3D domain, ScalarField3D<T> &scalar, ScalarField3D<T> &avg);
     virtual OwnUpdateAveScalarTransientStatistics3D<T>* clone() const;
     virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
     virtual BlockDomain::DomainT appliesTo() const;
 private:
     plint n;
 };

template<typename T>
ScalarField2D<T> extractMiddleLayer(MultiScalarField3D<T>& field3D, int dir=2);


#endif
