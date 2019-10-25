#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace std;
typedef double T;


int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    const int N = 10;
    ImageWriter<T> imageWriter("leeloo");

    MultiScalarField3D<T> xCoord(N+1,N+1,N+1);
    MultiScalarField3D<T> yCoord(N+1,N+1,N+1);
    MultiScalarField3D<T> zCoord(N+1,N+1,N+1);

    addInPlace(xCoord,T(1));
    addInPlace(yCoord,T(2));
    addInPlace(zCoord,T(3));

    MultiTensorField3D<T,3> tensor(N+1,N+1,N+1);
    insertComponent(tensor, xCoord, 0);
    insertComponent(tensor, yCoord, 1);
    insertComponent(tensor, zCoord, 2);

    //PLB_ASSERT( tensor.get(0,0,0)[0] == T(1) );
    //PLB_ASSERT( tensor.get(0,0,0)[1] == T(2) );
    //PLB_ASSERT( tensor.get(0,0,0)[2] == T(3) );

    pcout << tensor.get(0,0,0)[0] << tensor.get(0,0,0)[1] << tensor.get(0,0,0)[2] << endl;

    imageWriter.writeGif("X", *extractComponent(*extractSubDomain(tensor, Box3D(0,N, 0,N, N/2,N/2)), 0), 0,3);
    imageWriter.writeGif("Y", *extractComponent(*extractSubDomain(tensor, Box3D(0,N, 0,N, N/2,N/2)), 1), 0,3);
    imageWriter.writeGif("Z", *extractComponent(*extractSubDomain(tensor, Box3D(0,N, 0,N, N/2,N/2)), 2), 0,3);
}
