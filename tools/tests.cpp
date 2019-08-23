#include "../externalLibraries/Catch2/single_include/catch2/catch.hpp"

#include <memory>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

// for creating necessary folder if it does not exist yet:
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// for current time and date
#include <time.h>

#include "../tools/tools.h"
#include "../tools/tools.hh"
#include "../tools/tools_geometric.h"
#include "../tools/tools_geometric.hh"

using namespace std;


const int iterations = 100;

TEST_CASE( "skip this test", "[!hide]" ) { /* ... */ }

TEST_CASE( "test_rand" ) {

    randNumber<double>   randomDouble(time(0));
    randNumber<int> randomInt(time(0));
    double randTMin = 0.;
    double randTMax = 10.;
    int randIntMin = 0;
    int randIntMax = 100;

    for( int i=0; i < iterations; i++ ) {
        REQUIRE( randomDouble(randTMin, randTMax) >= randTMin );
        REQUIRE( randomDouble(randTMin, randTMax) <= randTMax );
        REQUIRE( randomInt(randIntMin, randIntMax) >= randIntMin );
        REQUIRE( randomInt(randIntMin, randIntMax) <= randIntMax );
    }
}

TEST_CASE( "test_rampVelocity_function" ) {

    randNumber<double>     randomDouble(time(0));
    randNumber<int> randomInt(time(0));
    double velMax;
    int startT;
    int rampT;

    for( int i=0; i<iterations; i++ ) {
        velMax         = randomDouble(0.1, 2.);
        rampT          = randomInt(5, 5000);
        startT         = randomInt(0, rampT);

     // velocity should be velMax*startT/rampT if startT<rampT
        REQUIRE( rampVelocity(velMax, startT, rampT) == velMax*startT/rampT );
     // velocity should be velMax if startT>=rampT
        REQUIRE( rampVelocity(velMax, int(rampT*1.5), rampT) == velMax );
    }

}

TEST_CASE( "test_poiseuille_function" ) {

    randNumber<double>     randomDouble(time(0));
    randNumber<int> randomInt(time(0));
    double velMax;
    int velPow;
    int inletRadius;

    for( int i=0; i<iterations; i++ ) {
        velMax         = randomDouble(0.1, 2.);
     // should work for arbitrary pows
        velPow         = randomInt(1, 20);
        inletRadius = randomInt(1, 100);
     // should work for pos or neg velocity
        if (i > iterations/2) velMax *= -1;

     // velocity in middle should be max
        REQUIRE( poiseuilleVelocity(0, velMax, velPow, inletRadius) == velMax );
     // velocity at corner should be 0
        REQUIRE( poiseuilleVelocity(inletRadius, velMax, velPow, inletRadius) == 0 );
     // velocity outside inletRadius should be 0
        REQUIRE( poiseuilleVelocity(inletRadius+1, velMax, velPow, inletRadius) == 0 );
        REQUIRE( poiseuilleVelocity(inletRadius*5, velMax, velPow, inletRadius) == 0 );
    }
}


TEST_CASE( "isInCircle" ) {

    SECTION( "radius ==0, ==-1" ) {
        REQUIRE( isInCircle(0,0, 0, 0,0) == false );
        REQUIRE( isInCircle(0,0,-1, 0,0) == false );
    }

    SECTION( "centre ==origin" ) {
        int cx = 0; // x-coordinate of centre
        int cy = 0; // y-coordinate of centre
        int r = 10; // radius
        int rSq = r*r;
        for (int x=-r; x<=r; x++) {
            int yR = sqrt(rSq-x*x);
            for (int y=-yR; y<=yR; y++) {
                REQUIRE( isInCircle(x,y,rSq,cx,cy) );
                if (isInCircle(x,y,rSq,cx,cy))
                    REQUIRE( distanceFromPoint(x,y,cx,cy) <= r );
                else
                    REQUIRE( distanceFromPoint(x,y,cx,cy) > r );
                REQUIRE( isInCircle(x,yR+1,rSq,cx,cy) == false );
            }
        }
    }

    SECTION( "centre !=origin" ) {
        int cx = 2; // x-coordinate of centre
        int cy = 4; // y-coordinate of centre
        int r = 10; // radius
        int rSq = r*r;
        for (int x=-r; x<=r; x++) {
            int yR = sqrt(rSq-x*x);
            for (int y=-yR; y<=yR; y++) {
                if (isInCircle(x,y,rSq,cx,cy))
                    REQUIRE( distanceFromPoint(x,y,cx,cy) <= r );
                else
                    REQUIRE( distanceFromPoint(x,y,cx,cy) > r );
            }
        }
    }
}


TEST_CASE( "isInSphere" ) {

    SECTION( "radius ==0, ==-1" ) {
        REQUIRE( isInSphere(0,0,0, 0,0,0,0)  == false );
        REQUIRE( isInSphere(0,0,0, -1,0,0,0) == false );
    }

    SECTION( "outside" ) {
        REQUIRE( isInSphere(1,1,0, 1,0,0,0) == false );
        REQUIRE( isInSphere(0,1,1, 1,0,0,0) == false );
        REQUIRE( isInSphere(1,0,1, 1,0,0,0) == false );
    }
    SECTION( "onSurface" ) {
        REQUIRE( isInSphere(0,0,0, 1,0,0,0) == true );
        REQUIRE( isInSphere(1,0,0, 1,0,0,0) == true );
        REQUIRE( isInSphere(0,1,0, 1,0,0,0) == true );
        REQUIRE( isInSphere(0,0,1, 1,0,0,0) == true );
    }
    SECTION( "inside" ) {
        REQUIRE( isInSphere(0,0,0, 10,0,0,0) == true );
        REQUIRE( isInSphere(1,0,0, 10,0,0,0) == true );
        REQUIRE( isInSphere(0,1,0, 10,0,0,0) == true );
        REQUIRE( isInSphere(0,0,1, 10,0,0,0) == true );
    }

}


TEST_CASE( "isInTube" ) {

    SECTION( "height ==-1" ) {
        int x=0, y=0, z=0;
        int r=1, cX=0, cY=0, Z0=5, Z1=1;
        REQUIRE( isInTube(x,y,z, r,cX,cY,Z0,Z1) == false );
        REQUIRE( Z0 > Z1 );
    }

    SECTION( "underneath" ) {
        int x=0, y=0, z=0;
        int r=1, cX=0, cY=0, Z0=1, Z1=2;
        REQUIRE( isInTube(x,y,z, r,cX,cY,Z0,Z1) == false );
        REQUIRE( z < Z0 );
    }
    SECTION( "above" ) {
        int x=0, y=0, z=3;
        int r=1, cX=0, cY=0, Z0=1, Z1=2;
        REQUIRE( isInTube(x,y,z, r,cX,cY,Z0,Z1) == false );
        REQUIRE( z > Z1 );
    }
    SECTION( "inside" ) {
        int x=0, y=0;
        int r=1, cX=0, cY=0, Z0=1, Z1=20;
        int z=2;
        REQUIRE( isInCircle(x,y, r,cX,cY) );
        REQUIRE( z>Z0 );
        REQUIRE( z<Z1 );
        REQUIRE( isInTube(x,y,z, r,cX,cY,Z0,Z1) == true );
    }
}


TEST_CASE( "isInParabola" ) {

    SECTION( "outside: point lies below apex" ) {
        int x=0,  y=0;
        int cX=0, cY=1;
        int steepness=1;
        REQUIRE( isInParabola(x,y, steepness,cX,cY) == false );
    }
    SECTION( "outside: point lies above apex, but too far to the side" ) {
        int x=10,  y=1;
        int cX=0, cY=0;
        int steepness=1;
        REQUIRE( isInParabola(x,y, steepness,cX,cY) == false );
    }
    SECTION( "outside: shifted centre" ) {
        int x=0,  y=2;
        int cX=2, cY=0;
        int steepness=1;
        REQUIRE( isInParabola(x,y, steepness,cX,cY) == false );
    }
    SECTION( "outside: rotated, point lies above" ) {
        int x=1,  y=0;
        int cX=2, cY=0;
        int steepness=1;
        REQUIRE( isInParabola(y,x, steepness,cY,cX) == false );
    }
    SECTION( "inside: rotated, point lies above" ) {
        int x=1,  y=0;
        int cX=0, cY=0;
        int steepness=1;
        REQUIRE( isInParabola(y,x, steepness,cY,cX) == true );
    }
    SECTION( "inside: point lies above" ) {
        int x=0,  y=2;
        int cX=0, cY=0;
        int steepness=1;
        REQUIRE( isInParabola(x,y, steepness,cX,cY) == true );
    }
    SECTION( "inside: different steepness" ) {
        int x=1,  y=4;
        int cX=0, cY=0;
        int steepness=4;
        REQUIRE( isInParabola(x,y, steepness,cX,cY) == true );
    }
    SECTION( "inside: shifted centre" ) {
        int x=1,  y=1;
        int cX=1, cY=0;
        int steepness=1;
        REQUIRE( isInParabola(x,y, steepness,cX,cY) == true );
    }
    SECTION( "inside: inverted " ) {
        int x=0,  y=0;
        int cX=0, cY=1;
        int steepness=-1;
        REQUIRE( isInParabola(x,y, steepness,cX,cY) == true );
    }
    SECTION( "outside: rescaled" ) {
        int x=5,  y=4;
        int cX=0, cY=0;
        int steepness=1;
        int width = 4;
        double rescaled = double(steepness) / double (width);
        REQUIRE( isInParabola(x,y, rescaled,cX,cY) == false );
    }
    SECTION( "inside: rescaled" ) {
        int x=4,  y=4;
        int cX=0, cY=0;
        int steepness=1;
        int width = 4;
        double rescaled = double(steepness) / double (width);
        REQUIRE( isInParabola(x,y, rescaled,cX,cY) == true );
    }
}

TEST_CASE( "isInParaboloid" ) {
    SECTION( "inside: above apex" ) {
        int  x=0,  y=0,  z=1;
        int cX=0, cY=0, cZ=0;
        int steepness=1;
        REQUIRE( isInParaboloid(x,y,z, steepness,cX,cY,cZ) == true );
    }
    SECTION( "outside: reversed" ) {
        int  x=0,  y=0,  z=1;
        int cX=0, cY=0, cZ=0;
        int steepness=-1;
        REQUIRE( isInParaboloid(x,y,z, steepness,cX,cY,cZ) == false );
    }
    SECTION( "outside: rotated, below" ) {
        int  x=1,  y=0,  z=0;
        int cX=0, cY=0, cZ=0;
        int steepness=-1;
        REQUIRE( isInParaboloid(y,z,x, steepness,cY,cZ,cX) == false );
    }
    SECTION( "inside: rotated, above" ) {
        int  x=1,  y=0,  z=0;
        int cX=0, cY=0, cZ=0;
        int steepness=1;
        REQUIRE( isInParaboloid(y,z,x, steepness,cY,cZ,cX) == true );
    }
    SECTION( "outside: below" ) {
        int  x=0,  y=0,  z=-1;
        int cX=0, cY=0, cZ=0;
        int steepness=1;
        REQUIRE( isInParaboloid(x,y,z, steepness,cX,cY,cZ) == false );
    }
    SECTION( "outside: point lies above apex, but too far to the side" ) {
        int  x=10,  y=0,  z=1;
        int cX=0, cY=0, cZ=0;
        int steepness=1;
        REQUIRE( isInParaboloid(x,y,z, steepness,cX,cY,cZ) == false );
    }
    SECTION( "outside: shifted centre" ) {
        int  x=0,  y=0,  z=1;
        int cX=10, cY=0, cZ=0;
        int steepness=1;
        REQUIRE( isInParaboloid(x,y,z, steepness,cX,cY,cZ) == false );
    }
    SECTION( "inside: different steepness" ) {
        int  x=1,  y=0,  z=4;
        int cX=0, cY=0, cZ=0;
        int steepness=4;
        REQUIRE( isInParaboloid(x,y,z, steepness,cX,cY,cZ) == true );
    }
    SECTION( "outside: rescaled" ) {
        int  x=4,  y=1,  z=4;
        int cX=0, cY=0, cZ=0;
        int steepness=1;
        int width = 4;
        double rescaled = double(steepness) / double (width);
        REQUIRE( isInParaboloid(x,y,z, rescaled,cX,cY,cZ) == false );
    }
    SECTION( "inside: rescaled" ) {
        int  x=4,  y=0,  z=4;
        int cX=0, cY=0, cZ=0;
        int steepness=1;
        int width = 4;
        double rescaled = double(steepness) / double (width);
        REQUIRE( isInParaboloid(x,y,z, rescaled,cX,cY,cZ) == true );
    }

}

TEST_CASE( "dotProduct" ) {
    SECTION( "float parallel" ) {
        float aX=0, aY=0, aZ=1.1;
        float bX=0, bY=0, bZ=10;
        REQUIRE( dotProduct(aX,aY,aZ, bX,bY,bZ) == 11. );
    }
    SECTION( "float perpendicular" ) {
        float aX=1, aY=1,  aZ=0;
        float bX=1, bY=-1, bZ=0;
        REQUIRE( dotProduct(aX,aY,aZ, bX,bY,bZ) == 0. );
    }
    SECTION( "float perpendicular (simple)" ) {
        float aX=0, aY=0, aZ=1.1;
        float bX=0, bY=.5, bZ=0;
        REQUIRE( dotProduct(aX,aY,aZ, bX,bY,bZ) == 0. );
    }
    SECTION( "integer perpendicular" ) {
        int aX=0, aY=0, aZ=1;
        int bX=0, bY=1, bZ=0;
        REQUIRE( dotProduct(aX,aY,aZ, bX,bY,bZ) == 0. );
    }
}

TEST_CASE( "normalizeVector" ) {
    SECTION( "float" ) {
        float x=0, y=0, z=1;
        std::vector<double> normedVec = {0,0,1};
        REQUIRE( normalizeVector(x,y,z) == normedVec );
    }
    SECTION( "int" ) {
        int x=0, y=0, z=1;
        std::vector<double> normedVec = {0,0,1};
        REQUIRE( normalizeVector(x,y,z) == normedVec );
    }
    SECTION( "int short length" ) {
        int x=0, y=1, z=1;
        double value = 1/sqrt(2);
        std::vector<double> normedVec = {0,value,value};
        REQUIRE( normalizeVector(x,y,z) == normedVec );
    }
    SECTION( "int long length" ) {
        int x=0, y=100, z=100;
        double value = 1/sqrt(2);
        std::vector<double> normedVec = {0,value,value};
        REQUIRE( normalizeVector(x,y,z) == normedVec );
    }
    SECTION( "int negative" ) {
        int x=-1, y=1, z=0;
        double value = 1/sqrt(2);
        std::vector<double> normedVec = {-value,value,0};
        REQUIRE( normalizeVector(x,y,z) == normedVec );
    }
}

TEST_CASE( "isBelowPlane" ) {
    SECTION( "below xy" ) {
        int nX=0, nY=0, nZ=1;
        double d = 0.;
        int  x=1,  y=1,  z=-1;
        REQUIRE( isBelowPlane(x,y,z, nX,nY,nZ, d) == true );
    }
    SECTION( "above xy" ) {
        int nX=0, nY=0, nZ=1;
        double d = 0.;
        int  x=1,  y=1,  z=1;
        REQUIRE( isBelowPlane(x,y,z, nX,nY,nZ, d) == false );
    }
    SECTION( "below arbitrary" ) {
        double nX=0, nY=sqrt(.5), nZ=sqrt(.5);
        double d = 1.;
        double  x=0, y=0, z=0;
        REQUIRE( isBelowPlane(x,y,z, nX,nY,nZ, d) == true );
    }
    SECTION( "below negative" ) {
        double nX=-sqrt(.5), nY=sqrt(.5), nZ=0;
        double d = 1.;
        double  x=1, y=0, z=0;
        REQUIRE( isBelowPlane(x,y,z, nX,nY,nZ, d) == true );
    }
}


TEST_CASE( "sign function" ) {

    SECTION( "positive number" ) {
        int x = 1;
        REQUIRE( 1 == sign(x) );
    }
    SECTION( "negative number" ) {
        int x = -1;
        REQUIRE( -1 == sign(x) );
    }
    SECTION( "zero" ) {
        int x = 0;
        REQUIRE( 0 == sign(x) );
    }
}

TEST_CASE( "gradient calculation" ) {

    SECTION( "slope of 1" ) {
        double x1=0, x2=1, step=1;
        double gradient = getGradient(x1, x2, step);

        REQUIRE( gradient == 1 );
    }

    SECTION( "slope of 2" ) {
        double x1=0, x2=2, step=1;
        double gradient = getGradient(x1, x2, step);

        REQUIRE( gradient == 2. );

        step=0.5;
        x2=1;
        gradient = getGradient(x1, x2, step);

        REQUIRE( gradient == 2. );
    }
}

TEST_CASE( "read last line from file" ) {
    SECTION( "short line" ) {
        std::string fileName   = "test_file.txt";
        std::string firstLine  = "first line of text";
        std::string lastLine   = "1. 2. 3.5 3333.223";

        std::stringstream text;
        text << firstLine << endl
             << lastLine  << endl;

        std::ofstream ofile;
        ofile.open(fileName.c_str(), std::ios_base::out | std::ios_base::trunc );
        ofile << text.str();
        ofile.close();

        std::string result = readLastLineInFile(fileName);
        REQUIRE( result == lastLine );
    }
    SECTION( "long line" ) {
        std::string fileName   = "test_file.txt";
        std::string firstLine  = "first line of text";
        std::string lastLine   = "1. 2. 3.5 3333.223";

        std::stringstream text;
        text << firstLine << endl
             << lastLine << lastLine
             << lastLine << lastLine
             << lastLine << lastLine
             << lastLine << lastLine
             << lastLine << lastLine
             << lastLine << lastLine
             << lastLine << lastLine
             << lastLine << lastLine
             << lastLine << lastLine
             << endl;

        std::ofstream ofile;
        ofile.open(fileName.c_str(), std::ios_base::out | std::ios_base::trunc );
        ofile << text.str();
        ofile.close();

        std::string result = readLastLineInFile(fileName);
        REQUIRE( result == lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine+lastLine );
    }
}

TEST_CASE( "convertStringToDouble" ) {
    std::vector<string> stringVector = {"1.", "2."};
    std::vector<double> doubleVector = {1., 2.};

    std::vector<double> result = convertStringToDouble(stringVector);
    REQUIRE( result == doubleVector );
}

TEST_CASE( "splitString" ) {
    std::string dataAsLine = "1. 2.2 3.33 44.44";
    std::vector<std::string> dataAsVector = {"1.", "2.2", "3.33", "44.44"};
    std::vector<string> result = splitString(dataAsLine);
    REQUIRE( result == dataAsVector );

}
TEST_CASE( "read elements of last line") {

    std::string fileName   = "test_file.txt";
    std::string firstLine  = "first line of text";
    std::string lastLine   = "1. 2. 3.5 3333.223";
    std::vector<double> lastLineDouble = {1., 2., 3.5, 3333.223};

    std::stringstream text;
    text << firstLine << endl
         << lastLine  << endl;

    std::ofstream ofile;
    ofile.open(fileName.c_str(), std::ios_base::out | std::ios_base::trunc );
    ofile << text.str();
    ofile.close();

    std::string lastLineTmp = readLastLineInFile(fileName);
    std::vector<string> lastLineElemtents = splitString(lastLineTmp);
    std::vector<double> recentData = convertStringToDouble(lastLineElemtents);

    REQUIRE( recentData == lastLineDouble );
}

TEST_CASE( "changeElementsInMap" ) {
    std::unordered_map<std::string, double> dataMap = {
      { "one", 1. },
      { "two", 2. }};

    std::unordered_map<std::string, double> updatedMap = {
      { "one", 3. },
      { "two", 4. }};

    std::vector<string> keys = {"one", "two"};
    std::vector<double> newData = {3., 4.};
    changeElementsInMap(dataMap, keys, newData);

    REQUIRE( dataMap["one"] ==  updatedMap["one"] );
    REQUIRE( dataMap["two"] ==  updatedMap["two"] );
}

TEST_CASE( "createStringFromVector" ) {

    std::vector<string> keys = {"one", "two"};
    std::string result = createStringFromVector(keys, " ");

    REQUIRE( result == "one two\n" );
}

TEST_CASE( "createStringFromMap" ) {
    std::unordered_map<std::string, double> dataMap;
    dataMap = {
      { "one", 1. },
      { "two", 2. }};

    std::vector<string> keys = {"one", "two"};
    std::string secondRow = createStringFromMap(dataMap, keys, " ");

    REQUIRE( secondRow == "1 2\n" );
}

TEST_CASE( "stationary") {
    double threshold = 0.1;
    double currentValue = 0.01;
    int characteristicCount = 10;

    belowThresholdCounter testChecker(threshold, characteristicCount);

    SECTION( "not at threshold long enough" ) {
        for (int i=1; i<characteristicCount; i++)
            REQUIRE( testChecker.repeatedlyBelowThrehold(currentValue) == false );
    }

    SECTION( "at threshold long enough" ) {
        for (int i=1; i<characteristicCount; i++)
            REQUIRE( testChecker.repeatedlyBelowThrehold(currentValue) == false );
        REQUIRE( testChecker.repeatedlyBelowThrehold(currentValue) == true );
    }

    SECTION( "negative values" ) {
        for (int i=1; i<characteristicCount; i++)
            REQUIRE( testChecker.repeatedlyBelowThrehold(currentValue) == false );
        REQUIRE( testChecker.repeatedlyBelowThrehold(-1*currentValue) == true );
    }

    SECTION( "above threshold" ) {
        for (int i=1; i<characteristicCount+1; i++)
            REQUIRE( testChecker.repeatedlyBelowThrehold(currentValue+1.) == false );
    }

}

TEST_CASE( "calculateCOM for single point" ) {

    std::vector<int> point{0,0,0};
    std::vector<std::vector<int>> voxels{ point };

    std::vector<float> center;
    center = calculateCOM<int,float>(voxels);

    REQUIRE( center.at(0) == (float) point.at(0) );
    REQUIRE( center.at(1) == (float) point.at(1) );
    REQUIRE( center.at(2) == (float) point.at(2) );
}

TEST_CASE( "calculateCOM for multiple points" ) {

    std::vector<std::vector<int>> voxels;
    for (int i=0; i<5; i++) {
        std::vector<int> point{i,i,i};
        voxels.push_back(point);
    }

    std::vector<float> center;
    center = calculateCOM<int,float>(voxels);
    std::vector<float> known_center{2,2,2};

    REQUIRE( center.at(0) == (float) known_center.at(0) );
    REQUIRE( center.at(1) == (float) known_center.at(1) );
    REQUIRE( center.at(2) == (float) known_center.at(2) );
}

TEST_CASE( "calculate diameterFromSphereVolume", "[!hide]" ) {

    SECTION( "volume=1" ) {
        int volume = 1;
        float diameter = diameterFromSphereVolume<int,float>(volume);
        float known_diameter = 1.2407;

        REQUIRE( diameter == known_diameter );
    }
}

