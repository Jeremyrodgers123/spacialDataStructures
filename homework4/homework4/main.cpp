//
//  main.cpp
//  homework4
//
//  Created by Jeremy Rodgers on 6/11/19.
//  Copyright © 2019 Jeremy Rodgers. All rights reserved.
//

#include <iostream>
#include <array>
#include <vector>
#include "Point.hpp"
#include "Bucket.hpp"
#include "BucketKNN.hpp"
#include "Generators.hpp"
#include "KDTree.hpp"
#include "QuadTree.hpp"
#include "csvBuilder.hpp"
#include <fstream>
#include "dumbKNN.hpp"
using namespace std;

template <int Dimension>
void printPoints(std::vector<Point<Dimension>> points){
    for( Point<2> p : points ){
        std:: cout << p << std::endl;
    }
    std::cout << std::endl;
}

void printResults(std::vector<std::vector<Point<2>>> results){
    for(std::vector<Point<2>> resultSet: results){
        printPoints(resultSet);
    }
}
void buildGaussian(long numPoints, long maxPoints, int n, std::ofstream &os){
    std::vector<Point<2>> gaussianPoints2D;
    std::vector<Point<3>> gaussianPoints3D;
    std::vector<Point<10>> gaussianPoints10D;
    
    while( numPoints <= maxPoints){
        GaussianGenerator<2> gaussian2D = GaussianGenerator<2>(-numPoints,+numPoints);
        GaussianGenerator<3> gaussian3D = GaussianGenerator<3>(-numPoints,+numPoints);
        GaussianGenerator<10> gaussian10D = GaussianGenerator<10>(-numPoints,+numPoints);
        
        for( int i = 0; i < numPoints; i++){
            Point<2> p2D = gaussian2D.generatePoint();
            Point<3> p3D = gaussian3D.generatePoint();
            Point<10> p10D = gaussian10D.generatePoint();
            gaussianPoints2D.push_back(p2D);
            gaussianPoints3D.push_back(p3D);
            gaussianPoints10D.push_back(p10D);
        }
        
        // put datastructures here!
        //DUMB KNN
        DumbKNN<2> dumbKNN2D = DumbKNN<2>(gaussianPoints2D);
        buildCSV<2>( gaussian2D, dumbKNN2D, 25, "gaussian", "Dumb", os);
        
        //BUCKETS
        BucketKNN<2> buckets = BucketKNN<2>(gaussianPoints2D, 5);
        buildCSV<2>( gaussian2D, buckets, 25, "gaussian", "Buckets", os);
        
        BucketKNN<3> buckets3D = BucketKNN<3>(gaussianPoints3D, 5);
        buildCSV<3>( gaussian3D, buckets3D, 25, "gaussian", "Buckets", os);
        
        BucketKNN<10> buckets10D = BucketKNN<10>(gaussianPoints10D, 5);
        buildCSV<10>( gaussian10D, buckets10D, 25, "gaussian", "Buckets", os);
        //QUAD TREE
        QuadTree<2> gaussianQT = QuadTree<2>(gaussianPoints2D, 5);
        buildCSV<2>( gaussian2D, gaussianQT, 25, "gaussian", "Quad Tree", os);
        
        //KD TREE
        KDTree<2> kdTree = KDTree<2>(gaussianPoints2D);
        buildCSV<2>( gaussian2D, kdTree, 25, "gaussian", "KD Tree", os);
        
        KDTree<3> kdTree3D = KDTree<3>(gaussianPoints3D);
        buildCSV<3>( gaussian3D, kdTree3D, 25, "gaussian", "KD Tree", os);
        
        KDTree<10> kdTree10D = KDTree<10>(gaussianPoints10D);
        buildCSV<10>( gaussian10D, kdTree10D, 25, "gaussian", "KD Tree", os);
        
        n++;
        numPoints = pow(2,n);
        std::cout << "n: " << n << std::endl;
        gaussianPoints2D.clear();
        gaussianPoints3D.clear();
        gaussianPoints10D.clear();
    }
}

void buildUniform(long numPoints, long maxPoints, int n, std::ofstream &os){
    std::vector<Point<2>> uniformPoints2D;
    std::vector<Point<3>> uniformPoints3D;
    std::vector<Point<10>> uniformPoints10D;
    
    while( numPoints <= maxPoints){
        UniformGenerator<2> uniform2D = UniformGenerator<2>(-numPoints,+numPoints);
        UniformGenerator<3> uniform3D = UniformGenerator<3>(-numPoints,+numPoints);
        UniformGenerator<10> uniform10D = UniformGenerator<10>(-numPoints,+numPoints);
        
        for( int i = 0; i < numPoints; i++){
            Point<2> p2D = uniform2D.generatePoint();
            Point<3> p3D = uniform3D.generatePoint();
            Point<10> p10D = uniform10D.generatePoint();
            
            uniformPoints2D.push_back(p2D);
            uniformPoints3D.push_back(p3D);
            uniformPoints10D.push_back(p10D);
        }
        //DUMB KNN
        DumbKNN<2> dumbKNN2D = DumbKNN<2>(uniformPoints2D);
        buildCSV<2>( uniform2D, dumbKNN2D, 25, "uniform", "Dumb", os);
        
        //BUCKETS
        BucketKNN<2> buckets = BucketKNN<2>(uniformPoints2D, 5);
        buildCSV<2>( uniform2D, buckets, 25, "uniform", "Buckets", os);
        
        BucketKNN<3> buckets3D = BucketKNN<3>(uniformPoints3D, 5);
        buildCSV<3>( uniform3D, buckets3D, 25, "uniform", "Buckets", os);
        
        BucketKNN<10> buckets10D = BucketKNN<10>(uniformPoints10D, 5);
        buildCSV<10>( uniform10D, buckets10D, 25, "uniform", "Buckets", os);
        //Quad TREE
        QuadTree<2> uniformQT = QuadTree<2>(uniformPoints2D, 5);
        buildCSV<2>( uniform2D, uniformQT, 25, "uniform", "Quad Tree", os);
        
        //KD TREE
        KDTree<2> kdTree = KDTree<2>(uniformPoints2D);
        buildCSV<2>( uniform2D, kdTree, 25, "uniform", "KD Tree", os);
        
        KDTree<3> kdTree3D = KDTree<3>(uniformPoints3D);
        buildCSV<3>( uniform3D, kdTree3D, 25, "uniform", "KD Tree", os);
        
        KDTree<10> kdTree10D = KDTree<10>(uniformPoints10D);
        buildCSV<10>( uniform10D, kdTree10D, 25, "uniform", "KD Tree", os);
        n++;
        numPoints = pow(2,n);
        std::cout << "n: " << n << std::endl;
        uniformPoints2D.clear();
        uniformPoints3D.clear();
        uniformPoints10D.clear();
    }
}


void test(int numPoints,std::ofstream &os){
     UniformGenerator<2> uniform2D = UniformGenerator<2>(-numPoints,+numPoints);
     std::vector<Point<2>> uniformPoints2D;
    for( int i = 0; i < numPoints; i++){
        Point<2> p2D = uniform2D.generatePoint();
        uniformPoints2D.push_back(p2D);
    }
    DumbKNN<2> dumbKNN2D = DumbKNN<2>(uniformPoints2D);
    buildCSV<2>( uniform2D, dumbKNN2D, 3, "uniform", "Dumb", os);
}

int main(int argc, const char * argv[]) {
    std::ofstream os;
   
    std::cout << argc << std::endl;
    if(argc != 2){
        std::cout << "only 1 argument, file name needed\n";
        return 0;
    }
    os.open(argv[1], ios::out);
    if(os.fail()){
        throw "file open failed";
    }
    os << "Type," << "Distribution,"<< "K," << "N," << "Dimensions," << "Time" << "\n";
    
    int start = 10;
    long minPoints = pow(2,start);
    long maxPoints = pow(2,21);
    
//    test(minPoints, os);
    buildGaussian(minPoints, maxPoints, start, os);
    buildUniform(minPoints, maxPoints, start, os);
    
    os.close();
    return 0;
}
