//
//  csvBuilder.hpp
//  homework4
//
//  Created by Jeremy Rodgers on 6/24/19.
//  Copyright © 2019 Jeremy Rodgers. All rights reserved.
//

#pragma once
#include "BucketKNN.hpp"
#include "KDTree.hpp"
#include <iostream>
#include <fstream>
#include "Generators.hpp"
#include "Stopwatch.hpp"
#include <vector>
#include "dumbKNN.hpp"

double getAvgTime(std::vector<double> &queryTimes){
    double sum = 0;
    for(int i = 0; i < queryTimes.size(); i++){
        sum += queryTimes[i];
    }
    return sum/ queryTimes.size();
}

template <int Dimension>
void clampToBounds(AABB<Dimension> &bounds, Point<Dimension> &center){
    for( int dim = 0; dim < Dimension; dim++){
        if(center[dim] > bounds.maxs[dim]){
            center.point[dim] = bounds.maxs[dim];
        }else if(center[dim] < bounds.mins[dim]){
            center.point[dim] = bounds.mins[dim];
        }
    }
}

template<int Dimension>
void buildCSV(Generator<Dimension> &generator, DataStructure<Dimension> &dataStructure, int kIncrement, string distroType, string structureType, std::ofstream &os){
//    int numQueries = 1;
//     std::cout << "Start" << std::endl;
    Stopwatch sw;
    
    int max = std::min(dataStructure.getSize()/2, 200);
    for(int k = 10; k <= max; k = k+= kIncrement){
        
        std::vector<double> queryTimes;
        // pick 10 different centers
        for(int i = 0; i < 3; i++){
            Point<Dimension> center = generator.generatePoint();
            AABB<Dimension> bounds = dataStructure.getBoxBounds();
            clampToBounds(bounds, center);
            sw.start();
            std::vector<Point<Dimension>> result = dataStructure.KNN(center, k);
            double queryTime = sw.stop();
            queryTimes.push_back(queryTime);
        }
        double avgQueryTime = getAvgTime(queryTimes);
        os << structureType << "," << distroType <<"," << k << "," << dataStructure.getSize() << "," << Dimension << "," << avgQueryTime << "\n";
//        std::cout << k <<std::endl;
        if(k == 10){
            k=0;
        }
//        std::cout << "k:" << k<< std::endl;
    }
//    std::cout << "Done" << std::endl;
};


template<int Dimension>
void buildCSV(Generator<Dimension> &generator, KDTree<Dimension> &dataStructure,  int kIncrement, string distroType, string structureType, std::ofstream &os){
    //    int numQueries = 1;
    //     std::cout << "Start" << std::endl;
    Stopwatch sw;
    
    int max = std::min(dataStructure.getSize(), 200);
    for(int k = 10; k <= max; k += kIncrement){
        
        std::vector<double> queryTimes;
        // pick 10 different centers
        for(int i = 0; i < 3; i++){
            Point<Dimension> center = generator.generatePoint();
            AABB<Dimension> bounds = dataStructure.getBoxBounds();
            clampToBounds(bounds, center);
            sw.start();
            std::vector<Point<Dimension>> result = dataStructure.KNN(center, k);
            double queryTime = sw.stop();
            queryTimes.push_back(queryTime);
        }
        double avgQueryTime = getAvgTime(queryTimes);
        os << structureType << "," << distroType <<"," << k << "," << dataStructure.getSize() << "," << Dimension << "," << avgQueryTime << "\n";
        if(k == 10){
            k=0;
        }
        //        std::cout << "j:" << j << std::endl;
    }
    //    std::cout << "Done" << std::endl;
};


