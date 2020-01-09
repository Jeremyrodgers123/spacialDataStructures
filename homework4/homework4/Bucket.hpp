//
//  BucketGrid.hpp
//  homework4
//
//  Created by Jeremy Rodgers on 6/11/19.
//  Copyright © 2019 Jeremy Rodgers. All rights reserved.
//

#ifndef Bucket_hpp
#define Bucket_hpp
#include <iostream>
#include "Point.hpp"
#include <stdio.h>
#include <vector>

template<int Dimension>
class Bucket{
    public:
    std::vector<Point<Dimension>> points;
    Point<Dimension> size; //size (bucket length) for each dimension.
    
    Bucket(const Point<Dimension> &sizes){
        setBucketSize(sizes);
    }
    
    Bucket(std::vector<Point<Dimension>> &p){
        points = p;
        //SET SIZE?
    }
    
    void setBucketSize(const Point<Dimension> &_size){
        size = _size;
    }
    
    /**
     Check if point is in bucket.
     params: Point<Dimension> PointToCheck - the point we are searching for
     return: true if the bucket contains the point
     **/
    
    bool contains(const Point<Dimension> &pointToCheck){
        if(points.size() <= 0){
            return false;
        }
        for( Point<Dimension> p : points){
            bool isMatch = true;
            for(int dim = 0; dim < Dimension; dim++){
                if(pointToCheck[dim] != p[dim]){
                    isMatch = false;
                    break;
                }
            }
            if(isMatch){
                return true;
            }
        }
        return false;
    }
    
    void print(){
        std::cout << "__________" << std::endl;
        for (int i = 0; i < points.size(); i++){
            std::cout << points[i] << std::endl;
        }
        std::cout << "__________" << std::endl;
    }
};


#endif /* Bucket_hpp */
