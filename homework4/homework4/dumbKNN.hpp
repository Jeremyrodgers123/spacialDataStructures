//
//  dumbKNN.hpp
//  homework4
//
//  Created by Jeremy Rodgers on 6/28/19.
//  Copyright © 2019 Jeremy Rodgers. All rights reserved.
//

#pragma once
#include <stdio.h>
#include <vector>
#include "Point.hpp"
#include <algorithm>
template <int Dimension>
class DumbKNN : public DataStructure<Dimension>{
public:
    std::vector<Point<Dimension>> points;
    AABB<Dimension> aabb;
    DumbKNN(std::vector<Point<Dimension>>& _points){
//         std::sort<Point<Dimension>, DistanceComparator>(_points.begin(), _points.end(), DistanceComparator<Dimension>(center));
        points = _points;
        aabb = getBounds(points);
    }
    
    std::vector<Point<Dimension>> KNN(Point<Dimension>& center, int k){
        std::vector<Point<Dimension>> ret;
        sort(center);
        for(int i = 0; i < k; i++){
            ret.push_back(points[i]);
        }
        return ret;
    }
    
    void sort(Point<Dimension>& center){
        std::sort<Point<Dimension>, DistanceComparator<Dimension>>(points.begin(), points.end(), DistanceComparator<Dimension>(center));
        aabb = getBounds(points);
    }
    
    int getSize(){
        return points.size();
    };
    
    AABB<Dimension> getBoxBounds(){
        return aabb;
    };
};
