//
//  QuadTree.hpp
//  homework4
//
//  Created by Jeremy Rodgers on 6/21/19.
//  Copyright © 2019 Jeremy Rodgers. All rights reserved.
//

#ifndef QuadTree_hpp
#define QuadTree_hpp

#include <stdio.h>
#include <vector>
#include <memory>
#include "Point.hpp"
#include <queue>

using namespace std;
template<int Dimension>
struct Node{
    bool isLeaf = false;
    std::vector<Point<Dimension>> points;
    std::unique_ptr<Node<Dimension>> NW, NE, SE,SW;
    AABB<Dimension> aabb;
    
    Node(typename std::vector<Point<Dimension>>::iterator begin, typename std::vector<Point<Dimension>>::iterator end, int maxBucketSize){
        typename std::vector<Point<Dimension>>::iterator i;
        aabb = getBounds<Dimension>(begin, end);
        int currentSize = end - begin;
        if(currentSize <= maxBucketSize){ //create leaf node
            isLeaf = true;
            for(i = begin; i != end; ++i){
                points.push_back(*i);
            }
        }else{
            Point<Dimension> middle = splitBox(aabb, 2);
            splitOnDimensions(begin, end, middle, 0);
            auto Q1 = partition_point(begin, end, [&middle](Point<Dimension> p){
                for(int dim = 0; dim < Dimension; dim++){
                    if(p[dim] < middle[dim]) return false;
                }
                return true;
            });
            auto Q2 = partition_point(Q1, end, [&middle](Point<Dimension> &p){
                if( (p[0] >= middle[0]) && (p[1] < middle[1])){
                    return true;
                }
                return false;
            });
            auto Q3 = partition_point(Q2, end, [&middle](Point<Dimension> p){
                if(p[0] < middle[0] && p[1] >= middle[1]){
                    return true;
                }
                return false;
            });
            NE = std::make_unique<Node<Dimension>> (begin, Q1, maxBucketSize);
            SE = std::make_unique<Node<Dimension>> (Q1, Q2, maxBucketSize);
            NW = std::make_unique<Node<Dimension>> (Q2, Q3, maxBucketSize);
            SW = std::make_unique<Node<Dimension>> (Q3, end, maxBucketSize);
        }
        
    }
    
    void splitOnDimensions(typename std::vector<Point<Dimension>>::iterator begin,
                           typename std::vector<Point<Dimension>>::iterator end,
                           Point<Dimension> middle, int dim){
        if(dim >= Dimension){
            return;
        }
        auto groupTwoIterator = std::partition(begin, end, [&middle, &dim](Point<Dimension> p)
        {
            return p[dim] >= middle[dim];
        });
        dim = dim + 1;
        splitOnDimensions(begin, groupTwoIterator, middle, dim);
        splitOnDimensions(groupTwoIterator, end, middle, dim);
    }
    
    void printNodes(){
        if(isLeaf){
            for(int i = 0; i < points.size(); i++){
                std::cout << points[i] << std::endl;
            }
            return;
        }
        NW->printNodes();
        NE->printNodes();
        SE->printNodes();
        SW->printNodes();
    }
    
    void nearRadius(Point<Dimension> center, float radius, std::vector<Node<Dimension>*> &ret){
        Point<Dimension> closestPoint = aabb.closestInBox(center);
        float d =0;
        d = distance(center, closestPoint);
        if(d > radius){
            return;
        }else if(isLeaf){
            ret.push_back(this);
            return;
        }else{
            NW->nearRadius(center, radius, ret);
            NE->nearRadius(center, radius, ret);
            SE->nearRadius(center, radius, ret);
            SW->nearRadius(center, radius, ret);
        }
    }
    
    void nearestNeighbors(Point<Dimension> center,
                          int k,
                          float radius,
                          std::priority_queue<Point<Dimension>, std::vector<Point<Dimension>>, DistanceComparator<Dimension>> &pq){
        Point<Dimension> closest = aabb.closestInBox(center);
        float dist = distance(closest, center);
        if((dist > radius) ){ //&& (pq.size() >= k)
            return;
        }
        if(isLeaf){
            for(Point<Dimension> p : points){
                if(pq.size() < k){
                    pq.push(p);
                }else{
                    float currentPDist = distance(p, center);
                    float pqFurthestDist = distance(pq.top(), center);
                    if(currentPDist < pqFurthestDist ){
                        pq.pop();
                        pq.push(p);
                        radius = distance(center, pq.top());
                    }
                }
            }
        }else{
            NW->nearestNeighbors( center, k, radius, pq);
            NE->nearestNeighbors( center, k, radius, pq);
            SE->nearestNeighbors( center, k, radius, pq);
            SW->nearestNeighbors( center, k, radius, pq);
        }
    }

    void printPoints(typename std::vector<Point<Dimension>>::iterator begin, typename std::vector<Point<Dimension>>::iterator end){
        std::cout << "start of points" <<std::endl;
        for(auto i = begin; i != end; i++){
            std::cout << *i << std::endl;
        }
        std::cout << std::endl;
    }
};


template<int Dimension>
class QuadTree : public DataStructure<Dimension>{
    Node<Dimension>* head;
    int size = 0;
public:
    QuadTree(std::vector<Point<Dimension> >& points, int maxPointSize){
        head = new Node<Dimension>(points.begin(), points.end(), maxPointSize);
//        head-> printNodes();
        size = points.size();
    }
    std::vector<Point<Dimension>> rangeQuery(Point<Dimension> center, float radius){
        std::vector<Node<Dimension>*> nodes;
        std::vector<Point<Dimension>> ret;
        head->nearRadius(center, radius, nodes);
        for(Node<Dimension> * n : nodes){
            for(int i = 0; i < n->points.size(); i++){
                float dist = distance(center,  n -> points[i]);
                if(dist <= radius){
                    ret.push_back(n -> points[i]);
                }
            }
        }
        return ret;
    }
    
    std::vector<Point<Dimension>> KNN(Point<Dimension>& center, int k){
        if(k > size){
            throw "k cannot be greater than the tree size";
        }
        std::vector<Node<Dimension>*> nodes;
        std::vector<Point<Dimension>> ret;
        auto pq = std::priority_queue<Point<Dimension>, std::vector<Point<Dimension>>, DistanceComparator<Dimension>> (DistanceComparator<Dimension>(center));
        head-> nearestNeighbors(center, k, std::numeric_limits<float>::max(), pq);
        
        while(!pq.empty()){
            ret.push_back(pq.top());
            pq.pop();
        }
        return ret;
    }
    
    int getSize(){
        return size;
    }
    
    AABB<Dimension> getBoxBounds(){
        return head->aabb;
    }
};
#endif /* QuadTree_hpp */
