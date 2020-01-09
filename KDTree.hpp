#pragma once
#include <stdio.h>
#include <vector>
#include <memory>
#include "Point.hpp"
#include <queue>
#include <algorithm>

template<int Dimension>
class KDTree{
    int size;
    /**
     The Node of our KD Tree
     **/
    template<int SplitDimension> //Don't store the split dimension explicitly
    struct Node{
        Point<Dimension> p;
        AABB<Dimension> aabb;
        std::unique_ptr<Node< (SplitDimension + 1)%Dimension> >  left, right;
        //    template<typename Iter>
        Node(typename std::vector<Point<Dimension>>::iterator begin,typename std::vector<Point<Dimension>>::iterator end, AABB<Dimension> _aabb){
            using ChildType = Node<(SplitDimension +1)%Dimension>;
            aabb = _aabb;
            auto midIter = findMedian(begin, end);
            p = *(midIter);
            int leftSize = midIter - begin;
            int rightSize = end - (++midIter);
            
            if(leftSize > 0){
                auto newEndPt = *(midIter - 1);
                AABB<Dimension> leftAabb = aabb;
                leftAabb.maxs[SplitDimension] = newEndPt[SplitDimension];
                left = std::make_unique<ChildType>(begin, --midIter, leftAabb);
            }else{
                left = nullptr; //leaf node
            }
            
            if(rightSize > 0){
                auto newStartPt = *(midIter + 1);
                 AABB<Dimension> rightAabb = aabb;
                rightAabb.mins[SplitDimension] = newStartPt[SplitDimension];
                
                right = std::make_unique<ChildType>(++midIter, end, rightAabb);
            }else{
                right = nullptr;
            }
        }

        /**
         partions the underlying data structure so that the nth element is in the place it would be if it were sorted.
         returns the size of the middle. Need to subtract one when indexing with an iterator cause begin counts as a position.
         **/
        typename std::vector<Point<Dimension>>::iterator findMedian(typename std::vector<Point<Dimension>>::iterator begin,typename std::vector<Point<Dimension>>::iterator end){
            if(end - begin < 2){
                return begin;
            }
            auto midIter = begin + ((end -begin) /2);
            CompareBy<SplitDimension> compare;
            std::nth_element(begin, midIter, end, compare);
            return midIter;
        }
        
        void RangeQuery(Point<Dimension> &center, float radius, std::vector<Point<Dimension>> &ret){
            float dist = distance(p, center);
            if(dist <= radius){
                ret.push_back(p);
            }
            
            if(left != nullptr){
                Point<Dimension> closestPoint = left->aabb.closestInBox(center);
                float dist = distance(closestPoint, center);
                if(dist <= radius){
                    left->RangeQuery(center, radius, ret);
                }
            }
            
            if(right != nullptr){
                Point<Dimension> closestPoint = right->aabb.closestInBox(center);
                float dist = distance(closestPoint, center);
                if(dist <= radius){
                    right->RangeQuery(center, radius, ret);
                }
            }
        }
        
        void nearestNeighbors(Point<Dimension> center,
                              int k,
                              float radius,
                              std::priority_queue<Point<Dimension>, std::vector<Point<Dimension>>, DistanceComparator<Dimension>> &pq){
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
            
            if(left != nullptr){
                Point<Dimension> closestPoint = left->aabb.closestInBox(center);
                float dist = distance(closestPoint, center);
                if(dist <= radius){
                    left->nearestNeighbors(center, k, radius, pq);
                }
            }
            
            if(right != nullptr){
                Point<Dimension> closestPoint = right->aabb.closestInBox(center);
                float dist = distance(closestPoint, center);
                if(dist <= radius){
                    right->nearestNeighbors(center, k, radius, pq);
                }
            }
        }
        
        void printNodePoints(typename std::vector<Point<Dimension>>::iterator begin,typename std::vector<Point<Dimension>>::iterator end){
            for( auto i = begin; i != end; ++i){
                std::cout << *i << std::endl;
            }
            std::cout << std::endl;
        }
        
        void print(){
            if(left != nullptr){
                left-> print();
            }
            std::cout << p << std::endl;
            if(right != nullptr){
                right-> print();
            }
        }
    };
    
public:
    
    Node<0>* head;
    int getSize(){
        return size;
    }
    KDTree(std::vector<Point<Dimension>>& points){
        size = points.size();
        AABB<Dimension> aabb = getBounds(points);
        head = new Node<0>(points.begin(), points.end(), aabb);
    }
    void print(){
        std::cout << "printing points" << std::endl;
        head->print();
    }

    std::vector<Point<Dimension>> RangeQuery(Point<Dimension> &center, float radius){
        std::vector<Point<Dimension>> ret;
        head->RangeQuery(center, radius, ret);
        return ret;
    }
    
    std::vector<Point<Dimension>> KNN(Point<Dimension> &center, int k){
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
    
    AABB<Dimension> getBoxBounds(){
        return head->aabb;
    }
};

