#pragma once
#include "Point.hpp"
#include <vector>
#include <algorithm>
#include "Bucket.hpp"
#include <math.h>
#include <queue>

template<int Dimension>
struct distanceCompare{
    distanceCompare(const Point<Dimension> _p) : center(_p) {}
    bool operator()(const Point<Dimension>& point1, const Point<Dimension>& point2){
        return distanceFrom(point1, center) < distanceFrom(point2, center);
    }
private:
    const Point<Dimension> center;
    float distanceFrom(Point<Dimension> currentPoint, Point<Dimension> centerPoint ){
        float squaredSum = 0;
        for( int dim = 0; dim < Dimension; dim++){
            squaredSum += pow(centerPoint[dim] - currentPoint[dim], 2);
        }
        return sqrt(squaredSum);
    }
};

template<int Dimension> //The dimension of the points.  This can be used like any other constant within.
class BucketKNN : public DataStructure<Dimension>{
    int size;
    int divisions;
    std::vector<Bucket<Dimension>> buckets;
    AABB<Dimension> bounds;
    Point<Dimension> bucketSizes;
public:
    // TEST HELPER FUNCTION
    BucketKNN(){}
    
    BucketKNN(std::vector<Point<Dimension> >& points, int divisions_)
    {
        if(divisions_ <= 0){
            throw "dimesions cannot be negative";
        }
        
        if(points.size() <= 0){
            std::cout << "size must be greater than 0";
        }
        
        size = points.size();
        divisions = divisions_;
        bounds = getBounds(points);
        
        calcBucketSizes();
        long totalBuckets = pow(divisions, Dimension);
        createNewBuckets(totalBuckets);
        
        for(Point<Dimension> &p: points){
            getCoordinates(p);
            long index = coordinatesToIndex(p.coordinates);
            buckets[index].points.push_back(p);
        }
    }
    
    /***
     computes the size of bucket in each Dimension
     ***/
    void calcBucketSizes(){
        if(divisions <= 0){
            throw "divisions can't be less than 1";
        }
        for( int dim = 0; dim < Dimension; dim++ ){
            float totalLen = bounds.maxs[dim] - bounds.mins[dim];
            float size = totalLen / divisions;
            bucketSizes.point[dim] = size;
        }
    }
   
    void createNewBuckets(int totalBuckets){
        buckets = std::vector<Bucket<Dimension>>{};
        for(int i = 0; i < totalBuckets; i++){
            buckets.push_back(Bucket<Dimension>(bucketSizes));
        }
    }
    
    /**
     Gets the coordinates of a bucket in multidimensional space (not a flattened array);
     Param: Point<Dimension> p - the actual, unaltered coordnates of a point
     return: the adapted point with coordinates
     **/
    Point<Dimension> getCoordinates(Point<Dimension> &p ){
        for(int i = 0; i < Dimension; i++){
            // i = point - minimum in the dimension / numDivisions; then clamp betweem 0 index and total size of flattened array.
            float distance =(p[i] - bounds.mins[i]);
            p.coordinates.coordinates[i] = std::clamp( (int) ( distance / bucketSizes[i]) , 0, divisions -1);
        }
        return p;
    }
    
    /**
     Takes coodinates and returns the 1 dimensional array index of a point.
     Param: Coordinates<Dimension> coords - multi dimensional bucket coordinates.
     return: 1 dimensional array index;
     **/
    
    long coordinatesToIndex(Coordinates<Dimension> &coords){
        int index = 0;
        for (int dim = 0; dim < Dimension; dim++){
            long jumpFactor = (Dimension -1) - dim;
            long pw = pow(divisions, jumpFactor);
            index += coords[dim] * pw;
        }
        return index;
    }
    /**
     Given a point and a radius, find all the points that fall within radius.
     Params:Point<Dimension>& p - a centeral point
            float radius - how big of a circle to draw.
     return vector<Points<Dimension> that fall within the range.
     **/
    std::vector<Point<Dimension>> rangeQuery(const Point<Dimension>& centerPoint, float radius){
        for (int dim = 0; dim < Dimension; dim++){
            if(centerPoint[dim] > bounds.maxs[dim] || centerPoint[dim] < bounds.mins[dim]){
                std::cout << "WARNING: center point must be within bounds" << std::endl;
                std::cout << "bounds" << bounds <<std::endl;
                std::cout << "center" << centerPoint <<std::endl;
                throw "out of bounds on range query";
                return {};
            }
        }
        std::vector<Point<Dimension>> ret;
        Point<Dimension> point = centerPoint;
        Coordinates<Dimension> minCoords = findMinBucketCoords(point, radius);
        Coordinates<Dimension> maxCoords = findMaxBucketCoords(point, radius);

        for( Coordinates<Dimension> currentCoords = minCoords;
            currentCoords != nextBucket(maxCoords, minCoords, maxCoords);
            currentCoords = nextBucket(currentCoords, minCoords, maxCoords)){
            
            long index = coordinatesToIndex(currentCoords);
           
            Bucket<Dimension> b =  buckets[index];
            for(int i = 0; i < b.points.size(); i++){
                float distance = distanceFrom(b.points[i], centerPoint);
                if(distance <= radius){
                    ret.push_back(b.points[i]);
                }
            }
        }
        return ret;
    }
    
    Coordinates<Dimension> findMinBucketCoords(Point<Dimension> p, float radius){
        p = p - radius; //overwritten operator to subtrack from point in each dimension
        getCoordinates(p);
        return p.coordinates;
    }
    
    Coordinates<Dimension> findMaxBucketCoords(Point<Dimension> p, float radius){
        p = p + radius;//overwritten operator to add to point in each dimension
        getCoordinates(p);
        return p.coordinates;
    }
    
    float distanceFrom(Point<Dimension> currentPoint, Point<Dimension> centerPoint ){
        float squaredSum = 0;
        for( int dim = 0; dim < Dimension; dim++){
            squaredSum += pow(centerPoint[dim] - currentPoint[dim], 2);
        }
        return sqrt(squaredSum);
    }
    
    Coordinates<Dimension> nextBucket(Coordinates<Dimension> currentBucketCoord, Coordinates<Dimension> minCoords, Coordinates<Dimension> maxCoords){
        currentBucketCoord.coordinates[Dimension -1]++; //increment the last (z) dimension
        
        for( int dim = Dimension -1; dim > 0; dim--){ //loop through the dimensions from last to first (z to x)
            if(currentBucketCoord.coordinates[dim] > maxCoords[dim]){ // if the point is past the max in the current dimension...
                currentBucketCoord.coordinates[dim] = minCoords[dim]; // reset the coordinate of the current dimension to the min
                currentBucketCoord.coordinates[dim - 1]++; // increment the next dimension
            }else{
                break;
            }
        }
        return currentBucketCoord;
        
    }
    
    std::vector<Point<Dimension>> KNN(Point<Dimension>& center, int k){
        if(size < k){
            std::cout << "k cannot be larger than the number of points" << std::endl;
            return {};
        }
        assert(size >= k);
        long radius = getMaxBucketSize(); //start with radius as bucket width in biggest dimension;
        bool searching = true;
        std::vector<Point<Dimension>> ret;
        while( searching ){
            ret = rangeQuery(center, radius);
            if( ret.size() == k ){
                return ret;
            }else if( ret.size() < k ){
                radius *= radius;
            }else{
                searching = false;
            }
        }
       return refineNeighbors(ret, center, k);
    }
    
    int getMaxBucketSize() const{
        int max = 0;
        for(int i = 0; i < Dimension; i++){
            if( bucketSizes[i] > max){
                max = bucketSizes[i];
            }
        }
        return max;
    }
    
    /**
     We have narrowed down to a "community" of points but we still don't have a neighborhood. Widdle down the result of the range query.
     **/
    std::vector<Point<Dimension>> refineNeighbors(std::vector<Point<Dimension>>& community,const Point<Dimension>& centerPoint,  int k){
        std::sort(community.begin(), community.end(), distanceCompare<Dimension>(centerPoint));
        if(community.size() < k ){
            throw "BucketKNN::refineNeighbors: community size must be greater than or equal to k";
        }
        int excessPoints = community.size() - k;
        for(int i = 0; i < excessPoints; i++){
            community.pop_back();
        }
        return community;
    }
    
    std::vector<Bucket<Dimension>> getBuckets(){
        return buckets;
    }
    
    Point<Dimension> getBucketSizes(){
        return bucketSizes;
    }
    
    int getDivisions(){
        return divisions;
    }
    
    AABB<Dimension> getBoxBounds(){
        return bounds;
    }
    
    int getDimensions(){
        return Dimension;
    }
    
    int getSize(){
        return size;
    }
};

