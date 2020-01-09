//
//  main.cpp
//  test
//
//  Created by Jeremy Rodgers on 6/13/19.
//  Copyright © 2019 Jeremy Rodgers. All rights reserved.
//

#include <iostream>
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "BucketKNN.hpp"
#include "QuadTree.hpp"
#include "KDTree.hpp"

template<int Dimension>
struct resultSort {
public:
    resultSort(){}
    bool operator()(const Point<Dimension>& lhs, const Point<Dimension>& rhs){
//        std::cout << "lhs: " << lhs << "rhs: " << rhs << std::endl;
        if(lhs[0] < rhs[0]){
            return true;
        }else if(lhs[0] == rhs[0]){
            if(lhs[1] < rhs[1]){
                return true;
            }else{
                return false;
            }
        }
        return false;
    }
};



template<int Dimension>
std::vector<Point<Dimension>> createPointVect(int size){
    std::vector<Point<Dimension>> pointsVect;
    for(int i = 0; i < size; i++){
        Point<Dimension> point = Point<Dimension>();
        point.point = std::array<float, Dimension> {(float)i,(float)i, (float) i};
        pointsVect.push_back(point);
    }
    return pointsVect;
}
/**
 creates a point with the save value for all dimensions
 **/
template <int Dimension>
Point<Dimension>generateSinglePoint(int val){
    Point<Dimension> p;
    for(int dim = 0; dim < Dimension; dim++){
        p.point[dim] = val;
    }
    return p;
}


template <int Dimension>
std::vector<Point<Dimension>> generatePoints(std::vector<int> pointVals){
    std::vector<Point<Dimension>> points;
    for(int i = 0; i < pointVals.size(); i++){
        points.push_back(generateSinglePoint<Dimension>(pointVals[i]));
    }
    return points;
}

template <int Dimension>
std::vector<Point<Dimension>> generatePoints(std::vector<std::vector<int>> pointVals){
    std::vector<Point<Dimension>> points;
    for(int i = 0; i < pointVals.size(); i++){
        Point<Dimension> p;
        for(int dim = 0; dim < Dimension; dim++){
            p.point[dim] = pointVals[i][dim];
        }
        points.push_back(p);
    }
    return points;
}

template<int Dimension>
void TestCalcBucketSize(){
    BucketKNN<3> buckets = BucketKNN<3>();
    Point<Dimension> bucketSize;
    AABB<Dimension> aabb; //= buckets.getTestBounds();
    
    for(int i = 0; i < Dimension; i++){
        aabb.mins[i] = 0;
    }
    
    for(int i = 0; i < Dimension; i++){
        aabb.maxs[i] = 10;
    }
    
    int divisions = 3;//= buckets.getDivisions();
    buckets.calcBucketSizes(aabb, divisions, bucketSize);
    std::cout << "Bucket Sizes" << bucketSize<< std::endl;
    std::vector<int> v = {0,1,2,3,4,5,6,7,8,9,10};
    std::vector<Point<3>> points = generatePoints<3>(v);
    for(int i = 0; i < points.size(); i++){
        Point<Dimension> point_coord = buckets.getCoordinates(aabb, bucketSize, divisions, points[i]);
        std::cout << points[i] << "index: " << point_coord.coordinates << std::endl;
    }
    buckets.coordinatesToIndex(points[2].coordinates);
}

template < int Dimension>
void printPoints(std::vector<Point<Dimension>> points){
    for(Point<Dimension> p : points){
        std::cout << p << std::endl;
    }
    std::cout << std::endl;
}

int getShortestResult(std::vector<long> resultSizes){
    int min =  std::numeric_limits<int>::max();
    for(int s: resultSizes){
        if(s < min ){
            min = s;
        }
    }
    return min;
}

void sortResults(std::vector<std::vector<Point<2>>> results){
    for(std::vector<Point<2>> resultSet: results){
        std::sort<Point<2>, resultSort<2>>(resultSet.begin(), resultSet.end(), resultSort<2>());
    }
}

TEST_CASE( "Test BucketKNN Constructor Methods", "[BucketKNN]" ) {
    //Specific Points
    std::vector<int> v = {0,1,2,3,4,5,6,7,8,9,10};
    std::vector<Point<3>> points = generatePoints<3>(v);
    int divisions = 3;//= buckets.getDivisions();
    BucketKNN<3> buckets1 = BucketKNN<3>(points, divisions);
    
    SECTION("Test calc bucket sizes"){
        float expectedBucketSizes = 10/3.0;
        REQUIRE(buckets1.getBucketSizes()[0] == expectedBucketSizes);
        REQUIRE(buckets1.getBucketSizes()[1] == expectedBucketSizes);
        REQUIRE(buckets1.getBucketSizes()[2] == expectedBucketSizes);
    }
    
    SECTION("Test init buckets"){
        REQUIRE(buckets1.getBuckets().size() == pow(divisions, buckets1.getDimensions()));
    }
    
    SECTION("Test bucket coordinates"){
        float expectedBucketSizes = 10/3.0;
        
        for( Point<3> p : points){
            for(int dim = 0; dim < buckets1.getDimensions(); dim++){
                if(p[dim] <= expectedBucketSizes){
                    REQUIRE( p.coordinates[dim] == 0);
                }else if(p[dim] > expectedBucketSizes && p[dim] <= expectedBucketSizes * 2){
                    REQUIRE( p.coordinates[dim] == 1);
                }else if(p[dim] > expectedBucketSizes * 2 && p[dim] <= expectedBucketSizes * 3){
                    REQUIRE( p.coordinates[dim] == 2);
                }else{
                    std::cout << "point would be out of points " << p << std::endl;
                    REQUIRE( p.coordinates[dim] == 2);
                }
            }
        }
    }
    
    SECTION("Test bucket sizes"){
        std::vector<Bucket<3>> buckets = buckets1.getBuckets();
        REQUIRE(buckets.at(0).points.size() == 4);
        REQUIRE(buckets.at(13).points.size() == 3);
        REQUIRE(buckets.at(26).points.size() == 4);
    }
    
}

TEST_CASE("Range Query Tests"){
    std::vector<int> v = {0,1,2,3,4,5,6,7,8,9,10};
    std::vector<Point<3>> points = generatePoints<3>(v);
    int divisions = 3;//= buckets.getDivisions();
    BucketKNN<3> buckets1 = BucketKNN<3>(points, divisions);
    
    SECTION("Test Min and Max"){
        Point<3> p;
        p.point[0] = 5;
        p.point[1] = 5;
        p.point[2] = 5;
        std::vector<Point<3>> result = buckets1.rangeQuery(p, 4);
        //        for( int i = 0; i < result.size(); i++){
        //            std::cout << "result item: " << result[i] << std::endl;
        //        }
        REQUIRE( result.size() == 5 );
    }
}

TEST_CASE("2D Range Query Tests"){
    std::vector<int> v = {0,1,2,3,4,5,6,7,8,9,10};
    std::vector<Point<2>> points = generatePoints<2>(v);
    int divisions = 3;//= buckets.getDivisions();
    BucketKNN<2> buckets1 = BucketKNN<2>(points, divisions);
    
    SECTION("Test Min and Max"){
        Point<2> p;
        p.point[0] = 5;
        p.point[1] = 5;
        std::vector<Point<2>> result = buckets1.rangeQuery(p, 4);
        REQUIRE( result.size() == 5 );
    }
}

TEST_CASE("2D KNN Tests"){
    std::vector<int> v = {0,1,2,3,4,5,6,7,8,9,10};
    std::vector<Point<2>> points = generatePoints<2>(v);
    int divisions = 3;//= buckets.getDivisions();
    BucketKNN<2> buckets1 = BucketKNN<2>(points, divisions);
    
    SECTION("KNN sorting"){
        Point<2> p;
        p.point[0] = 5;
        p.point[1] = 5;
        std::vector<Point<2>> result = buckets1.KNN(p, 4);
        Bucket<2> resultBucket = Bucket<2>(result);
        REQUIRE( result.size() == 4 );
        REQUIRE( resultBucket.contains(generateSinglePoint<2>(v.at(3))) );
        REQUIRE( resultBucket.contains(generateSinglePoint<2>(v.at(4))) );
        REQUIRE( resultBucket.contains(generateSinglePoint<2>(v.at(5))) );
        REQUIRE( resultBucket.contains(generateSinglePoint<2>(v.at(6))) );
    }
}

//TEST_CASE("Oct Trees"){
//    std::vector<std::vector<int>> v = {
//        {1,1,1}, {2,2,2}, {3,3,3},  //back NE
//        {1,-1,1}, {2,-2,2}, {3,-3,3},//back SE
//        {1,1,-1}, {2,2,-2}, {3,3,-3},//front NE
//        {1,-1,-1}, {2,-2,-2}, {3,-3,-3}, //front SE
//        {-1,1,1}, {-2,2,2}, {-3,3,3},//back NW
//        {-1,-1,1}, {-2,-2,2}, {-3,-3,3},//back SW
//        {-1,1,-1}, {-2,2,-2}, {-3,3,-3},//front NE
//        {-1,-1,-1}, {-2,-2,-2}, {-3,-3,-3}, //front SE
//        
//    };
//    std::vector<Point<3>> points = generatePoints<3>(v);
//    std::cout << "start points" << std::endl;
//    for(int i = 0; i < points.size(); i++){
//        std::cout << points[i] << std::endl;
//    }
//    std::cout << "end of points\n" << std::endl;
//    QuadTree<3> qt = QuadTree<3>(points, 3);
//    SECTION("Constructor works"){
//        
//    }
//}

TEST_CASE("Quad Trees Split"){
    std::vector<std::vector<int>> v = {
        {1,1}, {1,-1}, {-1,-1}, {-1,1},
        {2,2}, {2,-2}, {-2,-2}, {-2,2},
        {3,3}, {3,-3}, {-3,-3}, {-3,3},
        {4,4}, {4,-4}, {-4,-4}, {-4,4},
        {5,5}, {5,5}, {-5,-5}, {-5,5},
        {6,6}, {6,-6}, {-6,-6}, {-6,6},
        {7,7}, {7,-7}, {-7,-7}, {-7,7},
    };
    std::vector<Point<2>> points = generatePoints<2>(v);
    //    std::cout << "start points" << std::endl;
    //    for(int i = 0; i < points.size(); i++){
    //        std::cout << points[i] << std::endl;
    //    }
    //    std::cout << "end of points\n" << std::endl;
    QuadTree<2> qt = QuadTree<2>(points, 3);
    SECTION("QT Range Query"){
        Point<2> center = generateSinglePoint<2>(0);
        std::vector<Point<2>> withinRadius = qt.rangeQuery(center, 3.5);
        //        std::cout << "Within radius: \n\n";
        //        printPoints(withinRadius);
        for(Point<2> p : withinRadius){
            REQUIRE(distance(p, center) <= 3);
        }
    }
    
    SECTION("QT Range Query"){
        Point<2> center = generateSinglePoint<2>(0);
        std::vector<Point<2>> nearestNeighbors = qt.KNN(center, 5);
        //        std::cout << "Nearest Neghbors: \n\n";
        //        printPoints(nearestNeighbors);
        REQUIRE(nearestNeighbors.size() == 5);
    }
}

TEST_CASE(){
    std::vector<std::vector<int>> v = {
        {1,1}, {1,-1}, {-1,-1}, {-1,1},
        {2,2}, {2,-2}, {-2,-2}, {-2,2},
        {3,3}, {3,-3}, {-3,-3}, {-3,3},
        {4,4}, {4,-4}, {-4,-4}, {-4,4},
        {5,5}, {5,5}, {-5,-5}, {-5,5},
        {6,6}, {6,-6}, {-6,-6}, {-6,6},
        {7,7}, {7,-7}, {-7,-7}, {-7,7},
    };
    std::vector<Point<2>> points = generatePoints<2>(v);
    KDTree<2> kdTree = KDTree<2>(points);
    std::vector<Point<2>> smlEvenPts = std::vector<Point<2>> {points[0], points[12], points[8], points[4]};
    std::vector<Point<2>> negSmlEvenPts = std::vector<Point<2>> {points[2], points[8], points[0], points[1]};
    std::vector<Point<2>> smlPts = std::vector<Point<2>> {points[0], points[12], points[8]};
    std::vector<Point<2>> negSmlPts = std::vector<Point<2>> {points[2], points[8], points[0]};
    //
    //    SECTION("TEST Evens INIT KD"){
    //
    //        KDTree<2> kdTree = KDTree<2>(smlEvenPts);
    //        KDTree<2> negKdTree = KDTree<2>(negSmlEvenPts);
    //    }
    //    SECTION("TEST Odds INIT KD"){
    //
    //        KDTree<2> kdTree = KDTree<2>(smlPts);
    //        KDTree<2> negKdTree = KDTree<2>(negSmlPts);
    //    }
    //    SECTION("TEST Print KD"){
    //        KDTree<2> kdTree = KDTree<2>(smlEvenPts);
    ////        kdTree.print();
    //        KDTree<2> negKdTree = KDTree<2>(negSmlPts);
    //    }
    SECTION("range query"){
        Point<2> center = generateSinglePoint<2>(0);
        std::vector<Point<2>> ret = kdTree.RangeQuery(center, 3);
        std::cout << std::endl;
        std::cout << "range query result"<< std::endl;
        printPoints(ret);
    }
    
    SECTION("range query"){
        Point<2> center = generateSinglePoint<2>(0);
        std::vector<Point<2>> ret = kdTree.KNN(center, 5);
        std::cout << std::endl;
        std::cout << "KNN query result"<< std::endl;
        printPoints(ret);
    }
    
    
}

TEST_CASE("compare results"){
    std::vector<std::vector<int>> v = {
        {1,1}, {1,-1}, {-1,-1}, {-1,1},
        {2,2}, {2,-2}, {-2,-2}, {-2,2},
        {3,3}, {3,-3}, {-3,-3}, {-3,3},
        {4,4}, {4,-4}, {-4,-4}, {-4,4},
        {5,5}, {5,5}, {-5,-5}, {-5,5},
        {6,6}, {6,-6}, {-6,-6}, {-6,6},
        {7,7}, {7,-7}, {-7,-7}, {-7,7},
    };
    std::vector<Point<2>> points = generatePoints<2>(v);
    KDTree<2> kdTree = KDTree<2>(points);
    QuadTree<2> qt = QuadTree<2>(points, 3);
    BucketKNN<2> buckets = BucketKNN<2>(points, 3);


    SECTION("query check"){
        Point<2> center = generateSinglePoint<2>(0);
        std::vector<Point<2>> KDKnnResult = kdTree.KNN(center, 5);
        std::vector<Point<2>> KDRqResult =  kdTree.RangeQuery(center, 3);

        std::vector<Point<2>> QTKnnResult = qt.KNN(center, 5);
        std::vector<Point<2>> QTRqResult = qt.rangeQuery(center, 3);

        std::vector<Point<2>> BKnnResult = buckets.KNN(center, 5);
        std::vector<Point<2>> BRqResult = buckets.rangeQuery(center, 3);
        std::vector<long> resultSizes;
        resultSizes.push_back(KDKnnResult.size());
        resultSizes.push_back(KDRqResult.size());
        resultSizes.push_back(QTKnnResult.size());
        resultSizes.push_back(QTRqResult.size());
        resultSizes.push_back(BKnnResult.size());
        resultSizes.push_back(BRqResult.size());

        std::vector<std::vector<Point<2>>> results;
        results.push_back(KDKnnResult);
        results.push_back(KDRqResult);
        results.push_back(QTKnnResult);
        results.push_back(QTRqResult);
        results.push_back(BKnnResult);
        results.push_back(BRqResult);
        
        sortResults(results);

        printPoints(KDKnnResult);
        printPoints(KDRqResult);
        printPoints(QTKnnResult);
        printPoints(QTRqResult);
        printPoints(BKnnResult);
        printPoints(BRqResult);

        int shortestLen = getShortestResult(resultSizes);
        std::cout << "comparing first " << shortestLen << "results" << std::endl;
        for(int i = 0; i < shortestLen; i++){
//            REQUIRE(KDKnnResult[i] == QTKnnResult[i]);
//            REQUIRE(KDRqResult[i] == QTRqResult[i]);
//            REQUIRE(BKnnResult[i] == QTKnnResult[i]);
//            REQUIRE(BRqResult[i] == QTRqResult[i]);
        }
    }
}

TEST_CASE("TEST sort"){
    std::vector<std::vector<int>> v = {
                {1,2},{1,1},
                {2,2},{2,1},
                {3,3},{3,3},
                {4,4},
                {5,5},{5,5},
                {6,6},
                {7,7},
            };
   std::vector<Point<2>> points  = generatePoints<2>(v);
    std::vector<std::vector<Point<2>>> results {points};
//    printPoints(results[0]);
    sortResults(results);
}
