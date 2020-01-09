#pragma once

#include <array>
#include <ostream>
#include <cmath>


template<int Dimension>
struct Coordinates{
    //Fixed size, inline array
    std::array<int, Dimension> coordinates;
    int dimensions = Dimension;
    //index points directly.  This is read only.
    float operator[](int index) const { return coordinates[index]; }
    Coordinates<Dimension>& operator=(const Coordinates<Dimension> &rhs){
        for(int dim = 0; dim < Dimension; dim++){
            coordinates[dim] = rhs[dim];
        }
        return *this;
    }
    
    bool operator ==(const Coordinates<Dimension> &rhs){
        for(int dim = 0; dim < Dimension; dim++){
            if(coordinates[dim] != rhs[dim]){
                return false;
            }
        }
        return true;
    }
    
    bool operator !=(const Coordinates<Dimension> &rhs){
        return !(*this == rhs);
    }
};

//Points are Read only
//There's no constructor but you can make one from a std::array of appropriate size
template<int Dimension>
struct Point{
    //Fixed size, inline array
    std::array<float, Dimension> point;
    Coordinates<Dimension> coordinates;
    int dimensions = Dimension;
    //index points directly.  This is read only.
    float operator[](int index) const { return point[index]; }
    
    Point<Dimension>& operator-(float num) {
        for(int dim = 0; dim < Dimension; dim++){
            point.at(dim) = point.at(dim) - num;
        }
        return *this;
    }
    
    Point<Dimension>& operator+(float num) {
        for(int dim = 0; dim < Dimension; dim++){
            point[dim] = point[dim] + num;
        }
        return *this;
    }
    
    Point<Dimension>& operator=(const Point<Dimension> &rhs){
        for(int dim = 0; dim < Dimension; dim++){
            point[dim] = rhs[dim];
        }
        coordinates = rhs.coordinates;
        return *this;
    }
    
    bool operator==(const Point<Dimension> &rhs)const{
        for(int dim = 0; dim < Dimension; dim++){
            if( point[dim] != rhs[dim]){
                return false;
            };
        }
        return true;
    }
    
};


//returns the distance between 2 points.
//Note: we actually don't care about the distance.  We could use the squared distance
//If you want to optimize your datastructures, add a squared distance function and modify your data structures to use it
template<int Dimension>
float distance(const Point<Dimension>& a, const Point<Dimension>& b){
  float dist = 0;
  for(int i = 0; i < Dimension; ++i){
	dist += (a[i] - b[i])*(a[i] - b[i]);
  }
  return std::sqrt(dist);
}

//print a point.  Useful for debugging
template<int Dimension>
std::ostream& operator<<(std::ostream& outs, const Point<Dimension>& p){
  outs << "Point(";
  for(int i = 0; i < Dimension -1; ++i){
	outs << p[i] << ", ";
  }
  outs << p[Dimension -1] << ")";
  return outs;
}

template<int Dimension>
std::ostream& operator<<(std::ostream& outs, const Coordinates<Dimension>& c){
    outs << "Coordinate(";
    for(int i = 0; i < Dimension -1; ++i){
        outs << c[i] << ", ";
    }
    outs << c[Dimension -1] << ")";
    return outs;
}


//This is a comparator that you can pass to std::sort, or nth_element, etc

//Use like CompareBy<0> compareByX;  CompareBy<1> compareByY; std::sort(points.begin(), points.end(), compareByY);
//The template paramater must be known at compile time!
//If you choose to implement your KD tree storing the splitting dimension explicitly, you'll need to
//rewrite this to take the dimension as a constructor parameter, rather than as a template parameter

template<int Dimension>
struct CompareBy{

  template<int PD>
  bool operator()(const Point<PD>& lhs, const Point<PD>& rhs){
	static_assert(Dimension < PD, "must sort by a dimension that exists!");
	//lexographic ordering, after the provided dimension
	if(lhs[Dimension] == rhs[Dimension]){ //equal in the dim we care about, compare other dimensions
	  for(int i = 0; i < PD; ++i){
		if(i != Dimension){
		  if(lhs[i] != rhs[i]){
			return lhs[i] < rhs[i];
		  }
		}
	  }
	  return false; //points are equal
	  
	} else { //easy case
	  return lhs[Dimension] < rhs[Dimension];
	}
  }
  
};


//Same as above, but compares based on distance to a point
//You shouldn't have to provide the template parameter (it should be inferred when you pass a point of a given dimension)

template<int Dimension>
struct DistanceComparator {
public:
  DistanceComparator(const Point<Dimension>& q_ )
	:p{q_} {}

  bool operator()(const Point<Dimension>& lhs, const Point<Dimension>& rhs){
	return distance(p, lhs) < distance(p, rhs);
  }

private:
  Point<Dimension> p;
};



//Axis Aligned Bounding Box
template<int Dimension>
struct AABB{
  std::array<float, Dimension> mins, maxs;
  //default contstruct to be infinite, which is useful for KDTree
  AABB()
  {
	for(int i = 0; i < Dimension; ++i){
	  mins[i] = -std::numeric_limits<float>::max();
	  maxs[i] = std::numeric_limits<float>::max();
	}
  }

  //returns the point in the bounding box (likely on the boundary) that is closest to p
  Point<Dimension> closestInBox(const Point<Dimension>& p ){

	std::array<float, Dimension> arr;
	for(int i = 0; i < Dimension; ++i){
	  arr[i] = std::clamp(p[i], mins[i], maxs[i]);
	}
	return Point<Dimension>{arr};
	
  }
};

/***
 computes the size of bucket in each Dimension
 ***/
template<int Dimension>
Point<Dimension> splitBox(AABB<Dimension> &aabb, int divisions){
    if(divisions <= 0){
        throw "divisions can't be less than 1";
    }
    Point<Dimension> ret;
    for( int dim = 0; dim < Dimension; dim++ ){
        float totalLen = aabb.maxs[dim] - aabb.mins[dim];
        float size = totalLen / divisions;
        ret.point[dim] = aabb.mins[dim] + size;
    }
    return ret;
}

//Return the AABB surrounding a set of points
template<int Dimension>
AABB<Dimension> getBounds(const std::vector<Point<Dimension> > & points){

  std::array<float, Dimension> mins, maxs;
  for(int i = 0; i < Dimension; ++i){
	mins[i] = std::numeric_limits<float>::max();
	maxs[i] = -std::numeric_limits<float>::max();
  }
  for(const auto& p : points){
	for(int i = 0; i < Dimension; ++i){
	  mins[i] = std::min(mins[i], p[i]);
	  maxs[i] = std::max(maxs[i], p[i]);
	}
  }

  AABB<Dimension> ret;
  ret.mins = mins;
  ret.maxs = maxs;
  return ret;
  
}
template<int Dimension>
AABB<Dimension> getBounds(typename std::vector<Point<Dimension>>::iterator begin, typename std::vector<Point<Dimension>>::iterator end){
    typename std::vector<Point<Dimension>>::iterator j;
    std::array<float, Dimension> mins, maxs;
    for(int i = 0; i < Dimension; ++i){
        mins[i] = std::numeric_limits<float>::max();
        maxs[i] = -std::numeric_limits<float>::max();
    }
    for( j = begin; j != end; ++j){
        for(int i = 0; i < Dimension; ++i){
            Point<Dimension> p = *j;
            mins[i] = std::min(mins[i], p[i]);
            maxs[i] = std::max(maxs[i], p[i]);
        }
    }
    
    AABB<Dimension> ret;
    ret.mins = mins;
    ret.maxs = maxs;
    return ret;
}
//Print out an AABB
template<int Dimension>
std::ostream& operator<<(std::ostream& outs, const AABB<Dimension>& aabb){
  outs << "AABB( " << Point<Dimension>{aabb.mins} << " : " << Point<Dimension>{aabb.maxs} << " ) ";
  return outs;
}


template <int Dimension>
class DataStructure{
public:
    virtual std::vector<Point<Dimension>> KNN(Point<Dimension>& center, int k) = 0;
    virtual int getSize() = 0;
    virtual AABB<Dimension> getBoxBounds()= 0;
};
