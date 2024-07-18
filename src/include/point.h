#ifndef __POINT_H__
#define __POINT_H__

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
// #include <highfive/H5File.hpp>
// #include <highfive/H5DataSet.hpp>
// #include <highfive/H5DataSpace.hpp>

#define NOISE -1
#define NON_NOISE -2

#define UH_PRIME_DEFAULT 4294967291U

extern std::vector<uint64_t> HASH;


class BasePoint
{
  bool normalized = false;
 public:
  size_t id = -1;
  size_t weight = 1;
  std::vector<float> features;

  BasePoint();
  BasePoint(std::istringstream &, size_t id_ = -1);
  BasePoint(std::vector<float>& vec, size_t id_ = -1)
    {
      std::copy(vec.begin(), vec.end(), std::back_inserter(features));
      this->id = id_;
    }

  void normalize(bool auxDim = true);

  BasePoint& operator+=(const BasePoint& rhs);
  BasePoint& operator-=(const BasePoint& rhs);
  BasePoint& operator/=(const double& rhs);
  
  double norm() const;
  double innerProduct(const BasePoint&) const;
  BasePoint operator+(const BasePoint & rhs) const;
  float squaredEuclideanDistance(const BasePoint &) const;
  std::ostream & print(std::ostream & stream, char deli) const;
};

class Hyperplane : public BasePoint
{
 public:
  using BasePoint::BasePoint;
  Hyperplane(){}
};

class point : public BasePoint
{
 private:
  point* parent = NULL;
  int label = NOISE;
  bool corePoint = false;
    
 public:
  using BasePoint::BasePoint;

  bool processed = false; /*used in vanilla DBSCAN*/
  void setLabel(size_t);  /*used in vanilla DBSCAN*/
  bool isNoise();  /*used in vanilla DBSCAN*/
  point* findRoot();
  const point* cFindRoot() const;
  void unsafe_compress();
  void link(point*);
  void reLabel();
  void reset();
  bool isCore();
  void setAsCorePoint();
  std::ostream & print(std::ostream & stream, char deli, bool onlyLabel = false) const;
  int print();
};

class HashedPoint : public BasePoint
{
 public:
  std::vector<int> features;
  bool operator == (const HashedPoint & p) const;

  uint32_t combine() const {
    uint64_t h = 0;

    for (int i = 0; i < features.size(); i++) {
      h = h + (uint64_t)features[i] * HASH[i];
      h = (h & 4294967295U) + 5 * (h >> 32);
      if (h >= UH_PRIME_DEFAULT) {
        h = h - UH_PRIME_DEFAULT;
      }
    }

    return h;
  }
};
/*
namespace std
{
  template<> struct hash<HashedPoint>
  {
    std::size_t operator() (HashedPoint const & p) const noexcept
    {
      std::ostringstream concat;
      concat.precision(1);

      for (const auto & elem : p.features)
	concat << elem << ",";
      //std::cout << concat.str() << std::endl;
      
      std::size_t h = std::hash<std::string>() (concat.str());

      return h;
    }
  };
}
*/
class dataset
{
 public:
  size_t numberOfDimensions = 0;  
  std::vector<point> points;

  std::string name;
  
  void readData(std::string);
  void readData(std::vector<std::vector<float>>& data);
  void relabelData();
  void resetData();
  void normalizeData();
  void meanRemoveData();
  std::ostream & printData(std::ostream & stream, char deli, bool onlyLabel = false) const;
};


#endif
