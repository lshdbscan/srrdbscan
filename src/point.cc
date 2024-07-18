#include <iostream>
#include <sstream>
#include <point.h>
#include <transformation.h>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <x86intrin.h>


std::vector<uint64_t> HASH;

#define AVX2

BasePoint::BasePoint()
{}

BasePoint::BasePoint(std::istringstream & str_stream, size_t id_)
{
  double val;
  while (str_stream >> val)
    features.push_back(val);

  this->id = id_;
}

void dataset::readData(std::vector<std::vector<float>>& data) 
{
  size_t id = 0;
  for (auto& vec: data) {
    points.push_back(point(vec, id++));
  }
  numberOfDimensions = data[0].size();
}


void dataset::readData(std::string file_name)
{
  std::cout << "Reading input dataset..." << std::endl;
  
  std::string line;
  size_t id = 0;

  std::ifstream f(file_name);
  
  if (!f)
    {
      std::cerr << "Could not open the input file" << std::endl;
      exit(EXIT_FAILURE);
    }
  
  while (std::getline(f, line))
    {
      std::istringstream str_stream(line);
      points.push_back(point(str_stream, id++));
    }

  if (points.size() == 0)
    {
      std::cerr << "Error while reading data: "
		<< "Input dataset is empty" << std::endl;
      exit(EXIT_FAILURE);
    }

  this->numberOfDimensions = points[0].features.size();
  
  for (const auto & point : points)
    if (this->numberOfDimensions != point.features.size())
      {
	std::cerr << "Error while reading data:"
		  << "Inconsistent dimensionality among data points.\n" << std::endl;
	exit(EXIT_FAILURE);
      }
}

std::ostream & dataset::printData(std::ostream & stream, char deli, bool onlyLabel) const
{
  for ( const auto & p : this->points)
    p.print(stream, deli, onlyLabel) << std::endl;

  return stream;
}

void dataset::relabelData()
{
  for (auto & point : points)
    {
      point.unsafe_compress();
      point.reLabel();
    }
}

void dataset::resetData()
{
  for (auto & point : points)
    point.reset();
}

bool HashedPoint::operator == (const HashedPoint & p) const
{
  if (this->features.size() != p.features.size())
    {
      std::cerr << "Can not compare two vectors with different dimensions" << std::endl;
      exit(EXIT_FAILURE);
    }

  for (size_t i = 0; i < this->features.size(); i++)
    {
      if (this->features[i] != p.features[i])
	return false;
    }

  return true;
}

std::ostream & point::print(std::ostream & stream, char deli, bool onlyLabel) const
{
  if (onlyLabel == false)
    BasePoint::print(stream, deli) << deli << (cFindRoot()->id)
				   << deli << std::setw(3) << (corePoint ? "(c) " : " ")
				   << deli << label;
  else
    // stream << this->id; 
    stream << this-> label;
  
  return stream;
}

int point::print(){
  return this-> label;
}

/*
point* point::findRoot()
{
  point* parent = this->parent;
  
  if (parent == NULL)
    {
      return this;
    }
  else
    {
      return (parent->findRoot());
    }
}
*/

point* point::findRoot()
{
  point* u = this;
  while (true)
    {
      point* v = u->parent;
      if ( v == NULL )
	return u;
      point* w = v->parent;
      if ( w != NULL)
	__sync_bool_compare_and_swap (&(u->parent), v, w);
      u = v;
    }
}

void point::unsafe_compress()
{
  if (this->parent != NULL)
    {
      this->parent = this->findRoot();
    }
}

const point* point::cFindRoot() const
{
  const point* p = this;
  const point* father = this->parent;
  
  while (father != NULL)
    {
      p = father;
      father = p->parent;
    }
  
  return p;
}

void point::link(point * p)
{
  point* parent_of_this = this->findRoot();
  point* parent_of_s2   = p->findRoot();

  if (parent_of_this->label == NOISE)
    parent_of_this->label = NON_NOISE;

  if (parent_of_s2->label == NOISE)
    parent_of_s2->label = NON_NOISE;
  
  if (parent_of_this->id == parent_of_s2->id)
  {
      assert(parent_of_this == parent_of_s2); // debug assert
      return;
  }

  if (parent_of_this->weight <= parent_of_s2->weight)
  {
      if (__sync_bool_compare_and_swap (&(parent_of_this->parent), NULL, parent_of_s2)==true)
	    {
        parent_of_s2->weight += parent_of_this->weight;
	      return;
	    }
  }

  if (parent_of_this->weight > parent_of_s2->weight)
  {
      if (__sync_bool_compare_and_swap (&(parent_of_s2->parent), NULL, parent_of_this)==true)
	    {
        parent_of_this->weight += parent_of_s2->weight;
	      return;
	    }
  }

  parent_of_this->link(parent_of_s2);
}

void point::reLabel()
{
  if (label == NON_NOISE)
    {
      label = findRoot()->id;

      assert(label != NOISE);
    }
}

void point::reset()
{
  this->parent = NULL;
  this->label = NOISE;
  this->corePoint = false;
}

bool point::isCore()
{
  return this->corePoint;
}

void point::setAsCorePoint()
{
  this->corePoint = true;
}

/*******************************BasePoint*******************************/
std::ostream & BasePoint::print(std::ostream & stream, char deli) const
{
  stream << std::setw(3) << id << deli;
  for (const auto & feature : features)
    stream << std::setprecision(3) << feature << deli;
  
  return stream;
}

BasePoint BasePoint::operator+(const BasePoint & rhs) const
{
  assert(features.size() == rhs.features.size());

  BasePoint res;
  res.features.reserve(rhs.features.size());
  
  for (size_t i = 0; i < rhs.features.size(); i++)
    {
      res.features[i] = features[i] + rhs.features[i];
    }
  
  return res;
}

double BasePoint::innerProduct(const BasePoint & h) const
{
  assert(features.size() == h.features.size());
  return std::inner_product(features.cbegin(), features.cend(), h.features.begin(), 0.0);
}

inline float hsum256_ps_avx(__m256 v) {
  const __m128 x128 = _mm_add_ps(_mm256_extractf128_ps(v, 1), _mm256_castps256_ps128(v));
  const __m128 x64 = _mm_add_ps(x128, _mm_movehl_ps(x128, x128));
  const __m128 x32 = _mm_add_ss(x64, _mm_shuffle_ps(x64, x64, 0x55));
  return _mm_cvtss_f32(x32);
}

float BasePoint::squaredEuclideanDistance(const BasePoint & p) const
{
#ifdef AVX2
  float result=0;
  auto x = &p.features[0];
  auto y = &features[0];
  auto f = features.size();

  if (f > 7) {
    __m256 d = _mm256_setzero_ps();
    for (; f > 7; f -= 8) {
      const __m256 diff = _mm256_sub_ps(_mm256_loadu_ps(x), _mm256_loadu_ps(y));
      d = _mm256_add_ps(d, _mm256_mul_ps(diff, diff)); // no support for fmadd in AVX...
      x += 8;
      y += 8;
    }
    // Sum all floats in dot register.
    result = hsum256_ps_avx(d);
  }
  // Don't forget the remaining values.
  for (; f > 0; f--) {
    float tmp = *x - *y;
    result += tmp * tmp;
    x++;
    y++;
  }
  return result;

#else
  assert(p.features.size() == features.size());

  double result = 0.0;

  for (size_t i = 0; i < p.features.size(); i++)
    {
      double val = features[i]-p.features[i];
      val = val * val;
      result += val;
    }
  
  return result;
#endif
}

/*

double BasePoint::innerProduct(const BasePoint & h) const
{
  assert(features.size() == h.features.size());

  double result = 0.0;
  __m256d sum = _mm256_setzero_pd();
  
  // double *a1 = &features[0];
  auto a1 = &features[0];
  auto a2 = &h.features[0];
  auto a1End = &features[features.size()];
  
  while(a1 != a1End)
    {
      __m256d v1 = _mm256_loadu_pd(a1);
      __m256d v2 = _mm256_loadu_pd(a2);
      sum = _mm256_add_pd(sum, _mm256_mul_pd(v1, v2));
      a1 += 4;
      a2 += 4;
    }
  double buffer[4];
  _mm256_storeu_pd(buffer, sum);
  return buffer[0] + buffer[1] + buffer[2] + buffer[3];

  // return std::inner_product(features.cbegin(), features.cend(), h.features.begin(), 0.0);
}

double BasePoint::squaredEuclideanDistance(const BasePoint & p) const
{
  assert(p.features.size() == features.size());
  __m256d sum = _mm256_setzero_pd();
  
  auto a1 = &features[0];
  auto a2 = &p.features[0];
  auto a1End = &features[features.size()];
  
  while (a1 != a1End)
    {
      __m256d v = _mm256_sub_pd(_mm256_loadu_pd(a1), _mm256_loadu_pd(a2));
      sum = _mm256_add_pd(sum, _mm256_mul_pd(v, v));
      a1 += 4;
      a2 += 4;
    }
  double buffer[4];
  _mm256_storeu_pd(buffer, sum);
  return buffer[0] + buffer[1] + buffer[2] + buffer[3];
}
*/

double BasePoint::norm() const
{
  if (normalized == true)
    return 1.0;
  
  double result = 0.0;
  for (size_t i = 0; i < features.size(); i++)
    {
      result += features[i]*features[i];
    }
  result = sqrt(result);
  return result;
}

void BasePoint::normalize(bool auxDim)
{
  if (normalized == true)
    {
      std::cerr << "point is already normalized and it should not be normalized again"
		<< std::endl;
      exit(EXIT_FAILURE);
    }

  if (auxDim == true)
    features.push_back(1.0);
  
  const double norm_ = this->norm();

  std::transform(features.begin(),
		 features.end(),
		 features.begin(),
		 [&norm_](double val)->double
		 {
		   return val/norm_;
		 });
  
  normalized = true;
}

BasePoint& BasePoint::operator+=(const BasePoint& rhs)
{
  assert(this->features.size() == rhs.features.size());

  for (size_t i = 0; i < rhs.features.size(); i++)
    {
      this->features[i] += rhs.features[i];
    }
  
  return *this;
}

BasePoint& BasePoint::operator-=(const BasePoint& rhs)
{
  assert(this->features.size() == rhs.features.size());

  for (size_t i = 0; i < rhs.features.size(); i++)
    {
      this->features[i] -= rhs.features[i];
    }
  
  return *this;
}

BasePoint& BasePoint::operator/=(const double& rhs)
{
  if (rhs == 0.0)
    {
      std::cerr << "Error in operator /=: "
		<< "Division by zero." << std::endl;
      exit(EXIT_FAILURE);
    }
  
  for (size_t i = 0; i < this->features.size(); i++)
    {
      this->features[i] /= rhs;
    }
  
  return *this;
}

void dataset::normalizeData()
{
  this->numberOfDimensions++;
  std::cout << "Normalizing Data ..." << std::endl;
  for (auto & point : points)
    point.normalize();
}

// void dataset::meanRemoveData()
// {
//   std::cout << "Mean Removing Data ..." << std::endl;
  
//   point mean(std::vector<float> (numberOfDimensions, 0));
//   for (const auto & p : points)
//     mean += p;

//   mean /= points.size();

//   for (auto & p : points)
//     p -= mean;
// }

void point::setLabel(size_t l)
{
  this->label = l;
}

bool point::isNoise()
{
  return (label==NOISE);
}
