#include <randomGen.h>

RandGenerator::RandGenerator():uniform_dist(1, 536870912U)
{
  dist = new std::normal_distribution<double>(0.0, 1.0);
}

double RandGenerator::getRandVal()
{
  return (*dist)(engine);
}

uint32_t RandGenerator::getHashCoeff() {
  return uniform_dist(engine);
}
