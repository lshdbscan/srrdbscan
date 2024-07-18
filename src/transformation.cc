#include <algorithm>
#include <iomanip>
#include <cassert>
#include <transformation.h>
#include <globals.h>
#include <PQ.h>
#include <cmath>
#include <unordered_set>
#include <unordered_map>


HashTable::HashTable(dataset* ds_,
		     size_t numberOfHyperplanes_,
		     RandGenerator* gen_)
{
  this->ds  = ds_;
  this->gen = gen_;
  numberOfHyperplanes = numberOfHyperplanes_;
  this->initializeHashTable(numberOfHyperplanes_);
}

HashTable::HashTable(dataset* ds_,
		     std::string& fileName)
{
  this->ds = ds_;
  initializeHashTable(fileName);
}

void HashTable::populateHashTable()
{
  assert(ds != NULL);
  
  HashedPoint hashedPoint;
  hashedPoint.features.reserve(hyperplanes.size());
  for (auto & point : ds->points)
    {

      std::transform(hyperplanes.cbegin(), hyperplanes.cend(),
		     std::back_inserter(hashedPoint.features),
		     [&point](const Hyperplane & h)->double
		     {
		       return hashFunc(point, h);
		     });
    //   myMap[hashedPoint].push_back(&point);
	  auto hash = hashedPoint.combine();
	  //std::cout << hash << std::endl;
	  //hashes[&point] = hash;
	  hashedPoint.features.clear();
	  hashTable[hash].push_back(&point);
    }

	//std::cout << "hash table has " << hashTable.size() << " buckets." << std::endl;

	// shrink-to-fit all buckets

	// for (auto& bucket: hashTable) {
	// 	bucket.second.shrink_to_fit();
	// }
}

void HashTable::populateHashTable(std::vector<point>::iterator begin,
				  std::vector<point>::iterator end)
{
  assert((end - begin) > 0);

  HashedPoint hashedPoint;
  hashedPoint.features.reserve(hyperplanes.size());

  for (std::vector<point>::iterator pointIter = begin;
       pointIter < end;
       pointIter++)
    {

      std::transform(hyperplanes.cbegin(), hyperplanes.cend(),
		     std::back_inserter(hashedPoint.features),
		     [pointIter](const Hyperplane & h)->int
		     {
		       return hashFunc(*pointIter, h);
		     });
    	auto hash = hashedPoint.combine();
	 	hashedPoint.features.clear();
    //   myMap[hashedPoint].push_back(&(*pointIter));
        hashTable[hash].push_back(&(*pointIter));
    }
}
//This is added by outsiders (not from the original IPLSH papar)
void HashTable::populateHashTable(tbb::concurrent_vector<point*>::iterator begin,
				  tbb::concurrent_vector<point*>::iterator end)
{
  assert((end - begin) > 0);

  HashedPoint hashedPoint; hashedPoint.features.reserve(numberOfHyperplanes);
  for (tbb::concurrent_vector<point*>::iterator pointIter = begin;
       pointIter < end;
       pointIter++)
    {

      std::transform(hyperplanes.cbegin(), hyperplanes.cend(),
		     std::back_inserter(hashedPoint.features),
		     [pointIter](const Hyperplane & h)->int
		     {
		       return hashFunc(**pointIter, h);
		     });

	    auto hash = hashedPoint.combine();
	  	//(*pointIter)->save_hash((uint64_t) this, hash);
	 	hashedPoint.features.clear();
      //myMap[hashedPoint].push_back(*pointIter);
        hashTable[hash].push_back(*pointIter);
    }
}

void HashTable::identifyCoreBuckets()
{
  for (auto & element : myMap)
    {
      if (element.second.size() >= minPts)
	{
	  std::vector<point*> coreBucketElements;
	  coreBucketElements.reserve(element.second.size());

	  std::copy(element.second.begin(), element.second.end(),
		    std::back_inserter(coreBucketElements));

	  point* core = coreBucketElements[getMedian(coreBucketElements)];
	  size_t cnt = 0;
	  
	  std::for_each(coreBucketElements.begin(),
			coreBucketElements.end(),
			[core, this, &cnt](point* p)
			{
			  if ( (p->findRoot() == p) &&
			       (distanceFunc(*core, *p)  <= epsilon)
			       )
			    {
			      p->link(core);
			      cnt++;
			    }
			}
			);

	  if (cnt >= minPts)
	    {
	      core->setAsCorePoint();
	      
	      CoreBucket coreBucket;
	      coreBucket.representative = core;
	      coreBucket.members = coreBucketElements;
	      coreBuckets.push_back(coreBucket);
	    }
	}
    }
}

size_t HashTable::getMedian(std::vector<point*> & vec) const
{
  assert(vec.size() != 0);
  std::vector<double> projVals;
  projVals.reserve(vec.size());

  Hyperplane ones;
  ones.features.resize(vec[0]->features.size(), 1.0);

  std::transform(vec.begin(), vec.end(),
		 std::back_inserter(projVals),
		 [&ones](const point* p)->double
		 {
		   return p->innerProduct(ones);
		 }
		 );

  std::vector<size_t> indices;
  
  for (size_t i = 0; i < vec.size(); i++)
    {
      indices.push_back(i);
    }
  
  std::sort(indices.begin(), indices.end(),
	    [&projVals](size_t leftIndex, size_t rightIndex)
	    {
	      return projVals[leftIndex] < projVals[rightIndex];
	    }
	    );

  return indices[vec.size()/2];
}

void HashTable::identifyMergeTasks()
{
  for (const auto & coreBucket : coreBuckets)
    {
      std::vector<point *> toGetMerged;

      for (const auto element : coreBucket.members)
	{
	  if (element->isCore())
	    {
	      toGetMerged.push_back(element);
	    }
	}

      assert(toGetMerged.size() >= 1);
      
      point* corePoint_0 = toGetMerged.back();
      toGetMerged.pop_back();
      
      while (toGetMerged.size() > 0)
	{
	  point* corePoint = toGetMerged.back();
	  toGetMerged.pop_back();

	  mergeTasks.push_back({corePoint_0, corePoint});
	}
    }
}

void HashTable::identifyMergeTasks_V2()
{
  for (const auto & coreBucket : coreBuckets)
    {
      std::vector<point *> toGetMerged;

      for (const auto element : coreBucket.members)
	{
	  if (element->isCore())
	    {
	      toGetMerged.push_back(element);
	    }
	}

      assert(toGetMerged.size() >= 1);

      size_t i, j;

      for ( i = 0; i < toGetMerged.size(); i++)
	{
	  for ( j = i+1; j < toGetMerged.size(); j++)
	    {
	      if (distanceFunc(*toGetMerged[i], *toGetMerged[j]) <= epsilon)
		{
		  mergeTasks.push_back({toGetMerged[i], toGetMerged[j]});
		}
	    }
	}
    }
}

// void HashTable::mergeCorePoints(){
// 	for(auto &bucket: hashTable){ //
// 		size_t n = bucket.second.size();

// 		if (n == 0) {
// 			continue;
// 		}
// 		long long m = 0;
// 		long long necessary = 0;
// 		std::cout << "Have to carry out merging task in bucket of size " << n << std::endl;
		
// 		std::unordered_set<point*> repr;
// 		std::unordered_map<point*, std::vector<point*>> map;


// 		// split up points in bucket by representatitve in UF data structure
// 		for (size_t i = 0; i < n; i++) {
// 			auto p = bucket.second[i]->findRoot();
// 			repr.insert(p);
// 			if (map.find(p) != map.end()) {
// 				map[p].push_back(bucket.second[i]);
// 			} else{
// 				map[p] = std::vector<point*> {bucket.second[i]};
// 			}
// 		}


// 		if (map.size() == 1) {
// 		 	std::cout << "no work needs to be done in this bucket, all points are in the same cluster." << std::endl;
// 		 	continue;
// 		}

// 		std::vector<point*> work;
// 		std::move(repr.begin(), repr.end(), std::back_inserter(work));

// 		// sort lists of equivalence classes to check by size, smaller first.
// 		std::sort(work.begin(), work.end(), [&map] (point* a, point* b) {
// 			return map[a].size() < map[b].size();
// 		});



// 		for (size_t i = 0; i < work.size(); i++) {
// 			std::cout << i << std::endl;
// 			int skipped = 0;
// 			for (size_t j = i + 1; j < work.size(); j++) {
// 				//std::cout << j << std::endl;
// 				// all-to-all between EC i and EC j
// 				auto list1 =  map[work[i]];
// 				auto list2 = map[work[j]];
// 				if (list1[0]->findRoot() == list2[0]->findRoot()) {
// 					skipped += 1;
// 					continue;
// 				}
// 				for (size_t k = 0; k < list1.size(); k++) {
// 					bool linked = false;
// 					for (size_t ell = 0; ell < list2.size(); ell++) {
// 						m += 1;
// 						auto p = list1[k];
// 						auto q = list2[ell];
// 						//std::cout << p->squaredEuclideanDistance(*q) << " " << epsilon << std::endl;
// 						if (p->squaredEuclideanDistance(*q) <= epsilon) {
// 							necessary +=1;
// 							p->link(q);
// 							linked = true;
// 							break;
// 						}
// 					}
// 					if (linked or (list1[0]->findRoot() == list2[0]->findRoot())) {
// 						break;
// 					}
// 				}
// 			}

// 			std::cout << skipped / (double) (work.size() - i) << std::endl;

// 			// check if more work needs to be done
// 			auto p = map[work[i]][0]->findRoot();
// 			bool done = true;
// 			std::unordered_set<point*> visited;
// 			for (size_t k = 0; k < work.size(); k++) {
// 				//std::cout << k << std::endl;
// 				visited.insert(map[work[k]][0]->findRoot());
// 				if (p != map[work[k]][0]->findRoot()) {
// 					//std::cout << i << "-th iteration done, still multiple ECs" << std::endl;
// 					done = false;
// 					//break;
// 				}
// 			}
// 			std::cout << visited.size() << " ECs left" << std::endl; 
// 			if (done) {
// 				break;
// 			}
// 		}
// 		// for(size_t i = 0; i < n; i++){
// 		// 	for(size_t j = i+1; j < n; j++){
// 		// 		if(bucket.second[i]->squaredEuclideanDistance(*bucket.second[j]) <= epsilon){
// 		// 			m += 1;
// 		// 			if (bucket.second[i]->findRoot() != bucket.second[j]->findRoot()) {
// 		// 				necessary += 1;
// 		// 			}
// 		// 			bucket.second[i]->link(bucket.second[j]);
// 		// 		}
// 		// 	}
// 		// }
// 		std::cout << "statistics for bucket with " << map.size() << " keys:" ;
// 		for (const auto& [key, value]: map) {
// 			std::cout << value.size() << " ";
// 		}
// 		std::cout << std::endl;
// 		int frac = double(m) / (n * (n + 1) / 2)  * 100;
// 		std::cout << frac <<"% merges (" << necessary << " necessary merges carried out), " <<  
// 		m << " merges have been carried out (" << n * (n + 1) / 2 << " pairs checked)" << std::endl;

// 	}
// }

void HashTable::mergeCorePoints(){
	for(auto &bucket: hashTable){ //
		size_t n = bucket.second.size();

		if (n == 0) {
			continue;
		}
		long long m = 0;
		long long necessary = 0;
#ifdef SRR_DEBUG
		std::cout << "Have to carry out merging task in bucket of size " << n << std::endl;
#endif
		
		std::unordered_set<point*> repr;
		std::unordered_map<point*, std::vector<std::pair<size_t, point*>>> map;


		// split up points in bucket by representatitve in UF data structure
		for (size_t i = 0; i < n; i++) {
			auto p = bucket.second[i]->findRoot();
			repr.insert(p);
			if (map.find(p) != map.end()) {
				map[p].push_back(std::make_pair(i, bucket.second[i]));
			} else{
				map[p] = std::vector<std::pair<size_t, point*>> {std::make_pair(i, bucket.second[i])};
			}
		}


		if (map.size() == 1) {
#ifdef SRR_DEBUG
		 	std::cout << "no work needs to be done in this bucket, all points are in the same cluster." << std::endl;
#endif
		 	continue;
		}

		auto pq = better_priority_queue::updatable_priority_queue<size_t, int64_t>();

		for (size_t i = 0; i < n; i++) {
			auto p = bucket.second[i];
			auto u = p->findRoot();
			pq.push(i, -map[u].size());
		}

		size_t different_clusters = map.size();
		int rounds = 0;

		while (pq.size() > 0 && different_clusters > 0) {
			rounds += 1;
#ifdef SRR_DEBUG
			std::cout << "Currently in round " << rounds << std::endl;
			//std::cout << "Working on point " << pq.top().key << " with prio " << pq.top().priority << std::endl;
			std::cout << "There are " << map.size() << " clusters left" << std::endl;
#endif
			auto p = bucket.second[pq.pop_value().key];
			auto u = p->findRoot();

			std::vector<point*> remove_list;

			int ecs_inspected = 0;

			for (const auto & [key, value]: map) {
				if (key != u) {
					for (const auto & point: map[key]) {
						auto q = point.second;

						if (p->squaredEuclideanDistance(*q) <= epsilon) {
							p->link(q);
							remove_list.push_back(key);
							break;
						}
					}
					//std::cout << ++ecs_inspected << std::endl;
				}
			}

			//std::cout << "carrying out merge jobs" << std::endl;

			if (remove_list.size() > 0) {
				for (const auto& key: remove_list) {
					// merge ECs
					different_clusters--;
					std::move(map[key].begin(), map[key].end(), std::back_inserter(map[u]));
					map.erase(key);
				}
				for (const auto& q: map[u]) {
					if (q.second != p) {
						pq.update(q.first, -map[u].size());
					}
				}

			}
			//std::cout << different_clusters << std::endl;


	
		}

#ifdef SRR_DEBUG
		std::cout << "statistics for bucket with " << map.size() << " keys:" ;
		for (const auto& [key, value]: map) {
			std::cout << value.size() << " ";
		}
		std::cout << std::endl;
		int frac = double(m) / (n * (n + 1) / 2)  * 100;
		std::cout << frac <<"% merges (" << necessary << " necessary merges carried out), " <<  
		m << " merges have been carried out (" << n * (n + 1) / 2 << " pairs checked)" << std::endl;
#endif
	}
}

void HashTable::identifyMergeTasks_V2(tbb::concurrent_vector<std::pair<std::pair<point*, point*>, bool>> &globalMergeTasks)
{
  for (auto & coreBucket : coreBuckets)
    {
      point* core = coreBucket.representative; assert(core != NULL);
      std::vector<point *> toGetMerged;

      for (const auto element : coreBucket.members)
	{
	  if (element->isCore())
	    {
	      if ( (distanceFunc(*core, *element) <= epsilon)
		   &&
		   (core->findRoot() != element->findRoot())
		   )
		globalMergeTasks.push_back({{&*core, &*element}, false});
	    }
	}
    }
}

void HashTable::identifyAndPerformMergeTasks()
{
  for (auto & coreBucket : coreBuckets)
    {
      point* core = coreBucket.representative; assert(core != NULL);

      for (const auto element : coreBucket.members)
	{
	  if ( element->isCore() ) 
	    {
	      if ( (core->findRoot() != element->findRoot())
		   &&
		   (distanceFunc(*core, *element) <= epsilon)
		   )
		element->link(core);
	    }
	}
    }
}


void HashTable::initializeHashTable(size_t numberOfHyperplanes_)
{
  	for (size_t i = 0; i < numberOfHyperplanes_; i++)
	{
      	Hyperplane h;

      	for (size_t d = 0; d < ds->numberOfDimensions; d++)
		{
	  		h.features.push_back(gen->getRandVal());
		}

      	//h.normalize(false);/*without aux dimension*/
      	hyperplanes.push_back(h);
    }
}

void HashTable::initializeHashTable(std::string& fileName)
{
  std::string line;
  
  std::ifstream f(fileName);

  if (!f)
    {
      std::cerr << "Could not initialize hash table from file" << std::endl;
      exit(EXIT_FAILURE);
    }
  
  while (std::getline(f, line))
    {
      std::istringstream strStream(line);
      Hyperplane h(strStream);
      hyperplanes.push_back(h);

      if (h.features.size() != ds->numberOfDimensions)
	{
	  std::cerr << "Error while reading hyperplanes from file: "
		    << "Inconsistent dimensionality between points and hyperplanes.\n";
	  exit(EXIT_FAILURE);
	}
    }
}

void HashTable::performMergeTasks()
{
  for (const auto & task : mergeTasks)
    {
      task.first->link(task.second);
    }
}

std::ostream & HashTable::print(std::ostream & stream, char deli) const
{
  std::cout << "Hashed dataset's size is: " << hashedPoints.size() << std::endl;

  for ( const auto & p : this->hashedPoints)
    p.print(stream, deli) << std::endl;

  return stream;
}

std::ostream& HashTable::printCoreBuckets(std::ostream& stream, char deli) const
{
  size_t num = 0;
  for (const auto & coreBucket : coreBuckets)
    {
      stream << "Corebucket " << std::setw(2) << num++ << ":" << deli;
      for (const auto & element : coreBucket.members)
	{
	  stream << std::setw(3) << element->id
		 << std::setw(3) << (element->isCore() ? "(c) " : " ");
	}
      stream << std::endl;
    }
  return stream;
}

std::ostream& HashTable::printTable(std::ostream& stream, char deli) const
{
  size_t num = 0;

  for (const auto & element : myMap)
    {
      stream << "Bucket " << std::setw(2) << num++ << ":" << std::endl;

      for (const auto & p : element.second)
	{
	  p->print(stream << '\t', ' ') << std::endl;
	}	  
    }
  return stream;
}

std::ostream& HashTable::printMergeTasks(std::ostream& stream, char deli) const
{
  size_t num = 0;

  for (const auto & element : mergeTasks)
    {
      stream << "Merge task:" << deli
	     << std::setw(4)  << element.first->id << "&"
	     << std::setw(4)  << element.second->id << std::endl;
    }
  return stream;
}

void HashTable::identifyCoreBuckets_densityStyle()
{
  for (auto & element : myMap)
    {
      if (element.second.size() >= minPts)
	{
	  std::vector<point*> candidates, coreBucketElements;
	  candidates.reserve(element.second.size()); coreBucketElements.reserve(element.second.size());

	  std::copy(element.second.begin(), element.second.end(),
		    std::back_inserter(candidates));

	  point* core = candidates[getClosestToMean(candidates)];
	  //point* core = candidates[getMedian(candidates)]; #alternative
	  size_t cnt = 0;
	  
	  std::for_each(candidates.begin(),
			candidates.end(),
			[core, this, &cnt, &coreBucketElements](point* p)
			{
			  if (distanceFunc(*core, *p)  <= epsilon)
			    {
			      cnt++;
			      coreBucketElements.push_back(p);
			    }
			}
			);
	  
	  if (cnt >= minPts)
	    {
	      CoreBucket coreBucket;
	      coreBucket.representative = core;
	      coreBucket.members = coreBucketElements;
	      
	      core->setAsCorePoint();
	      std::for_each(coreBucketElements.begin(),
			    coreBucketElements.end(),
			    [core, this, &cnt](point* p)
			    {
			      if ( (p->findRoot() == p) )
				{
				  p->link(core);
				}
			    }
			    );
	      coreBuckets.push_back(coreBucket);
	    }
	}
    }
}

size_t HashTable::getClosestToMean(std::vector<point*> & vec) const
{
  assert(vec.size() != 0);
  std::vector<float> p(ds->numberOfDimensions, 0);
  point mean(p, 0);

  for (const point * p : vec)
    mean += *p;
  mean /= vec.size();

  std::vector<double> distanceToMean;
  distanceToMean.reserve(vec.size());

  std::transform(vec.begin(), vec.end(),
		 std::back_inserter(distanceToMean),
		 [&mean](const point* p)->double
		 {
		   return distanceFunc(*p, mean);
		 }
		 );

  size_t result = 0;
  double smallestVal = distanceToMean[0];

  for (size_t i = 1; i < vec.size(); ++i)
    {
      if (distanceToMean[i] < smallestVal)
	{
	  smallestVal = distanceToMean[i];
	  result = i;
	}
    }
  return result;
}

std::vector<point*> HashTable::getEpsNeighbours(point &query) /* used in ValillaDBSCAN+LSH */
{
  std::vector<point*> result;

  HashedPoint hashedPoint;

  std::transform(hyperplanes.cbegin(), hyperplanes.cend(),
		 std::back_inserter(hashedPoint.features),
		 [&query](const Hyperplane & h)->double
		 {
		   return hashFunc(query, h);
		 });

  
  return getEpsNeighbours(query, hashedPoint);
}


std::vector<point*> HashTable::getEpsNeighbours(point &query, HashedPoint &hashedPoint) /* used in ValillaDBSCAN+LSH */
{
  std::vector<point*> result;

  for (auto neighbour : myMap[hashedPoint])
    {
      if (distanceFunc(query, *neighbour) <= epsilon)
	{
	  result.push_back(neighbour);
	}
    }
  
  return result;
}
