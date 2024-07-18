#include <SRR_LSHDBSCAN.h>

long g_finished = 0;

//#define SRR_PERFORMANCE

Statistics counters;

SRR_LSHDBSCAN::SRR_LSHDBSCAN(dataset *ds_, 
                            double delta_,
                            double memoConstraint_,
                            bool benchmark_,
                            std::string benchName,
                            size_t numberOfThreads_,
                            int level_,
                            double shrinkageFactor_):delta(delta_),memoConstraint(memoConstraint_), 
                              numberOfThreads(numberOfThreads_) ,benchmark(benchmark_), benchStream(benchName), 
                              level(level_), shrinkageFactor(shrinkageFactor_)
{
  	for (int i = 0; i < 30; i++) {
		  HASH.push_back(gen.getHashCoeff());
	  }
    setDataset(ds_);
    setParams();  

    pid = new pthread_t[numberOfThreads];

    //initialize hashes

    //std::cerr << "depth: " << maxDepth << std::endl;
    initLevels();
    size_t numPoints = ds->points.size();
    std::cerr << "data points " << numPoints << std::endl;
    size_t noBatches = 100;
    assert(numPoints>= noBatches);
    size_t stepsize = (numPoints + noBatches - 1) / noBatches; //ceil  
    cpIdentifycationTasks.reserve(noBatches);
    std::vector<point>::iterator start     =  this->ds->points.begin();

    for(size_t i  = 0; i < (noBatches-1); i++){ //current idea is just to split the work equally over all available threads.
      cpIdentifycationTasks.emplace_back(start, start + stepsize);
      start     += stepsize;
    }
    cpIdentifycationTasks.emplace_back(start, this->ds->points.end());
}

SRR_LSHDBSCAN::~SRR_LSHDBSCAN(){
  //std::cerr << "Destructor being called" << std::endl;
  for(size_t k = 0; k < levels.size(); k++){
    for(auto T: *levels[k]){
      delete T;
    }
    delete levels[k];
  }
  delete pid;
  //std::cerr << "Destructor finished" << std::endl;
}

void SRR_LSHDBSCAN::setDataset(dataset* d)
{
  this->ds = d;
  if (ds == NULL)
    {
      std::cerr << "Invalid Dataset." << std::endl;
      exit(EXIT_FAILURE);
    }
}

double SRR_LSHDBSCAN::memoCost(size_t limit){
  double hashCost = ((1 - std::pow(1/p1, limit+1.0)) / (1.0 - 1/p1));
  std::cout << "hashCost: " << hashCost << std::endl;
  double r_ = 1.0/p1; // current this is the same as e as (e^-1)^-1 = e. However, this is more future proof as if we decide to normalise differently in the hashfunctions.
  double HPCost = ((r_ - (limit+1) * std::pow(M_E, limit+1) + limit * std::pow(M_E, limit + 2))/pow(1-M_E, 2));
  double total_cost = ((4 * ds->points.size() * std::log(1/delta) * hashCost) /1e9) + ((8 * std::log(1/delta)* HPCost)  / 1e9); 
  std::cerr << "level: " << limit << " cost: " << total_cost << std::endl;
  //std::cerr << "\t hashCost: " << hashCost << " HPCost: " << HPCost << std::endl;
  return total_cost;
}

//This should not be done in K_max rounds
void SRR_LSHDBSCAN::setParams(){
    p1 = .8; //1/M_E; //.8 is setting for 4 * epsilon choice.
    p2 = 0.0; //This values does not really matter until we care what the approximation factor(C) is.  
    size_t K_max = 0; 

    while(memoCost(K_max) < memoConstraint){
        std::cerr << "Level: " << K_max << " Cost: " << memoCost(K_max) << std::endl;
        HashTable ht(ds, K_max, &gen);
        ht.populateHashTable();
        if (ht.hashTable.size() > .1 * ds->points.size()) {
          std::cerr << "not building more levels since buckets are too small" << std::endl;
          break;
        }
        K_max++;
    }
    if(K_max > 0) this->maxDepth =  --K_max; //size_t could be unsigned int and (--) would underflow
    benchStream << "MemoConsumption: " <<  maxDepth + 1 <<  " levels will consume " << memoCost(maxDepth) << std::endl; // maxDepth + 1 since we have level 0

    return;
}


void SRR_LSHDBSCAN::initLevels(){
  std::cerr << "initting levels" << std::endl;
  for(size_t k = 0; k <= maxDepth; k++){
    std::cerr << shrinkageFactor << std::endl;
    size_t repetitions = ceil(reps(k) * log(1/delta) * shrinkageFactor);
    if(k == 0) repetitions = 1;
    std::cerr << k  << " " << reps(k) << " " << log(1/delta) << " " << repetitions; 
    std::vector<HashTable*> *levelK = new std::vector<HashTable*>();
    for(size_t T = 0; T < repetitions; T++){
      levelK->push_back(new HashTable(ds, k, &gen));
      populationTasks.emplace_back(*(levelK->back()));
    }
    std::cerr << " " << levelK->size() << std::endl;
    levels.push_back(levelK);
    std::cerr << "Done with level: " << k << " with reps(level): " << repetitions << std::endl; 
  }
}

void SRR_LSHDBSCAN::initLevelsCorePoints(){
  //std::cerr << "initting levels" << std::endl;
  for(size_t k = 0; k <= maxDepth; k++){
    size_t repetitions = ceil(reps(k) * log(1/delta) * shrinkageFactor);
    if(k == 0) repetitions = 1;
      //std::cerr << k  << " " << reps(k) << " " << log(1/delta) << " " << repetitions; 
      std::vector<HashTable*> *levelK = new std::vector<HashTable*>();
      for(size_t T = 0; T < repetitions; T++){
        levelK->push_back(new HashTable(ds, k, &gen));
        CPPopulationTasks.emplace_back(*(levelK->back()));
    }
    //std::cerr << " " << levelK->size() << std::endl;
    levels.push_back(levelK);
    //std::cerr << "Done with level: " << k << " with reps(level): " << repetitions << std::endl; 
  }
}

void SRR_LSHDBSCAN::populateHashTables(){
  for(auto level: levels){
    for(auto table: *level){
      table->populateHashTable();
    }
    std::cerr << "Done with level: " << level << " with: " << (*level).size() << " tables" << std::endl;
  }
}

void SRR_LSHDBSCAN::populateHashTables_threaded()
{
  benchStream << "Populating the concurrent way! \t"; // << std::endl;
  std::cout << "Populating the hash tables!" << std::endl;

  for (size_t threadID = 0; threadID < numberOfThreads; ++threadID)
    {
      pthread_create(&pid[threadID], NULL, populateHashTables_thread, this);
    }

  for (size_t threadID = 0; threadID < numberOfThreads; ++threadID)
    {
      pthread_join(pid[threadID], NULL);
    }
}

void SRR_LSHDBSCAN::identifyCorePoints_threaded(){
  benchStream << "Identifying the core points the concurrent way! \t";// << std::endl;
  for (size_t threadID = 0; threadID < numberOfThreads; ++threadID)
    {
      pthread_create(&pid[threadID], NULL, cpIdentifycation_thread, this);
    }

  for (size_t threadID = 0; threadID < numberOfThreads; ++threadID)
    {
      pthread_join(pid[threadID], NULL);
    }
}

void SRR_LSHDBSCAN::populateCorePoints_threaded(){
 benchStream << "Populating Core Points the concurrent way! \t";// << std::endl;
  for (size_t threadID = 0; threadID < numberOfThreads; ++threadID)
    {
      pthread_create(&pid[threadID], NULL, populateCorePoints_thread, this);
    }

  for (size_t threadID = 0; threadID < numberOfThreads; ++threadID)
    {
      pthread_join(pid[threadID], NULL);
    }
}

void SRR_LSHDBSCAN::identifyBorderPoints_threaded(){
  benchStream << "Identifying the border points the concurrent way! \t";// << std::endl;
  for (size_t threadID = 0; threadID < numberOfThreads; ++threadID)
    {
      pthread_create(&pid[threadID], NULL, bpIdentifycation_thread, this);
    }

  for (size_t threadID = 0; threadID < numberOfThreads; ++threadID)
    {
      pthread_join(pid[threadID], NULL);
    }
}

void SRR_LSHDBSCAN::CPMerging_threaded(){
  benchStream << "Merging Core points the concurrent way! \t";// << std::endl;
  for (size_t threadID = 0; threadID < numberOfThreads; ++threadID)
    {
      pthread_create(&pid[threadID], NULL, CPMerge_thread, this);
    }

  for (size_t threadID = 0; threadID < numberOfThreads; ++threadID)
    {
      pthread_join(pid[threadID], NULL);
    }
}



void* SRR_LSHDBSCAN::populateHashTables_thread(void* inputArg){
  auto input = (SRR_LSHDBSCAN*) inputArg;
  for ( auto & task : input->populationTasks){
    if (task.booked == false){
      if (__sync_bool_compare_and_swap ( &(task.booked),
                false,
                true)){
        task.table.populateHashTable();
      }
    }
  }  
  return NULL;
}

void* SRR_LSHDBSCAN::populateCorePoints_thread(void* inputArg){
  auto input = (SRR_LSHDBSCAN*) inputArg;
  for ( auto & task : input->CPPopulationTasks){
    if (task.booked == false){
      if (__sync_bool_compare_and_swap ( &(task.booked),
                false,
                true)){
        task.table.populateHashTable(input->corePoints.begin(), input->corePoints.end());
      }
    }
  }  
  return NULL;
}

void* SRR_LSHDBSCAN::bpIdentifycation_thread(void* inputArg){
  auto input = (SRR_LSHDBSCAN*) inputArg;
  for (auto & task: input->bpIdentificationTasks){
    if (task.booked == false){
      if (__sync_bool_compare_and_swap ( &(task.booked),
                false,
                true)){
        
        for (auto iter = task.begin; iter < task.end; iter++){
          //Performs the CP identityfication for each point that the thread is responsible for.
          point p = **iter;
          bool outerBreak = false;
          int bestK = input->findBestLevel(p);
          for(auto T = input->levels[bestK]->begin(); T != input->levels[bestK]->end(); T++){
            //Hash q with the hyperplanes of T
            HashedPoint hq;
            if(outerBreak)break;
            std::transform(
            (*T)->hyperplanes.cbegin(), 
            (*T)->hyperplanes.cend(),
            std::back_inserter(hq.features),
            [&p](const Hyperplane &h) -> double
            {
                  return hashFunc(p, h);
            });      
            for(auto q : (*T)->hashTable[hq.combine()]){
              if(p.squaredEuclideanDistance(*q) <= epsilon){//Since only CP are part of the tables we know p is a borderPoint
                (*iter)->link(q);
                outerBreak = true;
                break;
              }
            }
          }
		    }
      }
    }
  }
 return NULL;
}

//we need a populate Hashtables with Corepoints only function
//safe plan is just to use the exact same property that we saw before
//we are bound by finding the CP since the new hashtable is at most N big (not tight, should be decently less work)

void* SRR_LSHDBSCAN::cpIdentifycation_thread(void* inputArg){
  auto input = (SRR_LSHDBSCAN*) inputArg;
  for (auto & task: input->cpIdentifycationTasks){
    if (task.booked == false){
      if (__sync_bool_compare_and_swap ( &(task.booked),
                false,
                true)){
        
        std::unordered_set<point*> result_set;
        result_set.reserve(minPts);
        for (auto iter = task.begin; iter < task.end; iter++){
          //Performs the CP identityfication for each point that the thread is responsible for.
          point p = *iter;
          result_set.insert(&p);
          int bestK = input->level;
          if (bestK == -1) {
            bestK = input->cpLevel;
          } 
          long comparisons = 0;
          long truepoints = 0;
          //std::cout << "identification on level " << bestK << std::endl;
          //std::unordered_set<point*> result_set; result_set.insert(&p);
          for(auto T = input->levels[bestK]->begin(); T != input->levels[bestK]->end(); T++){
            //Hash q with the hyperplanes of T
            HashedPoint hq;
            std::transform(
                (*T)->hyperplanes.cbegin(), 
                (*T)->hyperplanes.cend(),
                std::back_inserter(hq.features),
                [&p](const Hyperplane &h) -> double
                {
                  return hashFunc(p, h);
                }); 
            for(auto q : (*T)->hashTable[hq.combine()]){     
              if (result_set.size() >= minPts)
                break;
              comparisons++;
              if(p.squaredEuclideanDistance(*q) <= epsilon){
                truepoints++;
                result_set.insert(q);
              }
            }
            if(result_set.size() >= minPts){
              iter->setAsCorePoint();
              input->corePoints.push_back(&(*iter)); //todo: is this still needed?
              result_set.clear();
              break;
            }
          }
#ifdef SRR_PERFORMANCE
          std::cout << "carried out " << comparisons << " comparisons on level " << bestK << " and found " << truepoints << " passing points. Resulted in " << result_set.size() << " points." << std::endl;
#endif
          g_finished++;
          if (g_finished % 10000 == 0) {
            std::cout << "finished roughly " << g_finished << " points" << std::endl;          
          }
          if(0 < result_set.size() && result_set.size() < minPts){
            input->possibleBorderPoints.push_back(&(*iter));
          }
		    }
      }
    }
  }
 return NULL;
}

void* SRR_LSHDBSCAN::CPMerge_thread(void* inputArg){
  auto input = (SRR_LSHDBSCAN*) inputArg;
  for (auto & task: input->CPMergingTasks){
    if (task.booked == false){
      if (__sync_bool_compare_and_swap ( &(task.booked),
                false,
                true)){
        task.table.mergeCorePoints();    
      }
    }    
  }
  return NULL;
}



size_t SRR_LSHDBSCAN::reps(size_t k){
    return std::ceil(1.0 /std::pow(p1, k));
}

size_t SRR_LSHDBSCAN::findBestLevel(point q){
    size_t k = 1, kbest = 0, w_kbest = ds->points.size();
    size_t L = reps(maxDepth); //we should be able to substitude L for max reps, since we dont get what L is;
    while(reps(k) <= std::min(L, w_kbest)){
        //std::cerr << "computing work for level: " << k << std::endl; 
        size_t w_k = 0;
        for(auto T = levels[k]->begin(); T != levels[k]->end(); T++){
          //Hash q with the hyperplanes of T
          HashedPoint hq;
          std::transform(
          (*T)->hyperplanes.cbegin(), 
          (*T)->hyperplanes.cend(),
          std::back_inserter(hq.features),
          [&q](const Hyperplane &h) -> double
          {
                return hashFunc(q, h);
          });
          
          //Add |T[hash_T(q)]| to w_k
          w_k += 1 + (*T)->hashTable[hq.combine()].size();
        }

#ifdef SRR_PERFORMANCE
        std::cout << "Work on level " << k << " is " << w_k << std::endl;
#endif

        if(w_k < w_kbest){
          kbest = k; w_kbest = w_k;
        }
        //std::cerr << "Done for level: " << k << std::endl;
        k++;
    }
    return kbest;
}

size_t SRR_LSHDBSCAN::findCPIdentificationLevel(){ //This can maybe be optimized to not look further at a given point
  
  size_t bestlevel = 0, minWork = ds->points.size() * ds->points.size(), index = 0;
  for(auto level = levels.begin(); level != levels.end(); level++){
    size_t totalWork_level = 0;
    for(auto table = (*level)->begin(); table != (*level)->end(); table++ ){
      for(auto bucket: (*table)->hashTable){
        size_t bucket_size = bucket.second.size(); 
        totalWork_level += 1 + (bucket_size * (bucket_size - 1)) / 2;
      }
      if(totalWork_level > minWork){
        break;
      }
    }
#ifdef SRR_PERFORMANCE
    std::cout << "Level " << index << " has total work " << totalWork_level << " per table" << std::endl;
#endif
    if(totalWork_level < minWork){
      minWork = totalWork_level;
      bestlevel = index;
    }
    index++;
  }
#ifdef SRR_PERFORMANCE
  std::cout << "Chosen level: " << bestlevel << std::endl; 
#endif
  return bestlevel;
}

double harmonic_mean(std::vector<uint64_t> arr) { 
  double res = 0.0;
  for (uint64_t& x: arr) {
    res += 1/(double)x;
  }
  return arr.size() / res;
}

uint64_t average(std::vector<uint64_t> arr) {
  uint64_t res = 0;
  for (auto& x: arr) {
    res += x;
  }
  return res / arr.size();
}

size_t SRR_LSHDBSCAN::findCPMergingLevel(){ //This can maybe be optimized to not look further at a given point
  
  size_t bestlevel = 0, minWork = corePoints.size() * corePoints.size(), index = 0;
  std::vector<std::vector<uint64_t>> totalWork; 
  for(auto level = levels.begin(); level != levels.end(); level++){
    size_t totalWork_level = 0;
    std::vector<uint64_t> workonlevel;
    for(auto table = (*level)->begin(); table != (*level)->end(); table++ ){
      uint64_t workintable = 0;
      for(auto bucket: (*table)->hashTable){
        size_t bucket_size = bucket.second.size(); 
        totalWork_level += 1 + (bucket_size * (bucket_size - 1)) / 2;
        workintable += 1 + (bucket_size * (bucket_size - 1)) / 2;
      }
      workonlevel.push_back(workintable);
      if(totalWork_level > minWork){
        break;
      }
    }

    std::sort(workonlevel.begin(), workonlevel.end());
    totalWork.push_back(workonlevel);
    // std::cout << "reps:" << workonlevel.size() << std::endl;

    // for (auto& x: workonlevel) { 
    //   std::cout << x << " ";
    // }
    // std::cout << std::endl;

    uint64_t median_totalwork = (*level)->size() * workonlevel[workonlevel.size() / 2];
    //uint64_t median_totalwork = totalWork_level;

    // std::cout << median_totalwork << std::endl;    


#ifdef SRR_PERFORMANCE
    std::cout << "Level " << index << " has total work " << totalWork_level << " per table" << std::endl;
    std::cout << "It has median work of " << workonlevel.size() * workonlevel[workonlevel.size() / 2] << std::endl; 
    //std::cout << "It has harmonic average work of " << workonlevel.size() * harmonic_mean(workonlevel) << std::endl;
    //std::cout << "It has average work of " << workonlevel.size() * average(workonlevel) << std::endl;
#endif
    // if(totalWork_level < minWork){
    //   minWork = totalWork_level;
    //   bestlevel = index;
    // }
    if(median_totalwork < minWork){
      minWork = median_totalwork;
      bestlevel = index;
    }
    index++;
  }
#ifdef SRR_PERFORMANCE
  std::cout << "Chosen level: " << bestlevel << std::endl; 
#endif

  // for (auto& w: totalWork) {
  //   std::sort(w.begin(), w.end());
  // }
  // for (auto& w: totalWork[bestlevel]) { 
  //   std::cout << " " << w;
  // }
  // std::cout << std::endl;
  // auto bestlevelwork = totalWork[bestlevel];
  // auto n_bestlevel = bestlevelwork.size();
  // std::cout << "Median-to-Max-Ratio:" << bestlevelwork[n_bestlevel - 1] / (double) bestlevelwork[n_bestlevel / 2] << std::endl;
  return bestlevel;
}


std::ostream& SRR_LSHDBSCAN::getWorkPoint(std::ostream& stream,char deli, point q){
 
  size_t k = 0, kbest = 0;

  size_t L = reps(maxDepth); //we should be able to substitude L for max reps, since we dont get what L is;
  //stream << w_kbest << deli;
  while(reps(k) <= L){
      //std::cerr << "computing work for level: " << k << std::endl; 
      size_t w_k = 0;
      for(auto T = levels[k]->begin(); T != levels[k]->end(); T++){
        //Hash q with the hyperplanes of T
        HashedPoint hq;
        std::transform(
        (*T)->hyperplanes.cbegin(), 
        (*T)->hyperplanes.cend(),
        std::back_inserter(hq.features),
        [&q](const Hyperplane &h) -> double
        {
              return hashFunc(q, h);
        });
        
        //Add |T[hash_T(q)]| to w_k
        w_k += 1 + (*T)->hashTable[hq.combine()].size();
      }
      stream << w_k << deli;
      k++;
    }
    return stream;
}

std::vector<size_t> SRR_LSHDBSCAN::getWorkPoint(point q){

  size_t k = 0, kbest = 0;
  size_t L = reps(maxDepth); //we should be able to substitude L for max reps, since we dont get what L is;
  //stream << w_kbest << deli;
  std::vector<size_t> totalWork;
  while(reps(k) <= L){
      //std::cerr << "computing work for level: " << k << std::endl; 
      size_t w_k = 0;
      for(auto T = levels[k]->begin(); T != levels[k]->end(); T++){
        //Hash q with the hyperplanes of T
        HashedPoint hq;
        std::transform(
        (*T)->hyperplanes.cbegin(), 
        (*T)->hyperplanes.cend(),
        std::back_inserter(hq.features),
        [&q](const Hyperplane &h) -> double
        {
              return hashFunc(q, h);
        });

                //Add |T[hash_T(q)]| to w_k
        w_k += 1 + (*T)->hashTable[hq.combine()].size();
      }
      totalWork.push_back(w_k);
      k++;
    }
    return totalWork;
}


void SRR_LSHDBSCAN::getWork(std::ostream& stream,char deli){ // this function is only for experiments - serves no purpose for actually solving the problem.
  populateHashTables_threaded();
  populationTasks.clear(); 
  
  std::cout << "done populating hashtables" << std::endl;
  
  for(auto p: ds->points){
    getWorkPoint(stream,deli,p) << std::endl; 
  }
}

void SRR_LSHDBSCAN::getWork(std::string fileName){
  //populateHashTables_threaded();
  //populationTasks.clear(); 
  //std::cout << "done populating hashtables" << std::endl;
  
  // HighFive::File file(fileName, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
  // HighFive::DataSet dataset = file.createDataSet<size_t>("work", {this->ds->points.size(), maxDepth+1}); 
  // size_t counter = 0;
  // for(auto p: ds->points){
  //   dataset.select({counter,0}, {1, maxDepth + 1}).write(getWorkPoint(p));
  //   counter++; 
  //   if (counter % 1000 == 0) {
  //     std::cout << counter << std::endl;
  //   }
  // }
}

void SRR_LSHDBSCAN::identifyCorePoints(){
  for(auto &p : ds->points){
    size_t bestK = findBestLevel(p);
    std::unordered_set<point*> result_set; result_set.insert(&p);

    for(auto T = levels[bestK]->begin(); T != levels[bestK]->end(); T++){
          //Hash q with the hyperplanes of T
          HashedPoint hq;
          std::transform(
          (*T)->hyperplanes.cbegin(), 
          (*T)->hyperplanes.cend(),
          std::back_inserter(hq.features),
          [&p](const Hyperplane &h) -> double
          {
                return hashFunc(p, h);
          });

          for(auto q : (*T)->hashTable[hq.combine()]){
            if(p.squaredEuclideanDistance(*q) <= epsilon){
              result_set.insert(q);
            }
          }
          //Unsure if breaking early will break the structure, as all points have to be
          //clustered to the CP
          //if(result_set.size() >= minPts){
            //p.setAsCorePoint();
          //  corePoints.push_back(&p);
          //  break;                        
          //}
        }
    //we should use UF here to keep track of density reachable points
    //Since all corePoints see their epsilon neighbors, a borderPoint will be assigned to a cluster
    //If it is a borderPoint. In short we only have to union the corePoints with their epsilon neighborhood
    //to ensure a true clustering of the dataset.
    if(result_set.size() >= minPts){
      p.setAsCorePoint();
      for(auto q: result_set){
        p.link(q);
      }
    }
  }
}


void SRR_LSHDBSCAN::labelPoints(){
  //std::cout << "starting labeling" << std::endl;
  for(auto &p: ds->points){
    //if(p.isCore()){ std::cout << p.id << std::endl;}
    p.reLabel();
  }
}

std::ostream& SRR_LSHDBSCAN::getCorePoints(std::ostream& stream, char deli) const {
  
  for(auto p: corePoints){
    p->print(stream, deli, true) << std::endl;
  }
  return stream;
}

std::ostream& SRR_LSHDBSCAN::getLabels(std::ostream& stream, char deli) const{
  for(auto p: ds->points){
    p.print(stream, deli, true) << std::endl;
  }
  return stream;
}

void SRR_LSHDBSCAN::introduceMe()
{
  benchStream << "Spherical range reporting LSHDBSCAN:\n"
	    << "\t#points: " << this->ds->points.size()
	    << "\t#dims: " << this->ds->numberOfDimensions
    	    << "\tmetric: " << ( (metric == angular) ? "angular" : "euclidean")
	    << std::endl
	    << "\tDelta: " << this->delta
	    << "\tMemory Constraint in GB: " << this->memoConstraint 
	    << std::endl
	    << "\tEps: " << epsilon_original
	    << "\tminPts: " << minPts
	    << std::endl;
}



void SRR_LSHDBSCAN::performClustering(){
  
  benchStream << "Populating hash tables:  "; 

  
  start = std::chrono::steady_clock::now();
  auto total_start = start;
  populateHashTables_threaded();
  populationTasks.clear(); //How to ensure this happens after above function terminates?
  // populateHashTables();
  stop  = std::chrono::steady_clock::now();
  duration_populatingHashTables = stop - start;
  benchStream << duration_populatingHashTables.count() << std::endl;
  counters.add_measurement("build ht", duration_populatingHashTables.count());

  std::cout << "Writing statistics about work" << std::endl;
  //getWork("test.hdf5");
  std::cout << "done writing statistics" << std::endl;

  benchStream << "Identifying core points:  " << std::flush; 
  start = std::chrono::steady_clock::now();
  //identifyCorePoints();
  cpLevel = findCPIdentificationLevel();
  std::cout << "best level is " << cpLevel; 
  identifyCorePoints_threaded();
  stop  = std::chrono::steady_clock::now();
  duration_identifyingCorePoints = stop - start;
  benchStream << duration_identifyingCorePoints.count() << std::endl;
  counters.add_measurement("identify cp", duration_identifyingCorePoints.count());
  //Delete the old multi level hash tables
  for(size_t k = 0; k < levels.size(); k++){
    for(auto T: *levels[k]){
      delete T;
    }
    delete levels[k];
    levels[k]->resize(0);
  }
  levels.resize(0);
  std::cout << "number of core points: " << corePoints.size() << std::endl;

  if (corePoints.size() > 0) {
  
    //construct new SRR with only the Core Points
    benchStream << "Populating tables with only CP:\t" << std::flush;
    
    start = std::chrono::steady_clock::now();  
    initLevelsCorePoints();
    benchStream << "initted" << std::endl;
    populateCorePoints_threaded();
    stop = std::chrono::steady_clock::now();
    duration_populatingHashTables_CP = stop - start;
    benchStream << duration_populatingHashTables_CP.count() << std::endl;

    benchStream << "Number of possibleBorderPoints: \t" << possibleBorderPoints.size() << std::endl;


    //see if any possible border point is within reach of a CP;
    //Link a border Point to only 1 CP (since it should not collapse 2 clusters into 1)
    if (possibleBorderPoints.size() > 0) { 
      benchStream << "Identifying border points: ";
      start = std::chrono::steady_clock::now();  
      size_t numPoints = possibleBorderPoints.size();
      size_t noBatches = std::min(numberOfThreads, numPoints);
      size_t stepsize = (numPoints + noBatches - 1) / noBatches; //ceil  
      cpIdentifycationTasks.reserve(noBatches);
      tbb::concurrent_vector<point*>::iterator starting     =  this->possibleBorderPoints.begin();

      for(size_t i  = 0; i < (noBatches-1); i++){ //current idea is just to split the work equally over all available threads.
        bpIdentificationTasks.emplace_back(starting, starting + stepsize);
        starting     += stepsize;
        //IFAIK starting going over end should be fine since this is checked in the threads work
      }
      bpIdentificationTasks.emplace_back(starting, this->possibleBorderPoints.end());
      
      //identifyBorderPoints_threaded();
      stop = std::chrono::steady_clock::now();  
      duration_identifyingBorderPoints = stop - start;
      benchStream << duration_identifyingBorderPoints.count() << std::endl;
      counters.add_measurement("identify bp", duration_identifyingBorderPoints.count());
    }

    //Union CP using the new multi level structure (only contains Core points)
    //Make tasks (select level based on sum of bucket sizes squared)

    benchStream << "Merging: ";
    start = std::chrono::steady_clock::now(); 
    size_t CPMergelevel;
    if (level != -1) {
      std::cout << "Using fixed level " << level << std::endl;
      CPMergelevel = level;
    } else {
      CPMergelevel = findCPMergingLevel(); 
      std::cout << "Merging on level " << CPMergelevel << std::endl;
    } 
    CPMergingTasks.reserve(reps(CPMergelevel)); 
    for(auto table = (*levels[CPMergelevel]).begin(); table != (*levels[CPMergelevel]).end(); table++){
      CPMergingTasks.emplace_back(**table);
    }  
    CPMerging_threaded();
    stop = std::chrono::steady_clock::now();
    duration_Merging = stop - start;
    benchStream << duration_Merging.count() << std::endl;
    counters.add_measurement("merging", duration_Merging.count());
  }
  
  benchStream << "labeling points: \t" << std::flush;
  start = std::chrono::steady_clock::now();
  labelPoints();
  stop = std::chrono::steady_clock::now();
  duration_relabelingData = stop - start;
  benchStream << duration_relabelingData.count() << std::endl;
  benchStream << "Total: " << (stop - total_start).count() / 1e9 << std::endl;
  counters.add_measurement("total", (stop - total_start).count() / 1e9);
}

//Depricated
size_t SRR_LSHDBSCAN::maxLevel(double P2_Tmp){
  return std::ceil(std::log(1.0 * minPts/ds->points.size())/std::log(P2_Tmp));
} 