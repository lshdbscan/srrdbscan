// #ifndef __SRR_LSHDBSCAN_H__
// #define __SRR_LSHDBSCAN_H__
#pragma once
#include <point.h>
#include <randomGen.h>
#include <transformation.h>
#include <statistics.h>
#include <unordered_set>
#include <chrono>
#include <globals.h>
#include <cassert>
//#include <highfive/H5File.hpp>

extern Statistics counters;

class PopulationTask
{
public:
    HashTable & table;
    bool booked = false;
    PopulationTask(HashTable & table_): table(table_)
    {}
};



class CorePointPopulationTask{
public:
    HashTable & table;
    bool booked = false;
    CorePointPopulationTask(HashTable & table_): table(table_)
    {}
};

class RangeTask
{
public:
    std::vector<point>::iterator begin, end;
    bool booked = false;

    RangeTask(std::vector<point>::iterator begin_,
		std::vector<point>::iterator end_):begin(begin_), end(end_)
    {}
};

class RangeTask_tbb
{
public:
    tbb::concurrent_vector<point*>::iterator begin, end;
    bool booked = false;

    RangeTask_tbb(tbb::concurrent_vector<point*>::iterator begin_,
		tbb::concurrent_vector<point*>::iterator end_):begin(begin_), end(end_)
    {}
};


class SRR_LSHDBSCAN{
private:
    //Given as input
    dataset* ds = NULL;
    double memoConstraint = 0.0; //Giga bytes
    double delta = 0.1; // between 0.0 and 0.5
    int level = -1; // fix the level? -1 means unfixed.
    double shrinkageFactor = 1.0;
    const bool benchmark = false;
    
    //Constructed at runtime
    double p1 = 0.0, p2 = 0.0;
    
    size_t maxDepth = 0;
    size_t EFP = 1; //Expected number of far away points for choosing lowest level  
    std::vector<std::vector<HashTable*>*> levels; // levels share hash functions, look at creating custom constructore that
                                                    // just concatenates 1 hash function to the previous one 
    RandGenerator gen;

    const size_t numberOfThreads; //#threads to be created
    pthread_t *pid; //threads reference
    
    void setDataset(dataset*);
    void setParams(); // Minimize C under memo contraints, as by notes from MATTI KARPPA (page 8) 
    void initLevels();
    void initLevelsCorePoints();
    double memoCost(size_t maxDepth);

    size_t reps(size_t k);
    std::ostream& getWorkPoint(std::ostream& stream, char deli, point q);
    std::vector<size_t> getWorkPoint(point q);
    size_t findBestLevel(point q); //as described in LSH_SRR paper, Algorithm 1  
    void populateHashTables();  //Intanciate all K levels and their hashtables
    void populateHashTables_threaded();
    
    void identifyCorePoints();  //Check all points if they are a CP using the best Level  
    void identifyCorePoints_threaded();

    void populateCorePoints_threaded();
    size_t findCPIdentificationLevel();
    size_t findCPMergingLevel();
    void CPMerging_threaded();
    void identifyBorderPoints_threaded();

    void labelPoints();
    
    std::ofstream benchStream;

    std::chrono::time_point<std::chrono::steady_clock> start, stop;
    
    std::chrono::duration<double> duration_initializingHashTables;
    std::chrono::duration<double> duration_populatingHashTables;
    std::chrono::duration<double> duration_identifyingCorePoints;
    std::chrono::duration<double> duration_populatingHashTables_CP;
    std::chrono::duration<double> duration_identifyingBorderPoints;
    std::chrono::duration<double> duration_Merging;
    std::chrono::duration<double> duration_relabelingData;

    
    
//Depricated
size_t maxLevel(double p2);  // Find largest k such that reps(k) <= L
                                //I need to know the sensitivy before finding this value is possible
    

public:
    std::vector<PopulationTask> populationTasks; //the different buckets over all levels, that the threads have to populate
    std::vector<RangeTask> cpIdentifycationTasks;
    std::vector<CorePointPopulationTask> CPPopulationTasks;
    std::vector<PopulationTask> CPMergingTasks;
    std::vector<RangeTask_tbb> bpIdentificationTasks;

    size_t cpLevel = 0;
    
    tbb::concurrent_vector<point*> corePoints; 
    tbb::concurrent_vector<point*> possibleBorderPoints; 
    static void* populateHashTables_thread(void*);
    static void* populateCorePoints_thread(void*);
    static void* cpIdentifycation_thread(void*);
    static void* CPMerge_thread(void *);
    static void* bpIdentifycation_thread(void *);
    

    void introduceMe();
    SRR_LSHDBSCAN(dataset *ds, double delta, double GBytes, bool benchmark, std::string benchName , size_t numberOfThreads = 2, 
        int level = -1, double shrinkageFactor=1.0); //epsilon might not be needed for the class, just set the value globally like they did
    void performClustering();
    std::ostream& getCorePoints(std::ostream&, char deli) const;
    std::ostream& getLabels(std::ostream&, char deli) const;
    void writeHDF5(std::string fileName, Statistics& counters);

    std::ostream& getBenchmarkResults(std::ostream&, char deli) const;
    void getWork(std::ostream& stream,char deli);
    void getWork(std::string fileName);
    ~SRR_LSHDBSCAN();
}; 

  

// #endif