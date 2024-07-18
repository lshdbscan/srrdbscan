//#pragma once
#include<globals.h>
#include<SRR_LSHDBSCAN.h>
#include<point.h>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

void read_hdf5(dataset& d, std::string fileName) {
  std::cout << "Reading HDF5 format input dataset..." << std::endl;
  HighFive::File file(fileName, HighFive::File::ReadOnly);
  auto hdf5Dataset = file.getDataSet("data");
  std::vector<std::vector<float>> data;
  hdf5Dataset.read(data);

  d.readData(data);
}

void write_result(dataset& ds, std::string fileName, Statistics& counters){
  HighFive::File file(fileName, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

  std::vector<int> labels;
  for(auto p: ds.points){
    labels.push_back(p.print());
  }

  for (auto &[key , val]: counters.stats) {
    file.createAttribute(key, val);
  }

  file.createDataSet("grp/labels", labels);
  return;
}

#define BENCHMARK true

int main(int argc, char* argv[])
{
  SRRParameters params;
  parseSRRArguments(argc, argv, params);

  dataset d;
  auto start = std::chrono::steady_clock::now();
  read_hdf5(d, params.fileName);
  auto stop = std::chrono::steady_clock::now();  
  std::chrono::duration<double> duration = (stop - start);
  std::cout << "Reading data took: " << duration.count() << std::endl;
  
  std::stringstream benchName;
  benchName << "Benchmark_file_" << d.name << "_eps_" << epsilon_original << "_minPts_" << minPts << "_delta_" << 
    params.delta << "_memoLimit_"<< params.memoConstraint << "_level_" << params.level << "_shrinkage_" << params.shrinkageFactor << ".txt";
  
  std::cout << "Saving in file " << benchName.str() << std::endl;

  if (metric == angular)
    {
      d.normalizeData();      //d.meanRemoveData();
    }

  std::cerr << "Memory constraint is " << params.memoConstraint << std::endl;
  SRR_LSHDBSCAN SRR_dbscan(&d,
        params.delta,
        params.memoConstraint,
				BENCHMARK,
        benchName.str(),      
        params.numberOfThreads,
        params.level,
        params.shrinkageFactor);
  SRR_dbscan.introduceMe();
  
  std::cout << "Running SSR_LSH clustering" << std::endl;
  SRR_dbscan.performClustering();
  std::cout << "Writing Labels to file" << std::endl;
  std::stringstream ss;
  ss << "Labels_file_" << d.name << "_eps_" << epsilon_original << "_minPts_" << minPts << "_delta_" << params.delta << "_memoLimit_"<< 
    params.memoConstraint << "_level_" << params.level << "_shrinkage_" << params.shrinkageFactor << ".h5"; 
  
  //SRR_dbscan.writeHDF5(ss.str(), counters);
  write_result(d, ss.str(), counters);

  std::cout << counters.stats["total"] << std::endl;



  // std::cout << "Writing CP to file" << std::endl;
  // std::stringstream ss;
  // ss << "CorePoints_eps_" << epsilon_original << "_minPts_" << minPts << "_delta_" << params.delta << ".txt"; 
  // std::ofstream CPStr(ss.str());
  // SRR_dbscan.getCorePoints(CPStr, '\t') << std::endl;
  // CPStr.close();


  // std::cout << "Writing work to file" << std::endl;
  // std::stringstream ss1;
  // ss1 << "Work_eps_" << epsilon_original << "_minPts_" << minPts << "_delta_" << params.delta << ".h5"; 
  // SRR_dbscan.getWork(ss1.str());


  // OLD STUFF  
  /*
  std::ofstream benchmarkStr("benchmark.txt", std::ofstream::out | std::ofstream::app);
  SRR_dbscan.getBenchmarkResults(benchmarkStr, '\t') << std::endl;
  benchmarkStr.close();

  std::cout << "Saving clustering labels..." << std::endl;
  std::string outputFile = params.fileName
    + "_"
    + std::to_string(params.numberOfHashTables)
    + "_"
    + std::to_string(params.numberOfHyperplanesPerTable)
    + "_"
    + std::to_string(params.numberOfThreads)
    + ".idx_concurrentlshdbscan";

  std::ofstream labelStr(outputFile);	  
  d.printData(labelStr, '\t', true);
  labelStr.close();

  if (!params.baselineFileName.empty())
    {
      std::cout << "Computing clustering accuracy against baseline: "
		<< params.baselineFileName
		<< " ..."
		<< std::endl;
      std::string accuracyCommand = std::string("./computeAccuracy_fast ")
	+ outputFile + " "
	+ params.baselineFileName + " "
	+ "accuracy.txt";

      system(accuracyCommand.c_str());
    }
  
  */

  return 0;
}
