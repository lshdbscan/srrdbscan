#include<globals.h>
#include<SRR_LSHDBSCAN.h>
#include<point.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace dbscan_srr {
namespace python {

namespace py = pybind11;


// int main(int argc, char* argv[])
// {
//   SRRParameters params;
//   parseSRRArguments(argc, argv, params);

//   dataset d;
//   auto start = std::chrono::steady_clock::now();
//   d.readData_HDF5(params.fileName);
//   auto stop = std::chrono::steady_clock::now();  
//   std::chrono::duration<double> duration = (stop - start);
//   std::cout << "Reading data took: " << duration.count() << std::endl;
  
//   std::stringstream benchName;
//   benchName << "Benchmark_file_" << d.name << "_eps_" << epsilon_original << "_minPts_" << minPts << "_delta_" << 
//     params.delta << "_memoLimit_"<< params.memoConstraint << "_level_" << params.level << "_shrinkage_" << params.shrinkageFactor << ".txt";
  
//   std::cout << "Saving in file " << benchName.str() << std::endl;

//   if (metric == angular)
//     {
//       d.normalizeData();      //d.meanRemoveData();
//     }

//   std::cerr << "Memory constraint is " << params.memoConstraint << std::endl;
//   SRR_LSHDBSCAN SRR_dbscan(&d,
//         params.delta,
//         params.memoConstraint,
// 				BENCHMARK,
//         benchName.str(),      
//         params.numberOfThreads,
//         params.level,
//         params.shrinkageFactor);
//   SRR_dbscan.introduceMe();
  
//   std::cout << "Running SSR_LSH clustering" << std::endl;
//   SRR_dbscan.performClustering();
//   std::cout << "Writing Labels to file" << std::endl;
//   std::stringstream ss;
//   ss << "Labels_file_" << d.name << "_eps_" << epsilon_original << "_minPts_" << minPts << "_delta_" << params.delta << "_memoLimit_"<< 
//     params.memoConstraint << "_level_" << params.level << "_shrinkage_" << params.shrinkageFactor << ".h5"; 
  
//   SRR_dbscan.writeHDF5(ss.str(), counters);

//   std::cout << counters.stats["total"] << std::endl;

//   return 0;
// }

class SRR {
    std::unique_ptr<SRR_LSHDBSCAN> dbscan;
    dataset* ds;

public:
    
    std::vector<int> fit_predict(
        std::vector<std::vector<float>>& data,
        double delta,
        double mem_constraint,
        bool benchmark,
        std::string benchName,
        size_t numberOfThreads = 2,
        int level=-1,
        double shrinkageFactor=1.0,
        double epsilon_=1.0,
        int minPts_=100
    ) {
        std::vector<int> labels;
        ds = new dataset();
        ds->readData(data);

        epsilon_original = epsilon_;
        epsilon = epsilon_ * epsilon_;
        minPts = minPts_;

        dbscan = std::make_unique<SRR_LSHDBSCAN>(ds, delta, mem_constraint,
            benchmark, benchName, numberOfThreads, level, shrinkageFactor);


        dbscan->performClustering();

        for(auto p: ds->points){
            //std::cout << p.print() << std::endl;
            labels.push_back(p.print());
        }


        return labels;

    }
};

PYBIND11_MODULE(dbscan_srr, m) {
    py::class_<SRR>(m, "SRR")
        .def(py::init<>())
        .def("fit_predict", &SRR::fit_predict);
}


} // namespace python
} // namespace puffinn


