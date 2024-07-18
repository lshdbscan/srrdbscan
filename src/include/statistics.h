#pragma once

#include<unordered_map>
#include<string>

class Statistics {
    public:
        std::unordered_map<std::string, double> stats;

        void add_measurement(std::string name, double val) {
            stats[name] = val;
        }
};

