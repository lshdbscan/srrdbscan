#pragma once
#include <vector>
#include <algorithm>

class UF
{
private:
    std::vector<int> *ids;
    std::vector<int> *weights;
    
public:
    UF(int n);
    ~UF();
    
    void mergePoint(int i, int j);
    bool sameCluster(int i, int j);
    bool isClustered(int i);
    int find(int i);

};

