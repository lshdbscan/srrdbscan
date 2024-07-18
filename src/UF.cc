#include <UF.h>

UF::UF(int size){
    ids = new std::vector<int>(size+1);
    weights = new std::vector<int>(size+1,1);
    int x = 0;
    std::generate(ids->begin(), ids->end(), [&]{return x++;});
}

UF::~UF(){
    delete ids;
    delete weights;
}

int UF::find(int i){
    while((*ids)[i] != i){
        i = (*ids)[i];
    }
    return i;
}

bool UF::sameCluster(int i, int j){
    return(find(i) == find(j));    
}

bool UF::isClustered(int i){
    return (*ids)[find(i)] != i;
}


void UF::mergePoint(int i, int j){ //prefers to merge j into i, such that i is the parrent when weight[i] >= weight[j]
    int ii = find(i);
    int jj = find(j);
    if(ii == jj) return;
    if((*weights)[ii] < (*weights)[jj]){
        (*ids)[ii] = jj;        
        (*weights)[jj] += (*weights)[ii];
    }
    else{
        (*ids)[jj] = ii;
        (*weights)[ii] += (*weights)[jj];
    } 
    return;
}

