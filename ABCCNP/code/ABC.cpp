
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <queue>
#include <limits>
#include <chrono>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "algorismes_comuns.h"
#include "Timer.h"

using namespace std;

using namespace std::chrono;


#define maxn 30100


//parameters
int candSize = 20;
double prob = 0.3;
int K = 150;
string input_file = "WS1000";
int SN = 10;
int MaxTrial = 20;
double deltaK = 0.5;
int lIters = 30;
double pwn = 0.7;
double tournSelprob = 0.5;
double time_limit = 3600;
int KBV = 37;
bool tuning_version = false;
double distrate = 0.5;


//actual solució sent tractada
vector<int> S, bestS;

//cost de la solució sent tractada
int SolCost, BestSolCost;

int vertex_num, num_edges;

bool deleted[maxn], visited[maxn];

vector<vector<int> > adjlist;

vector<int> nodes_degree;

vector<int> V;

int vertex_num_visited;

int vertex_list[maxn];

Conf C, bestC;

int weight[maxn];

bool articulation[maxn];
int artcost[maxn];

Timer timer;
double t;



void read_parameters(int argc, char **argv) {

    int iarg=1;
    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) input_file = argv[++iarg];
        else if (strcmp(argv[iarg],"-lIters")==0) lIters = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-cs")==0) candSize = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-p")==0) prob = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-SN")==0) SN = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-MaxTrial")==0) MaxTrial = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-deltaK")==0) deltaK = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pwn")==0) pwn = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-tsp")==0) tournSelprob = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-tuning")==0) tuning_version = true;
        else if (strcmp(argv[iarg],"-distrate")==0) distrate = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-h")==0) {
            cout << "Usage: " << argv[0] << " [-i input_file] [-lIters lIters] [-cs candSize] [-p prob] [-SN SN] [-MaxIter MaxIter] [-MaxTrial MaxTrial] [-deltaK deltaK] [-pwn pwn] [-tsp tournSelprob]" << endl;
            exit(0);
        }
        else {
            cout << "Unknown option " << argv[iarg] << endl;
            exit(0);
        }
        iarg++;
    }
}

void setK(){
    string file = input_file.substr(input_file.find_last_of("/\\") + 1);
    if(file == "BA500.txt") K = 50, KBV = 195;
    else if(file == "BA1000.txt") K = 75, KBV = 558;
    else if(file == "BA2500.txt") K = 100, KBV = 3704;
    else if(file == "BA5000.txt") K = 150, KBV = 10196;
    else if(file == "ER235.txt") K = 50, KBV = 295;
    else if(file == "ER466.txt") K = 80, KBV = -1;
    else if(file == "ER941.txt") K = 140, KBV = -1;
    else if(file == "ER2344.txt") K = 200, KBV = -1;
    else if(file == "FF250.txt") K = 50, KBV = 194;
    else if(file == "FF500.txt") K = 110, KBV = 257;
    else if(file == "FF1000.txt") K = 150, KBV = 1260;
    else if(file == "FF2000.txt") K = 200, KBV = 4545;
    else if(file == "WS250.txt") K = 70, KBV = -1;
    else if(file == "WS500.txt") K = 125, KBV = -1;
    else if(file == "WS1000.txt") K = 200, KBV = -1;
    else if(file == "WS1500.txt") K = 265, KBV = -1;

}


void shuffle(vector<int> &vec){
    int l = vec.size();
    int randomIdx = 0;
    for(int i = 0; i < l; ++i){
        randomIdx = rand()%l;
        swap(vec[i], vec[randomIdx]);
    }

}

void feasibleSolution(){
    vector<int> shuffled;
    for(int i = 0; i < vertex_num; ++i){
        shuffled.push_back(i);
    }

    shuffle(shuffled);

    fill(deleted, deleted + vertex_num, false);
    S.clear();

    for(int i = 0; i < K; ++i) {
        deleted[shuffled[i]] = true;
        S.push_back(shuffled[i]);
    }

    SolCost = computeConf();

}

int compDifference(vector<int> &S2) {
	int sum = 0;
	for (auto &num:S2) sum += deleted[num];
	return sum;
}





int nonDetMaxDegree(const vector<int> &ndegree, const vector<bool> &cover){

    //random number between 0 and 1
    double r = (double)rand() / (RAND_MAX);

    int returnNode = -1;

    if (r < prob){

        int maxdeg = -1;
        int maxnode = 0;

        for(int i = 0; i < vertex_num; ++i){
            if(!cover[i] and ndegree[i] > maxdeg){
                maxdeg = ndegree[i];
                maxnode = i;
            }
        }

        returnNode = maxnode;

    } else {

        vector<int> maxdeg, maxnode;

        for(int i = 0; i < vertex_num; ++i){
            if(!cover[i]){
                if(maxdeg.size() < candSize){
                    maxdeg.push_back(ndegree[i]);
                    maxnode.push_back(i);
                    break;
                } else {
                    for(int k = 0; k < candSize; ++k){
                        if(maxdeg[k] < ndegree[i]){
                            maxdeg[k] = ndegree[i];
                            maxnode[k] = i;
                            break;
                        }
                    }
                }
            }
        }

        int d = rand() % maxnode.size();
        returnNode = maxnode[d];
    }

    return returnNode;
}


vector<bool> vertex_cover_sorted(){

    vector<bool> cover(vertex_num, false);

    int nedg = num_edges;
    vector<int> ndegree = nodes_degree;


    while(nedg > 0){

        int v = nonDetMaxDegree(ndegree, cover);
        cover[v] = true;

        for (int i = 0; i < adjlist[v].size(); i++)
        {
            if(!cover[adjlist[v][i]]) {
                nedg--;
            }
            ndegree[adjlist[v][i]]--;
        }
    }

    return cover;
}


void VCoverOrdered(){

    vector<bool> vcover = vertex_cover_sorted();

    int n = 0;

    fill(deleted, deleted + vertex_num, false);
    S.clear();

    for(int i = 0; i < vertex_num; ++i){
        if (vcover[i]){
            ++n;
            deleted[i] = true;
            S.push_back(i);
        }
    }

    SolCost = computeConf();

    while (n > K)
    {
        int i = probWnodei(S, pwn);
        int w = S[i];
        S.erase(S.begin() + i);
        restoreVertex(w);
        --n;
    }
}


void greedy3d(){
    int n = 0;
    BestSolCost = vertex_num*(vertex_num-1)/2;
    bestS = S;
    do
    {
        while(S.size() > K - int(deltaK*K)){
            int iu = probWnodei(S, pwn);
            int u = S[iu];
            restoreVertex(u);
            S.erase(S.begin()+iu);
            if(S.size() == K){
                if(SolCost <= BestSolCost){
                    BestSolCost = SolCost;
                    bestS = S;
                    bestC = C;
                }
            }
        }
        ++n;
        while(S.size() < K + int(deltaK*K)){
            int cc = randomLargeCC();
            int u, bestDecr;
            biggesDecrement(C.connectedComponents[cc][0], u, bestDecr);
            S.push_back(u);
            deleteVertex(u);
            if(S.size() == K){
                if(SolCost <= BestSolCost){
                    BestSolCost = SolCost;
                    bestS = S;
                    bestC = C;
                }
            }
        }
        ++n;
    } while (n<=lIters);
}


int ArtificialBeeColony(){

    vector<int> food_sources[SN];
    Conf configurations[SN];
    int sources_costs[SN];
    int trials[SN];

    for(int i = 0; i < SN; ++i){
        VCoverOrdered();
        while(true){
            bool trobat = true;
            for(int j = 0; j < i; ++j){
                if(compDifference(food_sources[j]) > K*distrate){
                    trobat = false;
                    break;
                }
            }

            if(trobat) break;
            int change = rand() % K;
            swap(S[change], S.back());
            restoreVertex(S.back());
            change = rand() % vertex_num;
            while(deleted[change]) change = rand() % vertex_num;
            S.back() = change;
            deleteVertex(change);

        }
        food_sources[i] = S;
        configurations[i] = C;
        sources_costs[i] = SolCost;
        trials[i] = 0;
    }

    int BestAbcCost = vertex_num*(vertex_num-1);
    vector<int> BestAbcS;

    for(int i = 0; i < SN; ++i){
        if(sources_costs[i] < BestAbcCost){
            BestAbcCost = sources_costs[i];
            BestAbcS = food_sources[i];
        }
    }

    t = timer.elapsed_time(Timer::VIRTUAL);
    while(t < time_limit and BestAbcCost > KBV){
        //**** EMPLOYED BEE PHASE *******
        for(int i = 0; i < SN; ++i){
            C = configurations[i];
            S = food_sources[i];
            SolCost = sources_costs[i];
            fill(deleted, deleted + vertex_num, false);
            for(auto &v: S){
                deleted[v] = true;
            }
            greedy3d();
            if(BestSolCost >= sources_costs[i]){
                trials[i]++;
            } else {
                trials[i] = 0;
                configurations[i] = bestC;
                food_sources[i] = bestS;
                sources_costs[i] = BestSolCost;
            }
            if(BestSolCost < BestAbcCost){
                BestAbcCost = BestSolCost;
                BestAbcS = bestS;
                if(!tuning_version) cout<<"new best cost: "<<BestAbcCost<<endl;
            }
        }
        //*****ONLOOKER BEE PHASE******

        for(int i = 0; i < SN; ++i){
            int source1 = rand() % SN;
            int source2 = rand() % SN;
            while(source1 == source2){
                source2 = rand() % SN;
            }

            int best_source, worst_source, source;
            if(sources_costs[source1] < sources_costs[source2]){
                best_source = source1;
                worst_source = source2;
            } else {
                best_source = source2;
                worst_source = source1;
            }

            //random number between 0 and 1
            double r = (double)rand() / RAND_MAX;
            if(r < tournSelprob){
                source = best_source;
            } else {
                source = worst_source;
            }

            S = food_sources[source];
            C = configurations[source];
            SolCost = sources_costs[source];

            fill(deleted, deleted + vertex_num, false);
            for(auto &v: S){
                deleted[v] = true;
            }
            greedy3d();
            if(BestSolCost >= sources_costs[i]){
                trials[i]++;
            } else {
                trials[i] = 0;
                configurations[i] = bestC;
                food_sources[i] = bestS;
                sources_costs[i] = BestSolCost;
            }
            if(BestSolCost < BestAbcCost){
                BestAbcCost = BestSolCost;
                BestAbcS = bestS;
                if(!tuning_version) cout<<"new best cost: "<<BestAbcCost<<endl;
            }
        }

        //**** SCOUT BEE PHASE *****
        for(int i = 0; i < SN; ++i){
            if(trials[i] == MaxTrial){
                VCoverOrdered();
                while(true){
                    bool trobat = true;
                    for(int j = 0; j < i; ++j){
                        if(compDifference(food_sources[j]) > (K*distrate)){
                            trobat = false;
                            break;
                        }
                    }

                    if(trobat) break;
                    int change = rand() % K;
                    swap(S[change], S.back());
                    restoreVertex(S.back());
                    change = rand() % vertex_num;
                    while(deleted[change]) change = rand() % vertex_num;
                    S.back() = change;
                    deleteVertex(change);
                }
                food_sources[i] = S;
                configurations[i] = C;
                sources_costs[i] = SolCost;
                trials[i] = 0;

                if(SolCost < BestAbcCost){
                    BestAbcCost = SolCost;
                    BestAbcS = S;
                    if(!tuning_version) cout<<"new best cost: "<<BestAbcCost<<endl;
                }
            }
        }
        t = timer.elapsed_time(Timer::VIRTUAL);
        //if (!tuning_version) cout<<"time: "<<t<<endl;
    }

    return BestAbcCost;
}


int main( int argc, char **argv ){

    Timer timer;

    read_parameters(argc, argv);

    setK();

    string filename;

    filename = input_file.c_str();

    if (tuning_version) read_graph(filename);
    else read_graph("../instances/"+filename);
    srand(time(NULL));

    int BestABCCost = ArtificialBeeColony();
    if(!tuning_version) cout<<"time: "<<timer.elapsed_time(Timer::VIRTUAL)<<endl;
    cout<<BestABCCost<<endl;
}
