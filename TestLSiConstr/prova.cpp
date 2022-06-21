
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
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
using namespace std;

using namespace std::chrono;

#define maxn 30100
#define MAXITERS 3000
#define K 200
#define EPOCH 10
#define L 750

int vertex_num, num_edges;

bool deleted[maxn], visited[maxn];

vector<vector<int> > adjlist;

vector<int> nodes_degree;

vector<int> V;

int vertex_num_visited;

int vertex_list[maxn];

Conf C;

int weight[maxn];

bool articulation[maxn];
int artcost[maxn];

vector<int> S, bestS;

int SolCost, BestSolCost;

void compBasedNeighbSearch(){

    int bestSolCost, SolCost;
    int iter = 0;
    int count = 0;
    int swaps = 0;

    fill(weight, weight+vertex_num, 0);
   
    SolCost = computeCost();    
    bestSolCost = SolCost;
    bestS = S;

    while(iter < MAXITERS){
        int cc = randomLargeCC();
      
        //int dv = selectNodeWithHigherWeigth(cc);
        int dv, bestDecr;
        biggesDecrement(C.connectedComponents[cc][0], dv, bestDecr);

        deleteVertex(dv);

        S.push_back(dv);

        int iav = wnodei(S);

        restoreVertex(S[iav]);

        swap(S[iav], S[S.size()-1]);
        S.pop_back();
        

        SolCost = computeCost();
  
        if(SolCost < bestSolCost){
            count = 0;
            //we update the best solution
            bestSolCost = SolCost;
            bestS = S;
            weight[dv] = 0;
            iter = 0;
          
        } else {
            
            iter++;
            count++;
            
        }

        if(count == L){   
           for(int i = 0; i < int(0.2*S.size()); ++i){
                int k = rand()%S.size();
                int v = S[k];
                swap(S[k], S[S.size()-1]);
                S.pop_back();
                int l = rand()%vertex_num;
                while(deleted[l]) l = rand()%vertex_num;
                S.push_back(l);
                update(v, l);
            }
            count = 0;
            SolCost = computeCost();   
        }
    }

    for(auto &v: S){
        deleted[v] = false;
    }
    for(auto &v: bestS){
        deleted[v] = true;
    }

    bestSwap();
}

void shuffle(vector<int> &vec){
    int l = vec.size();
    int temp = 0;
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

}


void ComponentBased(){
    feasibleSolution();
    computeConf();
    compBasedNeighbSearch();
}

int main(){

    string filename, file;

    cout<<"specify network"<<endl;
    cin>>file;

    filename = "../cnd/" + file + ".txt";
    
    read_graph(filename);

    srand(time(0));

    int costs = 0;
    double exectime = 0;

    int costiter;

    int bestCost = vertex_num*(vertex_num-1)/2;

    for(int j = 0; j < EPOCH ; ++j){
        
        auto start = high_resolution_clock::now();
        ComponentBased();
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(stop - start);
        exectime += double(duration.count())/1000;
        costiter = compute_cost_func();
        if(costiter < bestCost) bestCost = costiter;
        costs += costiter;

        cout<<"epoch "<<j<<" done"<<endl;
    }

    cout << "avg Time taken by function: " <<exectime/EPOCH << " seconds" << endl;

    cout<<endl<<"---------------------------------"<<endl;
    cout<<"avg Cost :"<<costs/EPOCH<<endl;
    cout<<"best cost: "<<bestCost<<endl;


}