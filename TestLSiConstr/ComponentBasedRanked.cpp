
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
#define ELITESETSIZE 10
#define K 265
#define EPOCH 10
#define I 1000

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


int minDegree(const vector<bool> &labeled, const vector<int> &degree){
    int min = vertex_num+1;
    int mini = 0;
    for(int i = 0; i < vertex_num; ++i){
        if(!labeled[i] and degree[i] < min){
            min = degree[i];
            mini = i;
        }
    }

    return mini;
}

vector<int> MaxIndepSet(){
    vector<bool> labeled(vertex_num, false);
    int nlabeled = 0;
    vector<int> MIS;
    vector<int> degree = nodes_degree;

    while(nlabeled != vertex_num){
        int u = minDegree(labeled, degree);
        labeled[u] = true;

        MIS.push_back(u);

        ++nlabeled;

        for(auto& elem: adjlist[u]){
            if(!labeled[elem]) ++nlabeled;
            labeled[elem] = true;
            for(auto& n: adjlist[elem]){
                degree[n]--;
            }
        }
    }

    return MIS;

}

// que passa amb el node 0?? (millorar greedy MIS) 
// que passa quan K > |MIS|??
void MISGreedy(){

    vector<int> MIS = MaxIndepSet();
    fill(deleted, deleted + vertex_num, true);
    
    for(auto& elem: MIS){
        deleted[elem] = false;
    }

    int cost = computeConf();

    while(MIS.size() != vertex_num - K){
        int i = wnode();
        MIS.push_back(i);
        restoreVertex(i);
    }
}


int selectNodeWithHigherWeigth(int cc){

    int max = -1;
    int idx = 0;
    vector<int>& list = C.connectedComponents[cc];
    for (int i = 0; i < list.size(); i++)
    {
        if(weight[list[i]] > max){
            idx = list[i];
            max = weight[list[i]];
        } else if (weight[list[i]] == max and degree(list[i]) > degree(idx)){
            idx = list[i];
        }
        weight[list[i]]++;
    }
    
    weight[idx]--;
    return idx;

}



void compBasedNeighbSearch(){

    int bestSolCost, SolCost;
    int iter = 0;
    int count = 0;

    bool bestSol[maxn];
    copy(deleted, deleted + vertex_num, bestSol);

    fill(weight, weight+vertex_num, 0);
   
    SolCost = computeCost();    
    //SolCost = computeConf();
    bestSolCost = SolCost;
    
    while(iter < MAXITERS){

        int cc = randomLargeCC();
      
        int dv = selectNodeWithHigherWeigth(cc);

        deleteVertex(dv);

        int av = wnode();

        restoreVertex(av);

        SolCost = computeCost();

        if (SolCost != compute_cost_func()){
            cout << "ERROR" << endl;
            exit(0);
        }
  
        if(SolCost < bestSolCost){
            //cout<<"improvement"<<endl;
            iter = 0;
            count = 0;
            //we update the best solution
            bestSolCost = SolCost;
            copy(deleted, deleted + vertex_num, bestSol);
            weight[dv] = 0;
          
        } else {
            iter++;
            count++;
        }


        // swap 20% of nodes if no improvement
        if(count == I){
            for(int i = 0; i < int(0.2*vertex_num); ++i){
                
                int k = rand()%vertex_num;
                int l = rand()%vertex_num;
                swap(deleted[k], deleted[l]);
            }
            cout<<"swap 20%"<<endl;

            SolCost = computeConf();
          
            count = 0;
        }
    }

    copy(bestSol, bestSol + vertex_num, deleted);
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

    for(int i = 0; i < K; ++i) deleted[shuffled[i]] = true;

}


void ComponentBased(){

    MISGreedy();
    //feasibleSolution();
    compBasedNeighbSearch();
}

int main(){


    string filename, file;

    cout<<"specify network"<<endl;
    cin>>file;

    filename = "../cnd/" + file + ".txt";
    
    read_graph(filename);
    //construction();
    
    //show_adjlist();
    //show_critical_nodes();

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