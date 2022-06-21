
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
using namespace std;

using namespace std::chrono;



#define maxn 30100
#define MAXITERS 100
#define ELITESETSIZE 10
#define K 50

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


//--------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------//

void show_adjlist() {
	for(int i = 0; i < vertex_num; ++i) {
		cout << i << " : ";
		for(int j = 0; j < adjlist[i].size(); ++j) {
			cout << adjlist[i][j] << " ";
		}
		cout << endl;
	}
}



void read_graph(string filename) {
	ifstream in;
	in.open(filename);
	string line;
	getline(in, line);
	istringstream iss(line);
	iss >> vertex_num;
	adjlist.resize(vertex_num);
    nodes_degree.resize(vertex_num);
	int tmp;
	char c;

    int degree = 0;

	for (int i = 0; i < vertex_num; ++i) {
        V.push_back(i);
		getline(in, line);
		istringstream iss(line);
		iss >> tmp >> c;
		while(iss >> tmp){
            adjlist[i].push_back(tmp);
            ++degree;
        }

        nodes_degree[i] = degree;
        degree = 0;

	}
}

int DepthFirst(int v, int &ans){
    visited[v] = true;
    ++ans;
    for(int i = 0; i < adjlist[v].size(); ++i) {
        int neibor = adjlist[v][i];
        if(!deleted[neibor] && !visited[neibor]){
            DepthFirst(neibor,ans);
        }
    }

}

void DepthFirstComponents(int v, vector<int> &vec) {
	visited[v] = true;
    vec.push_back(v);
    for(int i = 0; i < adjlist[v].size(); ++i) {
        int neibor = adjlist[v][i];
        if(!deleted[neibor] && !visited[neibor]){
            DepthFirstComponents(neibor, vec);
        }
    }
}



int compute_cost_func(){
    int connectivity = 0;
	fill(visited, visited + vertex_num, false);
	for(int i = 0; i < vertex_num; ++i) {
		if(!deleted[i] and !visited[i]) {
            int size = 0;
			DepthFirst(i,size);
			connectivity += size * (size - 1) / 2;
		}
	}
	return connectivity;
}

int worst_node(){
    int cost_ini = compute_cost_func();
    int worst_value = vertex_num*(vertex_num-1)/2;
    int worst_node = 0;

    for (int i = 0; i < vertex_num; i++)
    {
        if(deleted[i]){
            deleted[i] = false;
            int cost = compute_cost_func();
            if ((cost-cost_ini) < worst_value){
                worst_value = cost-cost_ini;
                worst_node = i;
            }
            deleted[i] = true;
        }
    }
    
    return worst_node;

}

int best_node(){
    int cost_ini = compute_cost_func();
    int best_value = 0;
    int best_node = 0;

    for (int i = 0; i < vertex_num; i++)
    {
        if(!deleted[i]){
            deleted[i] = true;
            int cost = compute_cost_func();
            if ((cost_ini-cost) > best_value){
                best_value = cost_ini-cost;
                best_node = i;
            }
            deleted[i] = false;
        }
    }
    
    return best_node;

}

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
void greedy1(){

    vector<int> MIS = MaxIndepSet();
    fill(deleted, deleted + vertex_num, true);
    for(auto& elem: MIS){
        deleted[elem] = false;
    }
    while(MIS.size() != vertex_num - K){
        int i = worst_node();
        MIS.push_back(i);
        deleted[i] = false;
    }
}


bool comp(const pair<int,int> &p1, const pair<int,int> &p2){
    return p1.second > p2.second;
}


//greedy vertex cover sorting the nodes by its degree
vector<int> vertex_cover_sorted(){

    vector<int> cover;

    int num_edges = 0;

    vector<vector<bool> > adjmat(vertex_num, vector<bool>(vertex_num, false));
    for(int i = 0; i < vertex_num; ++i){
        for(int j = 0; j < adjlist[i].size(); ++j){
            adjmat[i][adjlist[i][j]] = true;
            adjmat[adjlist[i][j]][i] = true;
            ++num_edges;
        }
    }

    vector<pair<int,int> > pairdegrees;

    for (int i = 0; i < vertex_num; i++)
    {
        pairdegrees.push_back(make_pair(i,nodes_degree[i]));
    }

    sort(pairdegrees.begin(), pairdegrees.end(), &comp);

    int j = 0;
    while(num_edges > 0 and j < vertex_num){

        int v = pairdegrees[j].first;

        cover.push_back(v);

        for (int i = 0; i < vertex_num; i++)
        {
            if(adjmat[v][i]){
                adjmat[v][i] = false;
                adjmat[i][v] = false;
                --num_edges;
            }
        }
        ++j;
    }
    return cover;
}

//vertex cover with node sorting 
void greedy2(){

    vector<int> vcover = vertex_cover_sorted();

    int n = vcover.size();

    fill(deleted, deleted + vertex_num, false);
    for(auto& elem: vcover){
        deleted[elem] = true;
    }

    while (n > K)
    {
        int w = worst_node();
        deleted[w] = false;
        --n;

    }

}

void shuffle(vector<int> &vec){
    int l = vec.size();
    int temp = 0;
    int randomIdx = 0;
    for(int i = 0; i < l; ++i){
        randomIdx = rand()%l;
        temp = vec[i];
        vec[i] = vec[randomIdx];
        vec[randomIdx] = temp;
    }

}


vector<int> vertex_cover_shuffled(){

    vector<int> cover;

    vector<int> shuffled;

    int oneedges = 0;

    for (int i = 0; i < vertex_num; i++)
    {
        if (nodes_degree[i] != 1 and nodes_degree[i] != 0 or K - shuffled.size() >= vertex_num - i) shuffled.push_back(i);
        if (nodes_degree[i] == 1) ++oneedges;
    }
    
    shuffle(shuffled);

    int num_edges = 0;

    vector<vector<bool> > adjmat(vertex_num, vector<bool>(vertex_num, false));
    for(int i = 0; i < vertex_num; ++i){
        for(int j = 0; j < adjlist[i].size(); ++j){

            if(nodes_degree[i] != 0 and nodes_degree[i] != 1){
                adjmat[i][adjlist[i][j]] = true;
                adjmat[adjlist[i][j]][i] = true;
                ++num_edges;
            }
        }
    }

    num_edges -= oneedges;

    int j = 0;
    while(num_edges > 0 and j < shuffled.size()){

        int v = shuffled[j];

        cover.push_back(v);

        for (int i = 0; i < vertex_num; i++)
        {
            if(adjmat[v][i]){
                adjmat[v][i] = false;
                adjmat[i][v] = false;
                --num_edges;
            }
        }
        ++j;
    }


    return cover;
}

void greedy3(){

    vector<int> vcover = vertex_cover_shuffled();

    int n = vcover.size();

    fill(deleted, deleted + vertex_num, false);

    for(auto& elem: vcover){
        deleted[elem] = true;
    }

    while (n > K)
    {
        int w = worst_node();
        deleted[w] = false;
        --n;

    }
}

vector<vector<int> > obtainListConnectedComponents(){
    fill(visited, visited + vertex_num, false);
    vector<vector<int> > list;

    int size;

    for(int i = 0; i < vertex_num; ++i){
        if(!visited[i] and !deleted[i]){
            vector<int> elements;
            DepthFirstComponents(i,elements);
            list.push_back(elements);
        }
    }

    return list;

}



void compBasedNeighbSearch(){
    
    int bestSolCost, SolCost;
    int iter = 0;

    bestSolCost = compute_cost_func();

    map<int,int> weights;
    for (int i = 0; i < vertex_num; i++)
    {
        weights.insert(make_pair(i,0));
    }
    
    while(iter < MAXITERS){

        vector<vector<int> > listcc = obtainListConnectedComponents();

        int maxCC = 0, minCC = vertex_num, size;

        for(int i = 0; i < listcc.size(); ++i){
            size = listcc[i].size();
            if (size < minCC) minCC = size;
            if (size > maxCC) maxCC = size;
        }

        vector<vector<int> > nodesLargeCC;

        for(int i = 0; i < listcc.size(); ++i){
            if(listcc[i].size() > (maxCC + minCC)/2){
                nodesLargeCC.push_back(listcc[i]);
                
                /*
                cout<<"nodesLargeCC"<<endl;
                    for(auto &elem: listcc[i]) cout<< elem<<" ";
                    cout<<"--------------------"<<endl;
                    */
            } 
        }

        vector<int> comp = nodesLargeCC[rand()%(nodesLargeCC.size())];

        int max = -1;
        int idx = 0;
        for (int i = 0; i < comp.size(); i++)
        {
            if(weights[comp[i]] > max){
                idx = comp[i];
                max = weights[comp[i]];
            } else if (weights[comp[i]] == max and nodes_degree[comp[i]] > nodes_degree[idx]){
                idx = comp[i];
            }
        }

        deleted[idx] = true;

        int w = worst_node();

        if(w != idx){

              for (int i = 0; i < comp.size(); i++)
            {
                if(comp[i] != idx) weights[comp[i]]++;
            }

            deleted[w] = false;

            SolCost = compute_cost_func();

            if(SolCost < bestSolCost){
                iter = 0;
                bestSolCost = SolCost;
                weights[w] = 0;
            } else {
                deleted[idx] = false;
                deleted[w] = true; 
                iter++;
            }
        } else {
            deleted[idx] = false;
            iter++;
        }

        /*
        int count = 0;
        for(int i = 0; i < vertex_num; ++i){
            if(deleted[i]) ++count;
        } cout << count<< endl;
        */
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

//
void greedy4(){

    feasibleSolution();
    compBasedNeighbSearch();
}



void executeApproach(int greedy){
    if (greedy == 0) greedy1();
    else if (greedy == 1) greedy2();
    else if (greedy == 2) greedy3();
    else if (greedy == 3) greedy4();
}



int main()
{

    string filename, file;

    cout<<"specify network"<<endl;
    cin>>file;

    filename = "cnd/" + file + ".txt";

    int greedy = -1;
    cout<<"speciy greedy algorithm"<<endl;
    cin>>greedy;



    read_graph(filename);

    //construction();
    
    //show_adjlist();
    //show_critical_nodes();

    srand(time(0));

    bool deterministic = true;

    if(greedy == 3) deterministic = true;


    if (deterministic){

        auto start = high_resolution_clock::now();
        executeApproach(greedy);
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(stop - start);
        cout.precision(5);
        cout << "Time taken by function: " <<float(duration.count())/1000 << " seconds" << endl;

        for(int i = 0; i < vertex_num; ++i){
            if(deleted[i]) cout<<i<<" ";
        }
        cout<<endl<<"---------------------------------"<<endl;
        cout<<"Cost :"<<compute_cost_func()<<endl;

    } else {
        int costs = 0;
        double exectime = 0;

        int costiter;

        int bestCost = vertex_num*(vertex_num-1)/2;

        for(int i = 0; i < 10 ; ++i){
            
            auto start = high_resolution_clock::now();
            executeApproach(greedy);
            auto stop = high_resolution_clock::now();

            auto duration = duration_cast<milliseconds>(stop - start);
            exectime += double(duration.count())/1000;
            costiter = compute_cost_func();
            if(costiter < bestCost) bestCost = costiter;
            costs += costiter;

            cout<<"epoch "<<i<<" done"<<endl;
        }

        cout << "avg Time taken by function: " <<exectime/10 << " seconds" << endl;

        cout<<endl<<"---------------------------------"<<endl;
        cout<<"avg Cost :"<<costs/10<<endl;
        cout<<"best cost: "<<bestCost<<endl;

    }

    return 0;
} 





