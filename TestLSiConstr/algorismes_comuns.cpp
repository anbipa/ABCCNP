
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
#include <unordered_map>
#include "algorismes_comuns.h"

using namespace std;

using namespace std::chrono;

#define maxn 30100

extern int vertex_num, num_edges;

extern bool deleted[maxn], visited[maxn];

extern vector<vector<int> > adjlist;

extern vector<int> nodes_degree;

extern vector<int> V;

extern int vertex_num_visited;

extern int vertex_list[maxn];

extern Conf C;

extern int weight[maxn];

extern bool articulation[maxn];
extern int artcost[maxn];

extern vector<int> S, bestS;
extern int SolCost, BestSolCost;


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
            ++num_edges;
        }

        nodes_degree[i] = degree;
        degree = 0;

	}

    num_edges /= 2;
}

int pairCon(int a){
    return a*(a-1)/2;
}


void DepthFirst(int v, int &ans){
    visited[v] = true;
    ++ans;
    for(int i = 0; i < adjlist[v].size(); ++i) {
        int neibor = adjlist[v][i];
        if(!deleted[neibor] && !visited[neibor]){
            DepthFirst(neibor,ans);
        }
    }

}

void DepthFirstComponents(int v, int c, vector<int> &vec, int &ans) {
	visited[v] = true;
    vec[v] = c;
    ++ans;
    for(int i = 0; i < adjlist[v].size(); ++i) {
        int neibor = adjlist[v][i];
        if(!deleted[neibor] && !visited[neibor]){
            DepthFirstComponents(neibor, c, vec, ans);
        }
    }
}

void DepthFirstComponents2(int v, vector<int> &vec, int &ans) {
	visited[v] = true;
    vec.push_back(v);
    ++ans;
    for(int i = 0; i < adjlist[v].size(); ++i) {
        int neibor = adjlist[v][i];
        if(!deleted[neibor] && !visited[neibor]){
            DepthFirstComponents2(neibor, vec, ans);
        }
    }
}

vector<vector<int> > obtainListConnectedComponents(){
    fill(visited, visited + vertex_num, false);
    vector<vector<int> > list;

    int size = 0;

    for(int i = 0; i < vertex_num; ++i){
        if(!visited[i] and !deleted[i]){
            vector<int> elements;
            DepthFirstComponents2(i,elements, size);
            list.push_back(elements);
        }
    }

    return list;

}

Configuration computeConfiguration(int &cost){
    Configuration out;
    fill(visited, visited + vertex_num, false);
    vector<int > list(vertex_num,-1);
    unordered_map<int,int> sizes;

    cost = 0;


    for(int i = 0; i < vertex_num; ++i){
        if(!visited[i] and !deleted[i]){
            int size = 0;
            DepthFirstComponents(i,i,list,size);
            sizes.insert(make_pair(i,size));
            cost += size*(size-1)/2;
        } else if (deleted[i]){
            list[i] = i;
        }
    }

    out.indices = list;
    out.sizes = sizes;
    return out;

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

//compute degree of node x
int degree(int x){
    int degree = 0;
    for (int i = 0; i < adjlist[x].size(); i++)
    {
        if(!deleted[adjlist[x][i]]){
            ++degree;
        }
    }
    return degree;
}


int computeObjectiveFunction(vector<int> &sizes){
    int cost = 0;
    for (int i = 0; i < sizes.size(); i++)
    {
        cost += sizes[i]*(sizes[i]-1)/2;
    }
    return cost;
}

int IEM(const Configuration &c, int u){
    Configuration out;
    out.sizes = c.sizes;
    out.indices = c.indices;

    int sum1 = 0;
    int sum2 = 0;

    for (int i = 0; i < adjlist[u].size(); i++)
    {
        int v = adjlist[u][i];
        if(!deleted[v]){
            int vComp = out.indices[v];
            auto it = out.sizes.find(vComp);
            if(it != out.sizes.end()){
                sum1 += it -> second;
                sum2 += it -> second *(it->second-1)/2;
                /*
                for(int j = 0; j < out.indices.size(); ++j){
                    if(out.indices[j] == vComp){
                        out.indices[j] = out.indices[u];
                    }
                }
                */
                out.sizes.erase(it);
            }
        }
    }

    return (sum1+1)*sum1/2 - sum2;
}

int ARN(Configuration &c){

    int worst_node = 0;
    int min_increase = vertex_num*(vertex_num-1)/2;

    for(int i = 0; i < vertex_num; ++i) {
        if(deleted[i]){
            int increment = IEM(c, i);
            if(increment < min_increase){
                min_increase = increment;
                worst_node = i;
            }
        }
    }

    return worst_node;
}

void restoreWorstNodes(int n){
    
    int c;
    Configuration C = computeConfiguration(c);

    for (int i = 0; i < n; i++)
    {
        int worst_node = ARN(C);
        int sizesSize = C.sizes.size();

        C.sizes[sizesSize] = 1;
        C.indices[worst_node] = sizesSize;

        for(int j = 0; j < adjlist[worst_node].size(); ++j){
            int v = adjlist[worst_node][j];
            int vComp = C.indices[v];
            if(!deleted[v]){
                if(vComp != sizesSize){
                    C.sizes[sizesSize] += C.sizes[vComp];
                    for(int k = 0; k < C.indices.size(); ++k){
                        if(C.indices[k] == vComp){
                            C.indices[k] = sizesSize;
                        }
                    }
                }
            }
        }
        deleted[worst_node] = false;
    }

}


// u està a S
// calculem l'increment del cos de moure u de S a V\S
int compute_incr(int u){

    int sum1 = 0;
    int sum2 = 0;
    fill(visited, visited + vertex_num, false);
    visited[u] = true;

    for (int i = 0; i < adjlist[u].size(); i++){
        int v = adjlist[u][i];
        if(!deleted[v] and !visited[v]){
            int size = 0;
            DepthFirst(v,size);
            sum1 += size;
            sum2 += size*(size-1)/2;
        }
    }
    return (sum1+1)*sum1/2 - sum2;
}

// u no està a S
// calculem el decrement del cost de moure u de V\S a S
int compute_decr(int u){

    int sum1 = 0;
    int sum2 = 0;
    fill(visited, visited + vertex_num, false);
    visited[u] = true;

    for (int i = 0; i < adjlist[u].size(); i++){
        int v = adjlist[u][i];
        if(!deleted[v] and !visited[v]){
            int size = 0;
            DepthFirst(v,size);
            sum1 += size;
            sum2 += size*(size-1)/2;
        }
    }
    return (sum1+1)*sum1/2 - sum2;
}

int fastWorstNode(){
    
    int minIncrement = vertex_num*(vertex_num-1)/2;
    int worst_node = 0;

    for (int i = 0; i < vertex_num; i++)
    {
        if(deleted[i]){
            int incr = compute_incr(i);
            if (incr <= minIncrement){
                minIncrement = incr;
                worst_node = i;
            }
        }
    }
    return worst_node;
}

int fastBestNode(const vector<int> &nodes){
    int maxDecrease = 0;
    int best_node = 0;

    for (int i = 0; i < nodes.size(); i++)
    {
        int v = nodes[i];
        if(!deleted[v]){
            int decr = compute_incr(v);
            if (decr >= maxDecrease){
                maxDecrease = decr;
                best_node = v;
            }
        }
    }
    return best_node;
}


int fastCompIncr(Configuration &c, int u){


    int sum1 = 0;
    int sum2 = 0;

    fill(visited, visited + vertex_num, false);
    visited[u] = true;
    vector<bool> visitedComp = vector<bool>(c.sizes.size(), false);

    for (int i = 0; i < adjlist[u].size(); i++){
        int v = adjlist[u][i];
        int cv = c.indices[v];
        
        if(!deleted[v] and !visitedComp[cv]){
            cout<<"cv: "<<cv<<endl;
            int size = c.sizes[cv];
            sum1 += size;
            sum2 += size*(size-1)/2;
            visitedComp[cv] = true;
        }
    }

    return (sum1+1)*sum1/2 - sum2;
}

int newFastWorstNode(Configuration &c){

    int minIncrement = vertex_num*(vertex_num-1)/2;
    int worst_node = 0;

    for (int i = 0; i < vertex_num; i++)
    {
        if(deleted[i]){
            
            
            int incr = fastCompIncr(c,i);
            
            if (incr <= minIncrement){
                minIncrement = incr;
                worst_node = i;
          
            }
        }
    }
    
    return worst_node;

}

void updateConfiguration(Configuration &c, int u, int v){
    cout<<"entra update"<<endl;

    c.sizes.insert(make_pair(u, 1));
    c.indices[u] = u;

    for(int j = 0; j < adjlist[u].size(); ++j){
        int n = adjlist[u][j];
        int nComp = c.indices[n];
        if(!deleted[n]){
            auto itn = c.sizes.find(nComp);
            if(itn != c.sizes.end()){
                c.sizes[u] += itn->second;
                for(int k = 0; k < vertex_num; ++k){
                    if(c.indices[k] == nComp){
                        c.indices[k] = u;
                    }
                }
                c.sizes.erase(itn);
            }
        }
    }

    int vComp = c.indices[v];
    fill(visited, visited + vertex_num, false);
    visited[v] = true;

    for(int j = 0; j < adjlist[v].size(); ++j){
        int n = adjlist[v][j];
        int nComp = c.indices[n];
        if(!deleted[n] and nComp == vComp){
            int size = 0;
            vector<int> vec = c.indices;
            DepthFirstComponents(n,n,vec,size);
            c.indices = vec;
            c.sizes.insert(make_pair(n, size));
        }
    }

    c.sizes.erase(vComp);
    cout<<"surt update"<<endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------

void dfs(int u){
    
    vertex_list[vertex_num_visited] = u;
    vertex_num_visited++;
    visited[u] = true;
    for (int i = 0; i < adjlist[u].size(); i++){
        int v = adjlist[u][i];
        if(!deleted[v] && !visited[v]){
            dfs(v);
        }
    }
}


void externalDfs(int n){
    vertex_num_visited = 0;
    dfs(n);

    int newcc = C.freeid.top();
    C.freeid.pop();
    C.usedid.insert(newcc);

    vector<int>& vecc = C.connectedComponents[newcc];
    vecc.resize(vertex_num_visited);

    C.sizes[newcc] = vertex_num_visited;

    for(int i = 0; i < vertex_num_visited; ++i){
        vecc[i] = vertex_list[i];
        C.indices[vertex_list[i]] = newcc;
    }
}

int computeConf(){
    fill(visited, visited + vertex_num, false);

    C.indices = vector<int>(vertex_num, -1);
    C.sizes = vector<int>(vertex_num, 0);
    C.connectedComponents = vector<vector<int> >(vertex_num, vector<int>());

    while(!C.freeid.empty()) C.freeid.pop();
    for (int i = 0; i < vertex_num; i++) C.freeid.push(i);
    C.usedid.clear();

    int cost = 0;
    for (int i = 0; i < vertex_num; i++){
        if(!deleted[i] and !visited[i]){
            externalDfs(i);
            cost += (vertex_num_visited-1)*vertex_num_visited/2;
        }
    }
    fill(visited, visited + vertex_num, false);
    return cost;

}

void deleteVertex(int v){
    int vcc = C.indices[v];

    SolCost -= pairCon(C.sizes[vcc]);

    deleted[v] = true;
    fill(visited, visited + vertex_num, false);
    visited[v] = true;


    for(int i = 0; i < adjlist[v].size(); ++i){
        int n = adjlist[v][i];
        int ncc = C.indices[n];
        if(!deleted[n] and ncc == vcc){
            externalDfs(n);
            SolCost += pairCon(vertex_num_visited);
        }
    }
    C.freeid.push(vcc);
    C.usedid.erase(vcc);
}


void restoreVertex(int v){
    deleted[v] = false;
    fill(visited, visited + vertex_num, false);
    int n, ncc;

    int newcc = C.freeid.top();
    C.freeid.pop();
    C.usedid.insert(newcc);

    vector<int>& vecc = C.connectedComponents[newcc];
    vector<int> list;
    list.push_back(v);
    int newsize = 1;

    for (int i = 0; i < adjlist[v].size(); i++){
        n = adjlist[v][i];
        ncc = C.indices[n];
        if(!deleted[n] && !visited[ncc]){
            list.insert(list.end(), C.connectedComponents[ncc].begin(), C.connectedComponents[ncc].end());
            newsize += C.sizes[ncc];
            SolCost -= pairCon(C.sizes[ncc]);
            visited[ncc] = true;
            C.freeid.push(ncc);
            C.usedid.erase(ncc);            
        }
    }

    C.sizes[newcc] = newsize;
    vecc = list;
    for(auto &elem: vecc){
        C.indices[elem] = newcc;
    }

    //update solCost
    SolCost += pairCon(newsize);

}

int compIncr(int u){
    int sum1 = 0;
    int sum2 = 0;

    fill(visited, visited + vertex_num, false);
    for(int i = 0; i < adjlist[u].size(); ++i){
        int v = adjlist[u][i];
        int cv = C.indices[v];
        if(!deleted[v] and !visited[cv]){
            int size = C.sizes[cv];
            sum1 += size;
            sum2 += size*(size-1)/2;
            visited[cv] = true;
        }
    }

    return (sum1+1)*sum1/2 - sum2;
}

int compDecr(int u){
    int sum1 = 0;
    int sum2 = 0;

    fill(visited, visited + vertex_num, false);
    visited[u] = true;
    for(int i = 0; i < adjlist[u].size(); ++i){
        int v = adjlist[u][i];
        if(!deleted[v] and !visited[v]){
            vertex_num_visited = 0;
            dfs(v);
            sum1 += vertex_num_visited;
            sum2 += vertex_num_visited*(vertex_num_visited-1)/2;        
        }
    }

    return (sum1+1)*sum1/2 - sum2;
}

int wnode(){
    int minIncrement = vertex_num*(vertex_num-1)/2+1;
    vector<int> wnodes;

    for (int i = 0; i < vertex_num; i++)
    {
        if(deleted[i]){
            int incr = compIncr(i);
            if (incr < minIncrement){
                minIncrement = incr;
                wnodes.clear();
                wnodes.push_back(i);
            } else if (incr == minIncrement){
                wnodes.push_back(i);
            }
        }
    }
    return wnodes[rand()%wnodes.size()];
}


int probWnodei(const vector<int> &S, double r){

    vector<double> incr;
    vector<int> wnodes;
    int minIncr = vertex_num*(vertex_num-1)/2+1;
    int maxIncr = -1;
    int totIncr = 0;
    for(int i = 0; i < S.size(); ++i){
        int u = S[i];
        int incr_u = compIncr(u);
        incr.push_back(incr_u);
        if (incr_u < minIncr){
            minIncr = incr_u;
            wnodes.clear();
            wnodes.push_back(i);
        } else if(incr_u == minIncr){
            wnodes.push_back(i);
        }
        if(incr_u > maxIncr){
            maxIncr = incr_u;
        }
    }

    maxIncr++;

    for(int i = 0; i < incr.size(); ++i){
        incr[i] = maxIncr - incr[i];
        totIncr += incr[i];
    }


    //ramdom number between 0 and 1
    double l = (double)rand() / (RAND_MAX);
    if (l < r){
        return wnodes[rand()%wnodes.size()];
    } else {
        //random number between 0 and totIncr
       
        int num = rand()%(totIncr);
        
        for(int i = 0; i < incr.size(); ++i){
            if(num < (incr[i])){
                return i;
            } else {
                num -= incr[i];
            }
        }
    }

    cout<<"lool"<<endl;
}

int wnodei(const vector<int> &S){
    int minIncrement = vertex_num*(vertex_num-1)/2+1;
    vector<int> wnodes;

    for (int i = 0; i < S.size(); i++)
    {
        int incr = compIncr(S[i]);
        if (incr < minIncrement){
            minIncrement = incr;
            wnodes.clear();
            wnodes.push_back(i);
        } else if(incr == minIncrement){
            wnodes.push_back(i);
        }
    }
    
    return wnodes[rand()%wnodes.size()];
}


int randomLargeCC(){
    int maxCC = 0;
    int minCC = vertex_num;
    int size;
    vector<int> Lcc;
    for (auto &cc : C.usedid){
        size = C.sizes[cc];
        if(size > maxCC) maxCC = size;
        if(size < minCC) minCC = size;
    }
    
    int L = (maxCC + minCC)/2;

    for (auto &cc : C.usedid){
        size = C.sizes[cc];
        if(size >= L){
            Lcc.push_back(cc);
        }
    }
    
    int r = rand() % Lcc.size();
    return Lcc[r];

}

int computeCost(){
    int cost = 0;
    for (auto &cc : C.usedid){
        cost += (C.sizes[cc]-1)*C.sizes[cc]/2;
    }
    return cost;
}

void update(int u, int v){
    deleted[v] = true;
    restoreVertex(u);
    deleteVertex(v);
    SolCost = computeCost();
    /*
    if(SolCost != compute_cost_func()){
        cout << "error" << endl;
    }
    */
}

pair<int,int> biggestDecr(int v){
    int maxDecr = 0;
    int best_node = -1;

    int ccv = C.indices[v];

    for(int k = 0; k < C.connectedComponents[ccv].size(); ++k){
        int n = C.connectedComponents[ccv][k];
        if(!deleted[n]){
            int decr = compDecr(n);
            if(decr > maxDecr){
                maxDecr = decr;
                best_node = n;
            }
        }
    }
         
    return make_pair(best_node, maxDecr);
    
}

//From IEM 
bool bestSwap(){

    vector<pair<int,int> > decrements;

    //obtenim els punts d'articulacio amb cost O(V + E)
    getArticulationPonints();

    for(int i = 0; i < vertex_num; i++){
        if(!deleted[i]){
            decrements.push_back(make_pair(i, artcost[i]));
        }
    }

    pair<int, int> bSwap;
    int bestDecr = 1;
    int bestIncr = 1;
    int best, bigDecr;
    Conf Caux = C;

    for(int i = 0; i < vertex_num; i++){
        if(deleted[i]){
            int incr = compIncr(i);
            restoreVertex(i);
            biggesDecrement(i, best, bigDecr);
            vector<pair<int,int> > Vset;

            for(auto &elem: decrements){
                if(elem.second > bigDecr && C.indices[elem.first] != C.indices[i]){
                    Vset.push_back(elem);
                }
            }

            int w;
            int maxDecr = -1;

            if(Vset.size() == 0){
                w = best;
                maxDecr = bigDecr;
            } else {
                for(auto &elem: Vset){
                    if(elem.second > maxDecr){
                        maxDecr = elem.second;
                        w = elem.first;
                    }
                }
            }
            deleted[i] = true;
            C = Caux;

            if(bestDecr - bestIncr < maxDecr - incr){
                bSwap = make_pair(i, w);
                bestDecr = maxDecr;
                bestIncr = incr;
            } 
        }
    }

    if (bestDecr - bestIncr > 0){
        for(int i = 0; i < S.size(); i++){
            if(S[i] == bSwap.first){
                swap(S[i], S[S.size()-1]);
                S.pop_back();
            }
        }
       
        S.push_back(bSwap.second);
        update(bSwap.first, bSwap.second);
        return true;
    } else {
        return false;
    }
}

// first is the articulation point, second the pairwise connectivity
void articulationPoints(int start, int d, vector<int> &low, vector<int> &depth, vector<int> &parent, int &treesize, int ccsize){
    visited[start] = true;
    low[start] = depth[start] = d;
    int size1 = 0;
    int size2 = 0;
    int children = 0;
    bool ap = false;
    for(int i = 0; i < adjlist[start].size(); ++i){
        int v = adjlist[start][i];
        if(!deleted[v] && !visited[v]){
            parent[v] = start;
            children++;
            int ts = 0;
            articulationPoints(v, d+1, low, depth, parent, ts, ccsize);
            treesize += ts;
            if(low[v] >= depth[start]){
                size1 += ts*(ts-1)/2;
                size2 += ts;
            }
            low[start] = min(low[start], low[v]);
            if(parent[start] == -1 && children > 1) ap = true;
            if(parent[start] != -1 && low[v] >= depth[start]) ap = true;
        } else if(!deleted[v] && v != parent[start]){
            low[start] = min(low[start], depth[v]);
        }
    }

    if(ap){
        articulation[start] = true;
        artcost[start] = pairCon(ccsize) - (size1 + pairCon(ccsize-size2-1));
    } else {
        artcost[start] = ccsize-1;
    }

    treesize++;
}

void getArticulationPonints(){
    fill(visited, visited + vertex_num, false);
    vector<int> low(vertex_num, vertex_num+1);
    vector<int> depth(vertex_num, 0);
    vector<int> parent(vertex_num, -1);
    fill(articulation, articulation + vertex_num, false);
    fill(artcost, artcost + vertex_num, 0);
    for(int i = 0; i < vertex_num; ++i){
        if(!deleted[i] && !visited[i]){
            int d = 0;
            int treesize = 0;
            articulationPoints(i, d, low, depth, parent, treesize, C.sizes[C.indices[i]]);
        }
    }
}

// donat un node pertanyent a S u, obtenim el node v de la mateixa component connexa de u tal que el decrement de afegir a S és màxim
void biggesDecrement(int u, int &best, int &bestDecr){
    fill(visited, visited + vertex_num, false);
    vector<int> low(vertex_num, vertex_num+1);
    vector<int> depth(vertex_num, 0);
    vector<int> parent(vertex_num, -1);
    fill(articulation, articulation + vertex_num, false);
    fill(artcost, artcost + vertex_num, 0);

    int cc = C.indices[u];
    int sizeCC = C.sizes[cc];

    int d = 0;
    int treesize = 0;
    articulationPoints(u, d, low, depth, parent, treesize, sizeCC);

    bestDecr = -1;
    best = u;

    vector<int> ties;

    for(int i = 0; i < vertex_num; ++i){
        if(!deleted[i] && artcost[i] > bestDecr){
            ties.clear();
            bestDecr = artcost[i];
            best = i;
            ties.push_back(i);
        } else if(!deleted[i] && artcost[i] == bestDecr){
            ties.push_back(i);
        }
    }
    
    // escollim un node random d'entre els millors nodes articulació
    if(ties.size() >= 1) best = ties[rand()%ties.size()];

}

int bnode(){

    int bestDecr = -1;

    getArticulationPonints();

    vector<int> ties;

    for(int i = 0; i < vertex_num; ++i){
        if(!deleted[i] && artcost[i] > bestDecr){
            ties.clear();
            bestDecr = artcost[i];
            ties.push_back(i);
        } else if(!deleted[i] && artcost[i] == bestDecr){
            ties.push_back(i);
        }
    }
    
    return ties[rand()%ties.size()];
    
}
