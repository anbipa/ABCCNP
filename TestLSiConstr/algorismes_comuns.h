#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

using namespace std;

struct Configuration
{
    unordered_map<int, int> sizes;
    vector<int> indices;
};

struct Conf 
{
    
    vector< vector<int> > connectedComponents;
    vector<int> sizes;
    vector<int> indices;
    stack<int> freeid;
    set<int> usedid;
};

void show_adjlist();
void read_graph(string filename);
int pairCon(int a);
void DepthFirst(int v, int &ans);
void DepthFirstComponents(int v, int c, vector<int> &vec,  int &ans);
void DepthFirstComponents2(int v, vector<int> &vec,  int &ans);
vector<vector<int> > obtainListConnectedComponents();
int compute_cost_func();
int worst_node();
int best_node();
int degree(int x);
int ARN(Configuration &c);
int IEM(const Configuration &c, int u);
int computeObjectiveFunction(vector<int> &sizes);
Configuration computeConfiguration(int &cost);
void restoreWorstNodes(int n);
int compute_incr(int u);
int compute_decr(int u);
int fastWorstNode();
int fastBestNode(const vector<int> &nodes);
int newFastWorstNode(Configuration &c);
int fastCompIncr(Configuration &c, int u);
void updateConfiguration(Configuration &c, int u, int v);
void dfs(int u);
void externalDfs(int n);
int computeConf();
void deleteVertex(int v);
void restoreVertex(int v);
int compIncr(int u);
int compDecr(int u);
int wnode();
int wnodei(const vector<int> &S);
int probWnodei(const vector<int> &S, double prob);
int randomLargeCC();
int computeCost();
void update(int u, int v);
pair<int,int> biggestDecr(int v);
bool bestSwap();
void articulationPoints(int start, int d, vector<int> &low, vector<int> &depth, vector<int> &parent, int &treesize, int ccsize);
void getArticulationPonints();
void biggesDecrement(int u, int &best, int &bestDecr);
int bnode();


