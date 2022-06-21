#include <iostream>
using namespace std;

vector<vector<int> > adjlist;
bool visited[100];

//ccompute the size of the tree using DFS
void treeSize(int v, int &size)
{
    visited[v] = true;
    size++;
    for(int i = 0; i < adjlist[v].size(); i++)
    {
        if(!visited[adjlist[v][i]])
            treeSize(adjlist[v][i], size);
    }
    
}

int main(){

}