#include <stdio.h>
#include <string.h>
#include <vector>
#include <queue>
#include <iostream> 
#include <algorithm> 
#include <utility> 
#include <map>
#include <fstream> 
#include <sstream>
#include <time.h>
#include <math.h>
using namespace std;

struct edge {
	int u;
	int v;
	
}; 

struct vertex {

	int index;
	std::vector<int> adj;
	
	int numEdges;
	
	int w;
	int root;
	int rank;
	int d;
	int f;
	int h;
	int color;
	int visited;
	int pai;
	
};

struct param
{
	int clusters;
	
};




