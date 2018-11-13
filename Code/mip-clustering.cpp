#include <ilcplex/ilocplex.h>
#include "estruturas.h"
#include "mersenneTwister.h"
#include <climits>
using namespace std;   

std::vector<vertex> V;
std::vector<edge> E;
int n;
int c;
std::vector<vector<int>> D;
param param;
#define EPS 10e-6

bool readInput(char *filename) {
	int i, u, v;
	string str1,str2,str3,line;
	ifstream inputFile (filename);
	
	
	if (inputFile.is_open()){

		getline (inputFile,line);
		getline (inputFile,line);
		getline (inputFile,line);

		inputFile>>str1;
		i = 0;
		while (str1.compare("node") == 0) {

			V.push_back(vertex());
			V[i].numEdges = 0;
			V[i].index = i;
			i++;
			//printf("%d\n", V.size());

			inputFile>>str2;
			while (str2.compare("]") != 0) {
				//cout << str2;
				inputFile>>str2;
			}
			inputFile>>str1;
		}

		//printf("%d\n", i);

		//inputFile>>str1;
		i = 0;
		while (str1.compare("edge") == 0) {
			
			inputFile>>str3;
			E.push_back(edge());
			inputFile>>str2>>u;
			inputFile>>str3>>v;
			u--;
			v--;

			V[u].numEdges++;
			V[v].numEdges++;

			
			//printf("%d %d\n", u, v);

			V[u].adj.push_back(v);
			V[v].adj.push_back(u);
			if (u > v) {
				E[i].u = v;
				E[i].v = u;
			} else {
				E[i].u = u;
				E[i].v = v;
			}
			i++;
			inputFile>>str3;
			inputFile>>str1;
			//cout << str1;
		}

	} else {
		cout << "Unable to open file"; 
		return false;	
	}

	return true;
}

string IntToString (int a)
{
    ostringstream temp;
    temp<<a;
    return temp.str();
}


void clusteringModel() {

	int i, j, l, k;
	string varName;

	IloEnv env;
	//IloNum start;

	IloModel model(env);

	//decision variables
	//n denotes the number of objects
	//c denotes the number of clusters
	IloArray<IloIntVarArray> X(env, n); 

	IloIntVar Dm(env); 
	varName = "Dm";
    Dm.setName(varName.c_str());
                        
	IloIntVar Z(env); 
	varName = "Z";
    Z.setName(varName.c_str());

	for (i = 0; i < n; i++) {
		X[i] = IloIntVarArray(env, c, 0, 1);
		for (int v = 0; v < c; v++) {
		varName = "X("+IntToString(i) +","+IntToString(v)+")";
                        	X[i][v].setName(varName.c_str());
                        }
	}


	//(1) Objective function
	model.add(IloMinimize(env, Dm + Z));

	//(2)
	int limi = E.size();
	for (k = 0; k < limi; k++) {
		i = E[k].u;
		j = E[k].v;
		//printf("%d %d %d\n",n, i, j);
		for (l = 0; l < c; l++) {
			model.add(Dm >= D[i][j]*(X[i][l] + X[j][l] - 1));
		}
	}

	//(3)
	for (i = 0; i < n; i++) {
		IloExpr sum(env);
		for (l = 0; l < c; l++) {
			sum += X[i][l];
		}
		model.add(sum == 1);
		sum.end();
	}

	//(4)
	for (l = 0; l < c; l++) {
		IloExpr sum(env);
		for (i = 0; i < n; i++) {
			sum += X[i][l];
		}
		model.add(sum >= 1);
		sum.end();
	}

	//(5)
	for (i = 0; i < n; i++) {
		for (l = 0; l < c; l++) {
			IloExpr sum(env);
			IloExpr sum2(env);
			for (j = 0; j < V[i].numEdges; j++) {
				sum += X[V[i].adj[j]][l];
			}
			sum2 = (V[i].numEdges)*X[i][l];
			sum2 /= 2;
			model.add(sum >= sum2);
			sum.end();
			sum2.end();
		}

	}

	//(6)
	for (i = 0; i < n; i++) {
		
		for (l = 0; l < c; l++) {
			IloExpr sum(env);
			for (j = 0; j < V[i].numEdges; j++) {
				sum += X[j][l];
			}
			model.add(Z >= (V[i].numEdges*X[i][l] - sum));
			sum.end();
		}

	}

	double elapsed3;

	IloCplex cplex2(model);
	cplex2.exportModel("model.lp");
	
	cplex2.setParam(IloCplex::TiLim, 180);
	cplex2.setParam(IloCplex::EpAGap, 0.9);
	time_t start3 = time(NULL);

	if ( !cplex2.solve() ) {
		env.error() << "Failed to optimize LP." << endl;
		throw(-1);
	}

	env.out() << "Solution status = " << cplex2.getStatus() << endl;
	env.out() << "Solution value = " << cplex2.getObjValue() << endl;
	time_t end3 = time(NULL);
	elapsed3 = difftime(end3,start3);
	env.out() << "Solution time= " <<  elapsed3 << endl;
}

void FloydWarshall() {
	int i, j, k;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			D[i][j] = 500000;
		}
	}

	for (i = 0; i < n; i++) {
		D[i][i] = 0;
	}

	for (i = 0; i < E.size(); i++) {
		D[E[i].u][E[i].v] = 1;
		D[E[i].v][E[i].u] = 1;
	}

	for (k = 0; k < n; k++) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (D[i][j] > (D[i][k] + D[k][j]))
					D[i][j] = D[i][k] + D[k][j];
			}
		}
	}
}

int main(int argc, char **argv) {

	int i, j;

	param.clusters = atoi(argv[2]);
	c = param.clusters;

    //clock_t start = clock();




   
	readInput(argv[1]);

	n = V.size();

	D.resize(n);
	for (i = 0; i < n; i++) {
		D[i].resize(n);
	}

	FloydWarshall();

	/*for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%d ", D[i][j]);
		}
		printf("\n");
	}*/

	

	/*for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			D[i][j] = random(1,10);
		}
	}*/

	//printf("%d %d\n", V.size(), E.size());
	
    //clock_t end = clock();
	//etime = (double) (end-start) / CLOCKS_PER_SEC;


    clusteringModel();
	    	
	    

	   

	return 0;	
}
