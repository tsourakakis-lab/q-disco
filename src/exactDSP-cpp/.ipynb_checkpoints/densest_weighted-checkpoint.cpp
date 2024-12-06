extern "C" {
#include "hi_pr.h"
}

//#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include<cassert>
#include <vector>   
#include <queue> 
#include<list>
#include <set>
#include<cstring>
#include<ctime>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <vector> 
#include <cmath> 
#include <chrono>
using namespace std;


const int INF = (int)1e9;

void add_arc(int u, int v, int capacity, node *nodes, arc *arcs, cType *cap, int *cur_arc) {
    arcs[cur_arc[u]].head = nodes + v;
    arcs[cur_arc[u]].rev = arcs + cur_arc[v];
    arcs[cur_arc[v]].head = nodes + u;
    arcs[cur_arc[v]].rev = arcs + cur_arc[u];
    cap[cur_arc[u]] = capacity;
    cap[cur_arc[v]] = 0;
    ++cur_arc[u];
    ++cur_arc[v];
}

bool nontrivial(int n_nodes, node *nodes) {
    int res = 0;
    for (int i = 0; i < n_nodes; ++i) {
        if (nodes[i].d < n_nodes) {
            ++res;
        }
    }
    return res > 1;
}

inline char GET_CHAR(){
	const int maxn = 131072;
	static char buf[maxn],*p1=buf,*p2=buf;
	return p1==p2&&(p2=(p1=buf)+fread(buf,1,maxn,stdin),p1==p2)?EOF:*p1++;
}
inline int getInt() {
	int res(0);
	char c = GET_CHAR();
	while(c < '0') c = GET_CHAR();
	while(c >= '0') {
		res = res * 10 + (c - '0');
		c = GET_CHAR();
	}
	return res;
}

int main(int argc, char **argv) {
//THIS NUMBER DETERMINES HOW MANY DECIMALS WE COMPUTE THINGS TO
//  HI_PR ONLY TAKES INTEGER CAPACITIES
//    SO WE CAN'T PUT FRACTIONALS, INSTEAD, JUST MULTIPLY EVERYTHING BY 100 OR 1000...

    // Input file start with n,m. Then each line contains an integer weight for a node. Finally each line contains u,v,w which represents the weight of each uv is w (integer as well).
    auto startio = chrono::steady_clock::now();
    int opt_density_estimate_upper;

    int ACCURACY = atoi(argv[1]);
    int lambda_acc = atoi(argv[2]);
    float theta = atof(argv[3]);
    char* edgelist=argv[4];

    FILE *file;
    file=fopen(edgelist,"r");
    string output_file;
    ofstream outfile;
    
    if (argc >= 6)
    {
        output_file = argv[5];
    }
    else
    {
        output_file = "soln.tmp";
    }

    outfile.open(output_file.c_str());

    cout << "input estimate optimal density:";
    cin >> opt_density_estimate_upper;
    cout<< "phase 0";
    
    int n, m, edge_weight_sum=0, node_weight_sum=0, min_node_weight=100000;
    int k = 2; ///EDGES CASE!
    fscanf(file,"%u %u", &n, &m);

    cout<< "phase 1";

    int *cliques = new int[m * k];
    int *weight = new int[m];
    int *node_weight = new int[n];

    int *weighted_deg = new int[n];

    cout<< "phase 1";
    for (int i = 0; i < n; ++i) {
      int w;
      fscanf(file,"%u", &w);
      node_weight[i] = w;
      min_node_weight = min(min_node_weight, abs(w));
      // max_node_weight = max(w, max_node_weight);
    }
    cout<< "phase 1";
    for (int i = 0; i < m; ++i) {
      int p, q, w;
      fscanf(file,"%u %u %u", &p, &q, &w);
      cliques[i*2]=p;
      cliques[i*2+1]=q;
      weight[i] = w;
      weighted_deg[p] += w;
      weighted_deg[q] += w;
      edge_weight_sum += w;
    }
    edge_weight_sum = min(opt_density_estimate_upper, edge_weight_sum);
    fclose(file);

    auto endio = chrono::steady_clock::now();
    
    cout << "Time for reading input: " << chrono::duration_cast<chrono::milliseconds>(endio - startio).count() << " ms" << endl;
   
    auto start = chrono::steady_clock::now();
    
    int source = n;
    int sink = n+1;
    int n_nodes = n+2;
    int *deg = new int[n_nodes];
    int max_deg=0, max_node_weight=0;
    deg[source] = n;
    deg[sink] = n;
    for (int i = 0; i < n; ++i) {
        deg[i] = 2;
        max_deg = max(max_deg, weighted_deg[i]);
        max_node_weight = max(max_node_weight, node_weight[i]);
    }
    for (int i = 0; i < m * k; ++i) {
        // bi-directional edges between nodes
        deg[cliques[i]] += 2;
    }
    for (int i = 1; i < n_nodes; ++i) {
        deg[i] += deg[i - 1];
    }
    min_node_weight = max(1, min_node_weight);
    cout << lambda_acc << "initialized. " << min_node_weight << endl;
    int lambda_l = 0, lambda_r = (int)(lambda_acc*edge_weight_sum/min_node_weight);
    cout << "lambda upper bound = " << lambda_r << endl;
    while (lambda_l<lambda_r) {
        int lambda = (lambda_l + lambda_r + 1) / 2;
        cout << "using lambda = " << lambda << endl;
        int n_arcs = deg[n_nodes - 1];

        int *cur_arc = new int[n_nodes];

        node *nodes_ = new node[n_nodes + 1];
        arc *arcs = new arc[n_arcs];
        cType *cap  = new cType[n_arcs];
        int additional_cap = (int)(lambda_acc*edge_weight_sum/min_node_weight)*(max_node_weight*2+max_deg+1)* ACCURACY;
        // cout << "cap added "<< additional_cap << ", max weight:" << max_node_weight<<", max deg:"<<max_deg << endl;
        for (int i = 0; i < n_nodes; ++i) {
            cur_arc[i] = (i == 0) ? 0 : deg[i - 1];
            nodes_[i].first = arcs + cur_arc[i];
        }

        for (int i = 0; i < n; ++i) {
            add_arc(i, sink, additional_cap, nodes_, arcs, cap, cur_arc);
            // cout << i << ", " << sink << ", " << additional_cap << endl;
        }

        for (int i = 0; i < n; ++i) {
            add_arc(source, i, 0, nodes_, arcs, cap, cur_arc);
        }
        for (int i = 0; i < m; ++i) {
            int p = cliques[k*i], q = cliques[k*i+1];
            add_arc(q, p, lambda_acc*weight[i]* ACCURACY, nodes_, arcs, cap, cur_arc);
            add_arc(p, q, lambda_acc*weight[i]* ACCURACY, nodes_, arcs, cap, cur_arc);
            // cout << p << ", " << q << ", " << weight[i]* ACCURACY << endl;
        }

        // upper and lower bound on density
        int l = 0, r = (edge_weight_sum*lambda_acc+max_node_weight*lambda) * ACCURACY;
        cout << "upper bound: " << r << endl;
        // cout << "start" << endl;
        vector<long> subg;
        node *j;
        while (l < r) {
            int c = (l + r + 1) / 2;
            // cout << "[" << l << ", " << r << "] " << c << endl;
            for (int i = 0; i < n; ++i) {
                // if (additional_cap + 2*c - ACCURACY*(lambda_acc*weighted_deg[i]+lambda*2*node_weight[i])<=0):
                //     cout << source << ", " << i << ", " << additional_cap << "!!!!!" << endl;
                cap[nodes_[source].first - arcs + i] = additional_cap + 2*c - ACCURACY*(lambda_acc*weighted_deg[i]+lambda*2*node_weight[i]);
            }
            min_cut(n_nodes, n_arcs / 2, nodes_, arcs, cap, nodes_ + source, nodes_ + sink, 0);
            if (nontrivial(n_nodes, nodes_)) {
                l = c;
                subg.clear();
                forAllNodes(j)
                    if (j->d < n_nodes && nNode(j) < n)
                        subg.push_back(nNode(j));
            } else {
                r = c - 1;
            }

        }
        
        auto end = chrono::steady_clock::now();
        
        cout << "Time for finding solution: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

        int res_node_weight=0;
        for (int i : subg)
        {
            res_node_weight += node_weight[i];
            if (subg.size()<20)
                cout << " " << i;
        }
            cout << endl;

        cout << "size:" << subg.size() << ", Density is " << (double)l/ACCURACY/lambda_acc << ", node weight density: " << float(res_node_weight)/subg.size() << ", edge density:" << (double)l/ACCURACY/lambda_acc - float(res_node_weight)/subg.size()*lambda/lambda_acc << endl;

        if ((subg.size()<=0) || (float(res_node_weight)/subg.size()>theta)){
            // decrease lambda
            lambda_r = lambda - 1;
        } else {
            // increase lambda
            lambda_l = lambda;
        }

        outfile << "lambda: " << (double)lambda/lambda_acc << ", size: " << subg.size() << ", node weight density: " << float(res_node_weight)/subg.size() << ", edge density:" << (double)l/ACCURACY/lambda_acc - float(res_node_weight)/subg.size()*lambda/lambda_acc << endl;
        if (n<=20) {
            for (int i : subg)
            {
                outfile << i << endl;
            }
        }
    }
    

    return 0;
}
