// Densest subgraph for weighted graphs
// Uses BBST to store intermediate degrees
// 6-10 x slower than using linked lists, but simpler to understand

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
#include <chrono>
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////
// Helper for fast input

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

////////////////////////////////////////////////////////////////////////////////////////

vector<double> deg;

struct classcomp {
  bool operator() (const int& lhs, const int& rhs) const
  {return deg[lhs]<deg[rhs] || (deg[lhs]==deg[rhs] && lhs<rhs);}
};

set<int,classcomp> deg_sorted; //BBST storing degrees

//////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {

  cout << "Finding maximum subgraph density (naive version using BBST)..." << endl;
  
  auto startio = chrono::steady_clock::now();

  float theta = atof(argv[1]);
  float accuracy = atof(argv[2]);
  float d2_bound = atof(argv[3]);
  char* edgelist=argv[4];

  cout << theta << edgelist << endl;
  FILE *file;
  file=fopen(edgelist,"r");
  int m, n;
  fscanf(file,"%u %u", &n, &m);
  fclose(file);

  double mm_density = 0.0;
  vector<bool> ans(n);

  double d2_l = 0, d2_r =d2_bound;
  while (d2_r-d2_l>accuracy) {
    double d2 = (d2_l + d2_r) / 2;
    cout << "d2 = " << d2 << endl;
    double * init_deg = new double[n];
    memset(init_deg, 0., sizeof(double) * n);
    
    double * w = new double[n];
    memset(w, 0, sizeof(double) * n);
    
    deg.reserve(n);

    set<pair<int, double>> * nbrs = new set<pair<int, double>>[n];

    double sum_wts = 0.0, ttl_node_weight = 0.0;

    // int *node_weight = new int[n];

    cout<< "read node weights"<< endl;
    file=fopen(edgelist,"r");
    fscanf(file,"%u %u", &n, &m);
    // cout << n << ", " <<  m << endl;
    for (int i = 0; i < n; ++i) {
      int input_w;
      fscanf(file,"%u", &input_w);
      double nw = input_w*d2;
      w[i] = nw;
      ttl_node_weight += nw;
      // cout << nw << endl;
    }
    cout << ttl_node_weight << endl;
    cout<< "read edge weights"<< endl;
    for (int i = 0; i < m; i++) {
      int p, q, input_w;
      double wt=0;
      fscanf(file,"%u %u %u", &p, &q, &input_w);
    	wt = input_w;
      nbrs[p].insert(make_pair(q,wt));
      nbrs[q].insert(make_pair(p,wt));
      // cout << "add deg" << endl;
      init_deg[p]+=wt;
      init_deg[q]+=wt;
      sum_wts += wt;
    }
    fclose(file);

    vector<bool> exists(n);
    vector<bool> iterans(n);

    auto endio = chrono::steady_clock::now();
    int init_time = chrono::duration_cast<chrono::milliseconds>(endio - startio).count();
    
    cout << "Time for reading input and initialization: " << init_time << " ms" << endl;
    
    int sum_iter_times = 0;

    // for (int tt = 0; tt < iters; tt++) {

  	auto startiter = chrono::steady_clock::now();

  	deg_sorted.clear();

    for (int i = 0; i < n; i++) {
      deg[i] = w[i] + init_deg[i]; //degree for this iteration is "vertex weight" + actual degree
      deg_sorted.insert(i);
    }
        
    double iter_max_density = 0;
    if (ttl_node_weight/d2/n > theta)
      iter_max_density = (double) sum_wts / n;
    double cur_sum_wts = sum_wts, removed_node_w = 0.0, scaled_nw=0, tight_deg=0, cur_size=0, best_size=0;
    int cur_n = n;
    cout << "edge sum " << cur_sum_wts << endl;
    fill(exists.begin(), exists.end(), true);
    iterans = exists;

    while (cur_n > 1) {
      cur_n--;
      int k = *(deg_sorted.begin()); //k = min degree vertex

      // w[k] = deg[k]; //increment vertex weight for the next iteration (self loops)
      deg_sorted.erase(k); //delete k
      removed_node_w += w[k];
      if (deg[k]>tight_deg){
        tight_deg = deg[k];
        scaled_nw = w[k];
      }
      for (pair<int, double> j : nbrs[k]) { //decrement degrees of k's neighbors
        int nbr = j.first;
    	  double nbrwt = j.second;
        if (exists[nbr]) {
          deg_sorted.erase(nbr);
          deg[nbr] -= nbrwt;
          cur_sum_wts -= nbrwt;
          deg_sorted.insert(nbr);
        }
        exists[k] = false;
      }
      if ((iter_max_density < (double) cur_sum_wts / cur_n) && ((double)(ttl_node_weight-removed_node_w)/d2/cur_n > theta)) {
        iter_max_density = (double) cur_sum_wts / cur_n;
        iterans = exists;
        cur_size = cur_n;
      }

      // if iter_max_density
    }
    
    if((iter_max_density > mm_density)) {
      mm_density = iter_max_density;
      ans = iterans;
      best_size = cur_size;
    }

    auto enditer = chrono::steady_clock::now();
    int elapsed = chrono::duration_cast<chrono::milliseconds>(enditer - startiter).count();
    cout << cur_sum_wts << ", node weight left = " << ttl_node_weight-removed_node_w << endl;
    cout << "Max density: " << iter_max_density << " with size:" << cur_size << endl;
    cout << "Max density so far" <<": " << mm_density << " with size:" << best_size << endl;
    // cout << "Avg time per iteration: " << sum_iter_times/(tt+1) << " ms" << endl;
    cout << "iteration time: " << elapsed << " ms" << endl;

    // condition idea 1, change based on c_u
    if ( (iter_max_density<=0) || (scaled_nw/d2)<theta ){
      d2_l = d2;
    } else {
      d2_r = d2;
    }

    // condition idea 2, based on 
  }

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
  for (int i=0;i<n;i++)
  {
    if (ans[i]) outfile << i << endl;
  }

  return 0;
}