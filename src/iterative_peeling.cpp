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

int n, m;

struct Edge {
  int y, next;
};

struct Node {
  int deg, next, prev, idx;
  inline void clear() {
    deg = next = prev = 0;
    idx = -1;
  }
};

Node * lists;

__inline void linklists (int x, int y) {
  if(y == 0) return;
  lists[x].next = y;
  lists[y].prev = x;
};

int * nxt, * prv, *itr;

__inline void linknodes (int x, int y) {
  if(y == -1) return;
  nxt[x] = y;
  prv[y] = x;
};

__inline void eraselist(int x) {
  lists[lists[x].prev].next = lists[x].next;
  if(lists[x].next != 0) lists[lists[x].next].prev = lists[x].prev;
};

__inline void erasenode (int x) {
  if(prv[x] == -1) {
    lists[itr[x]].idx = nxt[x];
  }
  if(prv[x] != -1) nxt[prv[x]] = nxt[x];
  if(nxt[x] != -1) prv[nxt[x]] = prv[x];

};

int l = 0;
Edge * edges;
int * idx;

__inline void build(int x, int y) {
    edges[++l].next = idx[x];
    edges[l].y = y;
    idx[x] = l;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Main
int main(int argc, char** argv) {

  auto startio = chrono::steady_clock::now();

  float theta = atof(argv[1]);
  int edge_acc = atoi(argv[2]);
  int node_acc = atoi(argv[3]);
  char* edgelist=argv[4];

  // cout << theta << edgelist << endl;
  FILE *file;
  file=fopen(edgelist,"r");
  int m, n;
  fscanf(file,"%u %u", &n, &m);
  fclose(file);

  
  vector<int> m_ans;
  double mm_density = 0;

  double d2_l = 0, d2_r =node_acc;
  while (d2_r > d2_l) {
    int d2 = (d2_l + d2_r + 1) / 2;
    cout << "d2 = " << d2 << endl;
    edges = new Edge[m * 2 + 10];
    idx = new int[n];
    memset(idx, 0, sizeof(int) * n);
    int * init_deg = new int[n];
    memset(init_deg, 0, sizeof(int) * n);
    l = 0;
    
    lists = new Node[n + 2 * m + 10];
    int n_list = 0;
    itr = new int[n];
    int * deg = new int[n], * w = new int[n], * pos = new int [n];
    memset(deg, 0, sizeof(int) * n);
    memset(w, 0, sizeof(int) * n); //initial vertex weights=0, i.e., no self loops at the start 
    memset(pos, 0, sizeof(int) * n);
    prv = new int[n]; nxt = new int[n];

    // cout<< "read node weights"<< endl;
    file=fopen(edgelist,"r");
    int ttl_node_weight = 0;
    fscanf(file,"%u %u", &n, &m);
    // cout << n << ", " <<  m << endl;
    for (int i = 0; i < n; ++i) {
      int input_w;
      fscanf(file,"%u", &input_w);
      input_w = input_w*d2;
      w[i] = input_w;
      ttl_node_weight += input_w;
      // cout << nw << endl;
    }
    // cout << ttl_node_weight << endl;
    // cout<< "read edge weights"<< endl;
    // for (int i = 0; i < m; i++) {
    //   int p, q, input_w;
    //   double wt=0;
    //   fscanf(file,"%u %u %u", &p, &q, &input_w);
    //   wt = input_w;
    //   nbrs[p].insert(make_pair(q,wt));
    //   nbrs[q].insert(make_pair(p,wt));
    //   // cout << "add deg" << endl;
    //   init_deg[p]+=wt;
    //   init_deg[q]+=wt;
    //   sum_wts += wt;
    // }
    // fclose(file);

    for (int i = 0; i < m; i++) {
      int p, q, input_w;
      fscanf(file,"%u %u %u", &p, &q, &input_w);
      build(p, q);
      build(q, p);
      init_deg[p] += edge_acc;
      init_deg[q] += edge_acc;
    }
    fclose(file);

    pair<int, int> * deg_sorted = new pair<int, int>[n];

    auto endio = chrono::steady_clock::now();
    int sum_iter_times = 0;
    int init_time = chrono::duration_cast<chrono::milliseconds>(endio - startio).count();
    
    cout << "Time for reading input and initialization: " << init_time << " ms" << endl;
    

  	auto startiter = chrono::steady_clock::now();

    for (int i = 0; i < n; i++) {
      nxt[i] = prv[i] = -1;
      pos[i] = 0;
      
      deg[i] = w[i] + init_deg[i]; //degree for this iteration is "vertex weight" + actual degree
      deg_sorted[i] = make_pair(deg[i], i);
    }
    sort(deg_sorted, deg_sorted + n);
    n_list = 0;

    for(int i = 0; i < n; i++) {
      int v = deg_sorted[i].second;
      if(n_list == 0 || lists[n_list].deg != deg_sorted[i].first) {
        ++n_list;
        lists[n_list].clear();
        linklists(n_list - 1, n_list);
        lists[n_list].deg = deg_sorted[i].first;
      }
      
      linknodes(v, lists[n_list].idx);
      lists[n_list].idx = v;
      itr[v] = n_list;
    }
    
    int max_size = 0;
    double max_density = 0, node_density = 0;
    if ((double)ttl_node_weight/d2/n > theta){
      max_density = (double) m / n;
      node_density = (double)ttl_node_weight/d2/n;
      max_size = n;
    }

    int cur_m = m, cur_n = n, scaled_nw=0, tight_deg=0, removed_node_w = 0;
    vector<int> ans;
    while(lists[0].next) {
      
      int i = lists[0].next;
      int k = lists[i].idx;
      
      if(nxt[k] == -1) {
        eraselist(i);
      }else {
        erasenode(k);
      }
      pos[k] = -1;
      removed_node_w += w[k];
      if (deg[k]>tight_deg){
        tight_deg = deg[k];
        scaled_nw = w[k];
      }
      // w[k] = deg[k]; //increment vertex weight for the next iteration (self loops)
      cur_n -= 1;
      ans.push_back(k);
      for (int p = idx[k]; p; p = edges[p].next) { //decrement degrees of k's neighbors
        int j = edges[p].y;
        if(pos[j] == -1) continue;
        cur_m -= 1;
        
        int i = itr[j];
        erasenode(j);
        int i1 = lists[i].prev;
        
        if(lists[i].idx == -1) eraselist(i);
        deg[j]-= edge_acc;       
        prv[j] = nxt[j] = -1;
        while (i1!=0 && lists[i1].deg > deg[j])
          i1 = lists[i1].prev;
        if(i1 == 0 || lists[i1].deg != deg[j]) {
          ++n_list;
          lists[n_list].clear();
          itr[j] = n_list;
          int i2 = lists[i1].next;
          lists[n_list].deg = deg[j];
          lists[n_list].idx = j;
          linklists(i1, n_list);
          if(i2) linklists(n_list, i2);
        }
        else {
          linknodes(j, lists[i1].idx);
          lists[i1].idx = j;
          itr[j] = i1;
        }
      }
      if(cur_n == 0) continue;
      if((max_density < (double)cur_m / cur_n) && ((double)(ttl_node_weight-removed_node_w)/d2/cur_n > theta)) {
        max_size = cur_n;
        max_density = (double)cur_m / cur_n;
        node_density = (double)(ttl_node_weight-removed_node_w)/d2/cur_n;
      }
      
    }
    reverse(ans.begin(), ans.end());
    ans.resize(max_size);
    if(max_density > mm_density) {
      m_ans = ans;
      mm_density = max_density;
    }

    auto enditer = chrono::steady_clock::now();
    int elapsed = chrono::duration_cast<chrono::milliseconds>(enditer - startiter).count();
    sum_iter_times += elapsed;

    // cout << cur_m << ", node weight left = " << ttl_node_weight-removed_node_w << endl;
    cout << "Max density: " << max_density << " with size:" << max_size << " and node density:" << node_density << endl;
    cout << "Max density so far" <<": " << mm_density << endl;
    // cout << "Avg time per iteration: " << sum_iter_times/(tt+1) << " ms" << endl;
    cout << "iteration time: " << elapsed << " ms" << endl;
  // condition idea 1, change based on c_u
    if ( (max_density<=0) || (double)(scaled_nw/d2)<theta ){
      d2_l = d2;
    } else {
      d2_r = d2 - 1;
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
  for (int i : m_ans)
  {
    outfile << i << endl;
  }

  return 0;
}