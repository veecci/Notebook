.. _graph:

*************
Graph
*************

.. _graphd:

Graph
=====

.. code-block:: cpp

  int const N = 100010;
  int const M = 1000010;
  struct Edges {
    int u, next;
  } e[M];
  int p[N], idx;
  void init() { clr(p, 0xff); idx = 0; }
  void addedge(int u, int v) { e[idx].u = v, e[idx].next = p[u], p[u] = idx++; }

.. _cut_vertex_and_bridge:

Cut-Vertex and Bridge
=====================

.. code-block:: cpp

	//+ DSU/Stack to maintain the size of each component
	//* Cut-Vertex: (u) if (sc[u] > 1)
	//* Cut-Bridge: (u, v) if (low[v] > dfn[u]) *warning: parallel_edges
	struct CutVertex {
	  int dfn[N], low[N], sc[N], cnt;
	  void dfs(int u, int ori, int pre) {
	    dfn[u] = low[u] = ++cnt;
	    for (int i = p[u]; ~i; i = e[i].next) {
	      int v = e[i].u;
	      if (!dfn[v]) {
	        dfs(v, ori, i);
	        low[u] = min(low[u], low[v]);
	        if (low[v] >= dfn[u]) ++sc[u];
	      }
	      else if ((i ^ 1) != pre) { // ensure (e[(u, v)] ^ e[(v, u)]) == 1
	        low[u] = min(low[u], dfn[v]);
	      }
	    }
	    if (u != ori) ++sc[u];
	  }
	  void solve() {
	    cnt = 0, clr(dfn, 0), clr(sc, 0);
	    Rep(i, n) if (!dfn[i]) dfs(i, i, -1);
	  }
	};

.. _bcc:

BCC
=====================

.. image:: _static/bcc.png

.. code-block:: cpp

	// BCC(Edge)
	int col[N], cc; // clr(col, 0), cc = 0;
	void dfs(int u) {
	  if (!col[u]) col[u] = ++cc;
	  for (int i = p[u]; ~i; i = e[i].next) {
	    int v = e[i].u;
	    if (!col[v]) {
	      if (cv.low[v] > cv.dfn[u]) col[v] = ++cc;
	      else col[v] = col[u];
	      dfs(v);
	    }
	  }
	}
	// BCC(Vertex)
	bool vis[N]; // clr(vis, 0);
	int col[M], cc; // clr(col, 0), cc = 0;
	void dfs(int u, int pc) {
	  vis[u] = 1;
	  for (int i = p[u]; ~i; i = e[i].next) {
	    if (!col[i]) {
	      int v = e[i].u;
	      if (cv.low[v] >= cv.dfn[u]) col[i] = col[i ^ 1] = ++cc;
	      else col[i] = col[i ^ 1] = pc;
	      if (!vis[v]) dfs(v, col[i]);
	    }
	  }
	}

.. _scc_tarjan:

SCC (tarjan)
=====================

.. code-block:: cpp

	struct SCC {
	  int top, cnt, cc, t;
	  int st[N], dfn[N], low[N], col[N]; bool vis[N];
	  void tarjan(int u) {
	    dfn[u] = low[u] = ++cnt;
	    st[++top] = u, vis[u] = 1;
	    for (int i = p[u]; ~i; i = e[i].next) {
	      int v = e[i].u;
	      if (!dfn[v]) {
	        tarjan(v);
	        low[u] = min(low[u], low[v]);
	      }
	      else if (vis[v]) low[u] = min(low[u], dfn[v]);
	    }
	    if (dfn[u] == low[u]) {
	      do {
	        t = st[top--];
	        col[t] = cc;
	        vis[t] = 0;
	      } while (t != u);
	      ++cc;
	    }
	  }
	  void solve() {
	    top = cnt = cc = 0;
	    clr(vis, 0), clr(col, 0), clr(dfn, 0);
	    Rep(i, n) if (!dfn[i]) tarjan(i);
	  }
	};

.. _floyd:

SCC (floyd)
=====================

.. code-block:: cpp

	int n, mp[N][N]; // clr(mp, 0x3f); mp[i][i] = 0;
	void floyd() {
	  rep(k, n) rep(i, n) rep(j, n)
	    mp[i][j] = min(mp[i][j], mp[i][k] + mp[k][j]);
	}


.. _dijkstra:

dijkstra
=====================

.. code-block:: cpp

	priority_queue<pair<int, int> > Q;
	int dis[N]; bool vis[N];

	void dijkstra(int s) {
	  int u, v, w;
	  while (!Q.empty()) Q.pop(); clr(vis, 0), clr(dis, 0x3f);
	  Q.push(make_pair(0, s)), dis[s] = 0;
	  while (!Q.empty()) {
	    pair<int, int> tmp = Q.top(); Q.pop();
	    if (vis[u = tmp.second]) continue;
	    else vis[u] = true;
	    for (int i = p[u]; ~i; i = e[i].next) {
	      v = e[i].u, w = e[i].w;
	      if (!vis[v] && dis[u] + w < dis[v]) {
	        dis[v] = dis[u] + w;
	        Q.push(make_pair(-dis[v], v));
	      }
	    }
	  }
	}

.. _spfa:

SPFA
=====================

.. code-block:: cpp

	int dis[N]; bool vis[N];
	int Q[N * N];

	void spfa(int s) {
	  int u, v, w, l(0), h(0);
	  clr(vis, 0), clr(dis, 0x3f);
	  Q[h++] = s, dis[s] = 0;
	  while (l < h) {
	    u = Q[l++];
	    vis[u] = 0;
	    for (int i = p[u]; ~i; i = e[i].next) {
	      v = e[i].u, w = e[i].w;
	      if (dis[u] + w < dis[v]) {
	        dis[v] = dis[u] + w;
	        if (!vis[v]) {
	          vis[v] = 1;
	          Q[h++] = v;
	        }
	      }
	    }
	  }
	}

.. _hungary:

Hungary
=====================

.. code-block:: cpp

	/* Matching
	|Minimum Vertex Cover| = |Maximum Matching|
	|Maximum Independent Set| = |V| - |Maximum Matching|
	|Minimum Path Cover| = |V| - |Maximum Matching| (Directed Acyclic Graph)
	|Minimum Edge Cover| = |V| - |Maximum Matching| / 2 (Undirected Graph)
	 * warning: isolate vertex */
	// O(|V|*|E|)
	int n, m; // |V(x)|, |V(y)|
	int mp[N][N], matx[N], maty[N]; bool fy[N];

	int path(int u) {
	  rep(v, m) if (mp[u][v] && !fy[v]) {
	    fy[v] = 1;
	    if (!~maty[v] || path(maty[v])) {
	      matx[u] = v, maty[v] = u;
	      return 1;
	    }
	  }
	  return 0;
	}
	int hungary() {
	  int ret = 0;
	  clr(matx, 0xff), clr(maty, 0xff);
	  rep(i, n) if (!~matx[i]) {
	    clr(fy, 0);
	    ret += path(i);
	  }
	  return ret;
	}

.. _prim:

Prim
=====================

.. code-block:: cpp

	int n; // |V|
	int mp[N][N], ml[N]; bool vis[N];

	int prim() {
	  int ret(0); clr(vis, 0), clr(ml, 0x3f); ml[0] = 0;
	  rep(i, n) {
	    int id(-1);
	    rep(j, n) if (!vis[j] && (!~id || ml[j] < ml[id])) id = j;
	    vis[id] = 1, ret += ml[id];
	    rep(j, n) if (!vis[j] && mp[j][id] < ml[j]) ml[j] = mp[j][id];
	  }
	  return ret;
	}

.. _isap:

ISAP
=====================

.. code-block:: cpp

	//* bfs in the beginning to accelerate (+5%)
	int gap[N], lev[N], cur[N], pre[N];
	int sap(int s, int t) {
	  clr(gap, 0), clr(lev, 0), memcpy(cur, p, sizeof p);
	  int u, v, ret(0), step(inf), mi;
	  gap[0] = n, u = pre[s] = s;
	  while (lev[s] < n) { loop:
	    for (int &i = cur[u]; ~i; i = e[i].next) {
	      v = e[i].u;
	      if (e[i].c && lev[u] == lev[v] + 1) {
	        step = min(step, e[i].c);
	        pre[v] = u;
	        u = v;
	        if (v == t) {
	          while (v != s) {
	            u = pre[v];
	            e[cur[u]].c -= step;
	            e[cur[u] ^ 1].c += step;
	            v = u;
	          }
	          ret += step;
	          step = inf;
	        }
	        goto loop;
	      }
	    }
	    mi = n;
	    for (int i = p[u]; ~i; i = e[i].next) {
	      v = e[i].u;
	      if (e[i].c && lev[v] < mi) {
	        mi = lev[v];
	        cur[u] = i;
	      }
	    }
	    if (!--gap[lev[u]]) break;
	    ++gap[lev[u] = mi + 1];
	    u = pre[u];
	  }
	  return ret;
	}

.. _dinic:

Dinic
=====================

.. code-block:: cpp

	int cur[N], lev[N], Q[N], pre[N], st[N];
	bool bfs(int s, int t) {
	  clr(lev, 0xff); lev[s] = 0, Q[0] = s;
	  for (int l = 0, h = 1, u, v; l < h; ) {
	    u = Q[l++]; if (u == t) break;
	    for (int i = p[u]; ~i; i = e[i].next) {
	      v = e[i].u;
	      if (e[i].c > 0 && !~lev[v]) {
	        lev[v] = lev[u] + 1;
	        Q[h++] = v;
	      }
	    }
	  }
	  return lev[t] != -1;
	}
	int dfs(int s, int t) {
	  int u, v, top = 0, ret = 0, step, pos;
	  st[++top] = s;
	  while (top) {
	    u = st[top];
	    if (u == t) {
	      pos = 1;
	      Rep(i, top - 1) if (e[pre[i]].c < e[pre[pos]].c) pos = i;
	      ret += (step = e[pre[pos]].c);
	      Rep(i, top - 1) e[pre[i]].c -= step, e[pre[i] ^ 1].c += step;
	      u = st[top = pos];
	    }
	    for (int &i = cur[u]; ~i; i = e[i].next) {
	      v = e[i].u;
	      if (e[i].c && lev[u] < lev[v]) {
	        pre[top] = i, st[++top] = u = v;
	        break;
	      }
	    }
	    if (!~cur[u]) lev[u] = -1, --top;
	  }
	  return ret;
	}
	int dinic(int s, int t) {
	  int ret = 0, c;
	  while (bfs(s, t)) {
	    memcpy(cur, p, sizeof p);
	    ret += dfs(s, t);
	  }
	  return ret;
	}

.. _min_cost_max_flow:

MinCostMaxFlow
=====================

.. code-block:: cpp
	
	//addedge(u, v, cap,  cost);
	//addedge(v, u,   0, -cost);
	//*warning: no Negative-Cycle
	struct SSP {
	  int mc, mf;
	  int pre[N][2], Q[N], dis[N]; bool vis[N];
	  bool spfa(int s, int t) {
	    clr(dis, 0x3f), clr(vis, 0); dis[s] = 0, Q[0] = s, vis[s] = 1;
	    int u, v, w;
	    for (int l = 0, h = 1; l != h; ) {
	      vis[u = Q[l++]] = 0; if (l == N) l = 0;
	      for (int i = p[u]; ~i; i = e[i].next) {
	        v = e[i].u, w = e[i].f;
	        if (e[i].c && dis[u] + w < dis[v]) {
	          dis[v] = dis[u] + w;
	          pre[v][0] = u, pre[v][1] = i;
	          if (!vis[v]) {
	            vis[v] = 1;
	            Q[h++] = v; if (h == N) h = 0;
	          }
	        }
	      }
	    }
	    return dis[t] != inf; // return dis[t] < 0: any flow
	  }
	  void solve(int s, int t) {
	    mc = mf = 0; int u, step;
	    while (spfa(s, t)) {
	      step = inf;
	      for (u = t; u != s; u = pre[u][0]) {
	        step = min(step, e[pre[u][1]].c);
	      }
	      for (u = t; u != s; u = pre[u][0]) {
	        e[pre[u][1]].c -= step;
	        e[pre[u][1] ^ 1].c += step;
	      }
	      mf += step;
	      mc += dis[t] * step;
	    }
	  }
	};

.. _min_cut:

MinCut(UndirectedGraph)
=====================

.. code-block:: cpp

	int v[N], mp[N][N], dis[N]; bool vis[N];
	int Stoer_Wagner(int n) { // 0 ~ n-1
	  int ret = inf, now, pre;
	  rep(i, n) v[i] = i;
	  while (n > 1) {
	    now = 1, pre = 0;
	    Rep(i, n - 1) {
	      dis[v[i]] = mp[v[0]][v[i]];
	      if (dis[v[i]] > dis[v[now]]) now = i;
	    }
	    clr(vis, 0), vis[v[0]] = 1;
	    Rep(i, n - 2) {
	      vis[v[now]] = 1, pre = now, now = -1;
	      Rep(j, n - 1) if (!vis[v[j]]) {
	        dis[v[j]] += mp[v[pre]][v[j]];
	        if (now == -1 || dis[v[now]] < dis[v[j]]) now = j;
	      }
	    }
	    ret = min(ret, dis[v[now]]);
	    rep(i, n) {
	      mp[v[pre]][v[i]] += mp[v[now]][v[i]];
	      mp[v[i]][v[pre]] = mp[v[pre]][v[i]];
	    }
	    v[now] = v[--n];
	  }
	  return ret;
	}