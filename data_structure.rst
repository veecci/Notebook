.. _data_structure:

***************
Data Structure
***************

.. _dsu:

DSU
=====

.. code-block:: cpp

  int parents[N];
  int find(int a) { return parents[a] < 0 ? a : parents[a] = find(parents[a]); }
  void uni(int a, int b) {
    a = find(a), b = find(b); if (a == b) return;
    if (parents[a] < parents[b]) { parents[a] += parents[b], parents[b] = a; }
    else { parents[b] += parents[a], parents[a] = b; }
  }
  void init() { clr(parents, 0xff); }

.. _dsu_complex:

DSU(Vector)
===========

.. code-block:: cpp

	int fa[N], v[N];
	int Find(int a) {
	  if (fa[a] < 0) return a;
	  else {
	    int t = fa[a];
	    fa[a] = Find(fa[a]);
	    v[a] += v[t];
	    return fa[a];
	  }
	}
	//addEdge(b, a, c) -> Union(ra, rb, v[a] - v[b] + c)
	void Union(int a, int b, int c = 0) {
	  if (fa[a] < fa[b]) {
	    fa[a] += fa[b], fa[b] = a, v[b] += c;
	  } else {
	    fa[b] += fa[a], fa[a] = b, v[a] -= c;
	  }
	}
	void init() { clr(fa, 0xff), clr(v, 0); }

.. _fenwick_tree:

FenwickTree
============

.. code-block:: cpp

	struct FenwickTree {
	  ll a[N];
	  inline void init() { clr(a, 0); }
	  inline int lowbit(int x) { return x & -x; }
	  void update(ll p, ll c) {
	    while (p < N) {
	      a[p] += c;
	      p += lowbit(p);
	    }
	  }
	  ll query(ll p) {
	    ll ret = 0;
	    while (p > 0) {
	      ret += a[p];
	      p -= lowbit(p);
	    }
	    return ret;
	  }
	  int get_kth(ll k) {
	    int now = 0;
	    for(int i = 20; ~i; --i) { // for N ~ 1e6
	      now |= (1 << i);
	      if (now >= N || a[now] >= k)
	        now ^= (1 << i);
	      else k -= a[now];
	    }
	    return now + 1;
	  }
	};

.. _rmq:

RMQ
============

.. code-block:: cpp

	//+ pos if needed
	struct RMQ {
	  int lg[N], dmx[M][N]; // M = lg[N] + 1
	  void init() {
	    lg[0] = -1;
	    Rep(i, n) {
	      lg[i] = lg[i - 1] + !(i & (i - 1));
	      dmx[0][i] = a[i]; // the original array
	    }
	    for (int i = 1; (1 << i) <= n; ++i) {
	      for (int j = 1; j + (1 << i) - 1 <= n; ++j) {
	        dmx[i][j] = max(dmx[i - 1][j], dmx[i - 1][j + (1 << i - 1)]);
	      }
	    }
	  }
	  int get_mx(int a, int b) { // a <= b
	    int k = lg[b - a + 1];
	    return max(dmx[k][a], dmx[k][b - (1 << k) + 1]);
	  }
	};

.. _rmq_2d:

RMQ(2D)
============

.. code-block:: cpp

	int a[N][N], n, m;
	struct RMQ2D {
	  int lg[N], dmx[M][M][N][N];
	  void init() {
	    lg[0] = -1;
	    Rep(i, N - 1) lg[i] = lg[i - 1] + !(i & (i - 1));
	    Rep(i, n) Rep(j, m) dmx[0][0][i][j] = a[i][j];
	    rep(_i, lg[n] + 1) rep(_j, lg[m] + 1) if (_i || _j) {
	      Rep(i, n + 1 - (1 << _i)) Rep(j, m + 1 - (1 << _j)) {
	        if (_i == 0) dmx[_i][_j][i][j] =
	          max(dmx[_i][_j - 1][i][j], dmx[_i][_j - 1][i][j + (1 << _j - 1)]);
	        else dmx[_i][_j][i][j] =
	          max(dmx[_i - 1][_j][i][j], dmx[_i - 1][_j][i + (1 << _i - 1)][j]);
	      }
	    }
	  }
	  int get_mx(int x1, int y1, int x2, int y2) {
	    int kx = lg[x2 - x1 + 1], ky = lg[y2 - y1 + 1];
	    int m1 = dmx[kx][ky][x1][y1];
	    int m2 = dmx[kx][ky][x2 - (1 << kx) + 1][y1];
	    int m3 = dmx[kx][ky][x1][y2 - (1 << ky) + 1];
	    int m4 = dmx[kx][ky][x2 - (1 << kx) + 1][y2 - (1 << ky) + 1];
	    return max(max(max(m1, m2), m3), m4);
	  }
	};

.. _lca:

LCA(online)
============

.. code-block:: cpp

	int dep[N], dp[21][N]; // N ~ 1e6

	void dfs(int u, int d, int pre) {
	  dep[u] = d, dp[0][u] = pre;
	  for (int i = 1; (1 << i) <= d; ++i)
	    dp[i][u] = dp[i - 1][dp[i - 1][u]];
	  for (int i = p[u]; ~i; i = e[i].next) {
	    int v = e[i].u;
	    if (v != pre) dfs(v, d + 1, u);
	  }
	}
	int lca(int a, int b) {
	  if (dep[a] < dep[b]) swap(a, b);
	  for (int t = dep[a] - dep[b], step = 0; t; ++step, t >>= 1) 
	    if (t & 1) a = dp[step][a];
	  if (a == b) return a;
	  for (int i = 20; ~i; --i) 
	    if (dp[i][a] != dp[i][b]) a = dp[i][a], b = dp[i][b];
	  return dp[0][a];
	}

.. _hash_map:

HashMap
============

.. code-block:: cpp

	struct HashMap {
	  int p[M], v[M], f[M], idx; ll a[M]; // ll v[M] if (ll)u
	  void init() { idx = 0, clr(p, 0xff); }
	  void insert(int u, ll t) { //add
	    int x = u % M;
	    for (int i = p[x]; ~i; i = f[i]) {
	      if (v[i] == u) {
	        a[i] += t;
	        return;
	      }
	    }
	    a[idx] = t;
	    v[idx] = u, f[idx] = p[x], p[x] = idx++;
	  }
	};

.. _tree_linear:

TreeLinear
============

.. code-block:: cpp

	int L[N], R[N], _;
	void dfs(int u, int pre) {
	  L[u] = ++_;
	  for (int i = p[u]; ~i; i = e[i].next) {
	    if (e[i].u != pre) dfs(e[i].u, u);
	  }
	  R[u] = _;
	}

.. _binary_search:

BinarySearch
============

.. code-block:: cpp

  // for a[]
  upper_bound(a, a + n, x) - a - 1; // a[res] <= x
  lower_bound(a, a + n, x) - a - 1; // a[res] < x
  lower_bound(a, a + n, x) - a; // a[res] >= x
  upper_bound(a, a + n, x) - a; // a[res] > x

  // for vector
  upper_bound(v.begin(), v.end(), x) - v.begin() - 1; // v[res] <= x
  lower_bound(v.begin(), v.end(), x) - v.begin() - 1; // v[res] <= x
  lower_bound(v.begin(), v.end(), x) - v.begin(); // v[res] <= x
  upper_bound(v.begin(), v.end(), x) - v.begin(); // v[res] <= x

  // for set/multiset
  *s.lower_bound(x); // res >= x
  *s.upper_bound(x); // res > x

.. _discrete:

Discrete
============

.. code-block:: cpp

  vector<int> vt(a, a + n);
  sort(vt.begin(), vt.end());
  vt.erase(unique(vt.begin(), vt.end()), vt.end());

.. _lis:

LIS
============

.. code-block:: cpp

  int lis(int a[], int n) {
    vector<int> dp;
    rep(i, n) {
      int t = lower_bound(dp.begin(), dp.end(), a[i]) - dp.begin();
      if (t >= dp.size()) dp.push_back(a[i]);
      else dp[t] = a[i];
    }
    return (int)dp.size();
  }

.. _matrix:

Matrix
============

.. code-block:: cpp

	//*warning: stackoverflow
	struct Matrix {
	  int n; ll a[N][N];
	  Matrix(int _n = 0) {
	    n = _n;
	    clr(a, 0);
	  }
	  Matrix operator* (Matrix const &t) {
	    Matrix r(n);
	    rep(i, n) rep(j, n) if (a[i][j]) rep(k, n) {
	      r.a[i][k] += a[i][j] * t.a[j][k];
	    }
	    return r;
	  }
	  Matrix operator^ (ll m) {
	    Matrix r(n); rep(i, n) r.a[i][i] = 1;
	    Matrix s(*this);
	    for (; m; m >>= 1) {
	      if (m & 1) r = r * s;
	      s = s * s;
	    }
	    return r;
	  }
	  void pr() {
	    rep(i, n) rep(j, n) cout << a[i][j] << (j == n - 1 ? '\n' : ' ');
	  }
	} ;

.. _cartesian_tree:

CartesianTree
=============

.. code-block:: cpp

	int const N = 5005000;
	int a[N];
	int tr[N][2], st[N], top;
	int cartesian(int n) {
	  top = 0; 
	  Rep(i, n) tr[i][0] = tr[i][1] = 0;
	  Rep(i, n) {
	    int t = top;
	    while (t > 0 && a[st[t]] > a[i]) --t;
	    if (t) tr[st[t]][1] = i;
	    if (t < top) tr[i][0] = st[t + 1];
	    st[++t] = i; 
	    top = t;
	  }
          return st[1];
	}

.. _kd_tree:

KDTree
============

.. code-block:: cpp

	//*warning: coincident points
	inline ll sqr(ll x) { return x * x; }
	int k, cur; //k: Dimension
	struct Point {
	  ll x[M]; //M: max_k
	  bool operator < (Point const &t) const {
	    return x[cur] < t.x[cur];
	  }
	} p[N];

	struct KD_Tree {
	  Point tp[N];
	  int sel; ll ret;
	  void build(int l, int r, int d = 0) {
	    if (l >= r) return;
	    if (d == k) d = 0;
	    int mid = (l + r) >> 1; cur = d;
	/* optimization
	  ll x[M][2];
	  rep(i, k) x[i][0] = Inf, x[i][1] = -Inf;
	  for (int i = l; i < r; ++i) rep(j, k) {
	    x[j][0] = min(x[j][0], p[i].x[j]);
	    x[j][1] = max(x[j][1], p[i].x[j]);
	  }
	  g[mid] = 0;
	  for (int i = 1; i < k; ++i) {
	    if (x[i][1] - x[i][0] > x[g[mid]][1] - x[g[mid]][0]) {
	      g[mid] = i;
	    }
	  }
	  cur = g[mid];
	*/
	    nth_element(p + l, p + mid, p + r);
	    tp[mid] = p[mid];
	    if (l + 1 == r) return;
	    build(l, mid, d + 1);
	    build(mid + 1, r, d + 1);
	  }
	  void query(int l, int r, Point o, int d = 0) {
	    if (l >= r) return;
	    if (d == k) d = 0;
	    int mid = (l + r) >> 1;
	    ll t = 0; rep(i, k) t += sqr(o.x[i] - tp[mid].x[i]);
	    if (t < ret && t) { // && t (ignore itself)
	      ret = t, sel = mid;
	    }
	    int l1 = l, r1 = mid, l2 = mid + 1, r2 = r;
	    if (o.x[d] > tp[mid].x[d]) swap(l1, l2), swap(r1, r2);
	    query(l1, r1, o, d + 1);
	    if (ret > sqr(o.x[d] - tp[mid].x[d])) {
	      query(l2, r2, o, d + 1);
	    }
	  }
	} ;
