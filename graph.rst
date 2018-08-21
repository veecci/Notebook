Graph
-----
.. code-block:: cpp

  int const N = 100010;
  int const M = 1000010;
  struct Edges {
    int u, next;
  } e[M];
  int p[N], idx;
  void init() { clr(p, 0xff); idx = 0; }
  void addedge(int u, int v) { e[idx].u = v, e[idx].next = p[u], p[u] = idx++; }
