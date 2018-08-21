DSU
---
.. code-block:: cpp

  int parents[N];
  int find(int a) { return parents[a] < 0 ? a : parents[a] = find(parents[a]); }
  void uni(int a, int b) {
    a = find(a), b = find(b); if (a == b) return;
    if (parents[a] < parents[b]) { parents[a] += parents[b], parents[b] = a; }
    else { parents[b] += parents[a], parents[a] = b; }
  }
  void init() { clr(parents, 0xff); }
