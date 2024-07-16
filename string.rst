.. _string:

*************
String
*************

.. _kmp:

KMP
============

.. code-block:: cpp

	int fail[N], len;
	void buildF(char *p) {
	  for (int i = 1, j = fail[0] = -1; i < len; ++i) {
	    while (~j && p[i] != p[j + 1]) j = fail[j];
	    fail[i] = j += p[i] == p[j + 1];
	  }
	}
	int kmp(char *s, char *p) {
	  len = strlen(p); buildF(p);
	  int ret (0);
	  for (int i = 0, j = -1; s[i]; ++i) {
	    while (~j && s[i] != p[j + 1]) j = fail[j];
	    if (s[i] == p[j + 1]) ++j;
	    if (j == len - 1) {
	      ++ret;
	      j = fail[j];
	    }
	  }
	  return ret;
	}
	// all repetends
	vector<int> r; //empty
	void get_r() {
	  int now = len - 1;
	  while (~now) {
	    r.push_back(len - 1 - fail[now]);
	    now = fail[now];
	  }
	}


.. _ext_kmp:

extKMP
============

.. code-block:: cpp

	// str = "aaaba", pat = "aba", then ex[] = {1, 1, 3, 0, 1}, ext[] = {3, 0, 1}
	int const N = 200200;
	int ext[N], ex[N]; 
	void extkmp(char *str, char *pat, int ext[], int ex[]) {
	  int la = strlen(str), lb = strlen(pat);
	  int p = 0, k = 1;
	  while (pat[p] == pat[p + 1]) ++p;
	  ext[0] = lb, ext[1] = p;
	  for (int i = 2; i < lb; ++i) {
	    int x = k + ext[k] - i, y = ext[i - k];
	    if (y < x) ext[i] = y;
	    else {
	      p = max(0, x);
	      while (pat[p] == pat[p + i]) ++p;
	      ext[i] = p;
	      k = i;
	    }
	  }
	  p = k = 0;
	  while (str[p] && str[p] == pat[p]) ++p;
	  ex[0] = p;
	  for (int i = 1; i < la; ++i) {
	    int x = k + ex[k] - i, y = ext[i - k];
	    if (y < x) ex[i] = y;
	    else {
	      p = max(0, x);
	      while (pat[p] && pat[p] == str[p + i]) ++p;
	      ex[i] = p;
	      k = i;
	    }
	  }
	}

.. _trie:

Trie
============

.. code-block:: cpp

	int const N = 100100; // size
	int const M = 32;     // length
	struct Trie {
	  int idx, cnt;
	  struct Trie_Node {
	    int id;
	    Trie_Node *next[26];
	    void init() {
	      id = -1;
	      clr(next, NULL);
	    }
	  } trie[N*M], root;

	  int insert(char* s) {
	    Trie_Node *p = &root;
	    for (int i = 0; s[i]; ++i) {
	      int j = s[i] - 'a';
	      if (p -> next[j] == NULL) {
	        trie[idx].init();
	        p -> next[j] = &trie[idx++];
	      }
	      p = p -> next[j];
	    }
	    if (p -> id == -1) p -> id = cnt++;
	    return p -> id;
	  }
	  void init() {
	    root.init();
	    idx = cnt = 0;
	  }
	};

.. _minimum_representation:

MinimumRepresentation
======================

.. code-block:: cpp

	int mins(char *s, int n) {
	  int i = 0, j = 1, len = 0, x, y;
	  while (i < n && j < n && len < n) {
	    x = i + len; if (x >= n) x -= n;
	    y = j + len; if (y >= n) y -= n;
	    if (s[x] == s[y]) ++len;
	    else if (s[x] < s[y]) j += len + 1, len = 0;
	    else i = j++, len = 0;
	  }
	  return i;
	}

.. _manacher:

Manacher
============

.. code-block:: cpp

	struct Manacher {
	  int p[N<<1], len; char str[N<<1];
	  int id, ret; // maxPalindrome_idx, maxPalindrome_length
	  void func() {
	    int mx (0);
	    id = 0;
	    rep(i, len) {
	      if (mx > i) p[i] = min(p[id + id - i], mx - i);
	      else p[i] = 1;
	      for (; str[i + p[i]] == str[i - p[i]]; ++p[i]);
	      ret = max(ret, p[i]);
	      if (p[i] + i > mx) {
	        mx = p[i] + i;
	        id = i;
	      }
	    }
	    --ret;
	  }
	  void init(char *s) {
	    // "aaa" -> "!#a#a#a#"
	    len = 0; str[len++]= '!', str[len++] = '#';
	    for (int i = 0; s[i]; ++i) {
	      str[len++] = s[i];
	      str[len++] = '#';
	    }
	    str[len] = 0;
	    ret = 0;
	    func();
	  }
	  bool check(int l, int r) {
	    if ((r - l) % 2 == 0) {
	      int mid = (r + l) / 2;
	      mid = mid * 2 + 2;
	      return p[mid] == r - l + 2;
	    } else {
	      int mid = (l + 1 + r + 1);
	      return p[mid] == r - l + 2;
	    }
	  }
	};

.. _suffix_array:

SuffixArray
============

.. code-block:: cpp

	int const N = 200200;
	int arr[3][N], cnt[N], mc[256], h[N], *sa, *ta, *r, *tr, sz;
	void sa_init(char *str, int len) {
	  sa = arr[0], ta = arr[1], r = arr[2], sz = 0;
	  rep(i, len) ta[i] = str[i];
	  sort(ta, ta + len);
	  Rep(i, len) {
	    if (ta[i] != ta[i - 1] || i == len) {
	      cnt[mc[ta[i - 1]] = sz++] = i;
	    }
	  }
	  for (int i = len - 1; i >= 0; --i) {
	    sa[--cnt[r[i] = mc[str[i]]]] = i;
	  }
	  for (int k = 1; k < len && r[sa[len - 1]] < len - 1; k <<= 1) {
	    for (int i = 0; i < len; i++) {
	      cnt[r[sa[i]]] = i + 1;
	    }
	    for (int i = len - 1; i >= 0; --i) {
	      if (sa[i] >= k) {
	        ta[--cnt[r[sa[i] - k]]] = sa[i] - k;
	      }
	    }
	    for (int i = len - k; i < len; ++i) {
	      ta[--cnt[r[i]]] = i;
	    }
	    tr = sa, sa = ta, tr[sa[0]] = 0;
	    for (int i = 1; i < len; ++i) {
	      tr[sa[i]] =
	        tr[sa[i - 1]] + (r[sa[i]] != r[sa[i - 1]] || sa[i - 1] + k >= len || r[sa[i] + k] != r[sa[i - 1] + k]);
	    }
	    ta = r, r = tr;
	  }
	
	  for (int i = 0, d = 0, j; i < len; ++i) {
	    if (str[i] == '#' || r[i] == len - 1) {
	      h[r[i]] = d = 0;
	    } else {
	      if (d) {
	        --d;
	      }
	      j = sa[r[i] + 1];
	      while (str[i + d] != '#' && str[j + d] != '#' && str[i + d] == str[j + d]) {
	        ++d;
	      }
	      h[r[i]] = d;
	    }
	  }
	}

.. _aho_corasick:

Aho-corasick(trie graph)
=========================

.. code-block:: cpp

	int root, idx;
	struct trie_node{
	    int next[size];
	    int fail;
	    bool flag;
	    void init(){
	        fail = -1, flag = false;
	        memset(next, 0, sizeof(next));
	    }
	}trie[maxn * leng];
	int q[maxn * leng];
	void trie_init(){
	    root = idx = 0;
	    trie[root].init();
	}
	void insert(char *s){
	    int i, j, p = root;
	    for(i=0;s[i];i++){
	        j = s[i] - 'A';
	        if(!trie[p].next[j]){
	            trie[++idx].init();
	            trie[p].next[j] = idx;
	        }
	        p = trie[p].next[j];
	    }
	    trie[p].flag = true;
	}
	void build(){
	    int j, p;
	    q[0] = root;
	    for(int l=0,h=1;l<h;){
	        p = q[l++];
	        for(j=0;j<size;j++){
	            if(trie[p].next[j]){
	                q[h++] = trie[p].next[j];
	                if(trie[p].fail == -1)
	                    trie[trie[p].next[j]].fail = root;
	                else{
	                    trie[trie[p].next[j]].fail =
	                        trie[trie[p].fail].next[j];

	                    trie[trie[p].next[j]].flag |=
	                        trie[trie[trie[p].fail].next[j]].flag;
	                }
	            }
	            else{
	                if(trie[p].fail != -1)
	                    trie[p].next[j] = trie[trie[p].fail].next[j];
	            }
	        }
	    }
	}
