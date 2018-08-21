.. _maths:

******
Maths
******

.. _gcdnlcm:

Gcd & Lcm
=========================

.. code-block:: cpp

  int gcd(int a, int b) { return b ? gcd(b, a % b) : a;}

  int lcm(int a, int b) { return a / gcd(a, b) * b; }

  db fgcd(db a, db b) {
    if(b > -eps && b < eps) return a;
    else return fgcd( b, fmod(a, b) );
  }

.. _powMod:

PowMod & MulMod
=================

.. code-block:: cpp
  
  //powMod
  ll powMod(ll a, ll b, ll c) {
    ll ret = 1 % c;
    for (; b; a = a * a % c, b >>= 1)
      if(b & 1) ret = ret * a % c;
    return ret;
  }

  //powMod plus
  ll mulMod(ll a, ll b, ll c) {
    ll ret = 0;
    for (; b; a = (a << 1) % c, b >>= 1)
      if (b & 1) ret = (ret + a) % c;
    return ret;
  }

  ll powMod(ll a, ll b, ll c) {
    ll ret = 1 % c;
    for (; b; a = mulMod(a, a, c), b >>= 1)
      if (b & 1) ret = mulMod(ret, a, c);
    return ret;
  }

  // mulMod plus for a, b, c <= 2^40
  ll mulMod(ll a, ll b, ll c) {
    return (((a*(b>>20)%c)<<20) + a*(b&((1<<20)-1)))%c;
  }

.. _ext_gcd:

ExtGcd
=================

.. code-block:: cpp

  ll ext_gcd(ll a, ll b, ll &x, ll &y) {
    ll t, ret;
    if (!b) {
      x = 1, y = 0;
      return a;
    }
    ret = ext_gcd(b, a % b, x, y);
    t = x, x = y, y = t - a / b * y;
    return ret;
  }

.. _inv:

Inv
=================

.. code-block:: cpp

  // (b / a) % c
  ll inv(ll a, ll b, ll c) {
    ll x, y;
    ext_gcd(a, c, x, y);
    return (1LL * x * b % c + c) % c;
  }

  //sieve inv
  int inv[N];
  void sieve_inv() {
    inv[1] = 1;
    for (int i = 2; i < N; ++i)
      inv[i] = inv[mod % i] * (mod - mod / i) % mod;
  }

.. _sieve_primes:

Sieve Primes
=================

.. code-block:: cpp

  //mark[i]: the minimum factor of i (for prime, mark[i] = i)
  int pri[N], mark[N], cnt;
  void sieve() {
    cnt = 0, mark[0] = mark[1] = 1;
    for (int i = 2; i < N; ++i) {
      if (!mark[i]) pri[cnt++] = mark[i] = i;
      for (int j = 0; pri[j] * i < N; ++j) {
        mark[i * pri[j]] = pri[j];
        if (!(i % pri[j])) break;
      }
    }
  }

.. _phi:

Phi
=================

.. code-block:: cpp

  //phi
  int phi(int n) {
    int ret = n;
    for (int i = 2; i * i <= n; i += (i != 2) + 1) {
      if (!(n % i)) {
        ret = ret / i * (i - 1);
        while (n % i == 0) n /= i;
      }
    }
    if (n > 1) ret = ret / n * (n - 1);
    return ret;
  }

  //phi plus (sieve() first & (n < N))
  int phi(int n) {
    int ret = n, t;
    while ((t = mark[n]) != 1) {
      ret = ret / t * (t - 1);
      while (mark[n] == t) n /= mark[n];
    }
    return ret;
  }

  //sieve phi
  int pri[N], phi[N], cnt;
  void sieve_phi() {
    cnt = 0, phi[1] = 1;
    for (int i = 2; i < N; ++i) {
      if (!phi[i]) pri[cnt++] = i, phi[i] = i - 1;
      for (int j = 0; pri[j] * i < N; ++j) {
        if (!(i % pri[j])) {
          phi[i * pri[j]] = phi[i] * pri[j];
          break;
        } else {
          phi[i * pri[j]] = phi[i] * (pri[j] - 1);
        }
      }
    }
  }

.. _divisors:

Divisors
=================

.. code-block:: cpp

  //the number of divisors
  int d_func(int n) {
    int ret = 1, t = 1;
    for (int i = 2; i * i <= n; i += (i != 2) + 1) {
      if (!(n % i)) {
        while (!(n % i)) ++t, n /= i;
        ret *= t, t = 1;
      }
    }
    return n > 1 ? ret << 1 : ret;
  }

  //sieve the number of divisors (O(nlogn))
  int nod[N];
  void sieve_nod() {
    for (int i = 1; i < N; ++i)
      for (int j = i; j < N; j += i)
        ++nod[j];
  }

  //sieve the number of divisors (O(n))
  int pri[N], e[N], divs[N], cnt;
  void sieve_nod() {
    cnt = 0, divs[0] = divs[1] = 1;
    for (int i = 2; i < N; ++i) {
      if (!divs[i]) divs[i] = 2, e[i] = 1, pri[cnt++] = i;
      for (int j = 0; i * pri[j] < N; ++j) {
        int k = i * pri[j];
        if (i % pri[j] == 0) {
          e[k] = e[i] + 1;
          divs[k] = divs[i] / (e[i] + 1) * (e[i] + 2);
          break;
        } else {
          e[k] = 1, divs[k] = divs[i] << 1;
        }
      }
    }
  }

  //the sum of all divisors
  int ds_func(int n) {
    int ret = 1, t;
    for (int i = 2; i * i <= n; i += (i != 2) + 1) {
      if (!(n % i)) {
        t = i * i, n /= i;
        while (!(n % i)) t *= i, n /= i;
        ret *= (t - 1) / (i - 1);
      }
    }
    return n > 1 ? ret * (n + 1) : ret;
  }

.. _miller_rabin:

Miller-Rabin
=================

.. code-block:: cpp

  bool suspect(ll a, int s, ll d, ll n) {
    ll x = powMod(a, d, n);
    if (x == 1) return  true;
    for (int r = 0; r < s; ++r) {
      if (x == n - 1) return  true;
      x = mulMod(x, x, n);
    }
    return false;
  }

  // {2,7,61,-1} is for n < 4759123141 (2^32)
  int const test[] = {2,3,5,7,11,13,17,19,23,-1}; // for n < 10^16
  bool isPrime(ll n) {
    if (n <= 1 || (n > 2 && n % 2 == 0)) return false;
    ll d = n - 1, s = 0;
    while (d % 2 == 0) ++s, d /= 2;
    for (int i = 0; test[i] < n && ~test[i]; ++i)
      if (!suspect(test[i], s, d, n)) return false;
    return true;
  }

.. _pollard_rho:

Pollard-Rho
=================

.. code-block:: cpp

  ll pollard_rho(ll n, ll c) {
    ll d, x = rand() % n, y = x;
    for (ll i = 1, k = 2; ; ++i) {
      x = (mulMod(x, x, n) + c) % n;
      d = gcd(y - x, n);
      if (d > 1 && d < n) return d;
      if (x == y) return n;
      if (i == k) y = x, k <<= 1;
    }
    return 0;
  }

.. _find_factors:

Find Factors
=================

.. code-block:: cpp

  //find factors
  int facs[N];
  int find_fac(int n) {
    int cnt = 0;
    for(int i = 2; i * i <= n; i += (i != 2) + 1)
      while (!(n % i)) n /= i, facs[cnt++] = i;
    if (n > 1) facs[cnt++] = n;
    return cnt;
  }

  //find factors plus (sieve() first & (n < N))
  int facs[N];
  int find_fac(int n) {
    int cnt = 0;
    while (mark[n] != 1)
      facs[cnt++] = mark[n], n /= mark[n];
    return cnt;
  }

.. _place_n_balls_into_m_boxes:

Place n Balls into m Boxes
============================

+-----------+-----------+-------------+----------------------------------------------------------------------------------------------------------------+
|  Balls    |  Boxes    | Empty Boxes |       Answer                                                                                                   |
+===========+===========+=============+================================================================================================================+
| Different | Different |     Yes     | :math:`m^n`                                                                                                    |
+-----------+-----------+-------------+----------------------------------------------------------------------------------------------------------------+
| Different | Different |     No      | :math:`m!S\left( {n,m} \right)`                                                                                |
+-----------+-----------+-------------+----------------------------------------------------------------------------------------------------------------+
| Different | Same      |     Yes     | :math:`S\left( {n,1} \right) + S\left( {n,2} \right) + \ldots + S\left( {n,\min \left( {n,m} \right)} \right)` |
+-----------+-----------+-------------+----------------------------------------------------------------------------------------------------------------+
| Different | Same      |     No      | :math:`S\left( {n,m} \right)`                                                                                  |                                        
+-----------+-----------+-------------+----------------------------------------------------------------------------------------------------------------+
| Same      | Different |     Yes     | :math:`C\left( {n + m - 1,n} \right)`                                                                          |
+-----------+-----------+-------------+----------------------------------------------------------------------------------------------------------------+
| Same      | Different |     No      | :math:`C\left( {n - 1,m - 1} \right)`                                                                          |
+-----------+-----------+-------------+----------------------------------------------------------------------------------------------------------------+
| Same      | Same      |     Yes     | :math:`F\left( {n,m} \right)`                                                                                  |
+-----------+-----------+-------------+----------------------------------------------------------------------------------------------------------------+
| Same      | Same      |     No      | :math:`F\left( {n - m,m} \right)`                                                                              |
+-----------+-----------+-------------+----------------------------------------------------------------------------------------------------------------+

.. code-block:: cpp

  //+ mod if needed
  ll C[N][N];
  void Cinit() {
    for (int i = 0; i < N; ++i) {
      C[i][0] = 1;
      for (int j = 1; j <= i; ++j) {
        C[i][j] = C[i - 1][j] + C[i - 1][j - 1];
      }
    }
  }

  ll S[N][N]; // Strling2[]
  void Sinit() {
    S[0][0] = 1;
    for (int i = 1; i < N; ++i) {
      S[i][1] = 1;
      for (int j = 2; j <= i; ++j) {
        S[i][j] = S[i - 1][j - 1] + j * S[i - 1][j];
      }
    }
  }

  ll F[N][N];
  void Finit() {
    for (int i = 0; i < N; ++i) F[i][1] = F[0][i] = 1;
    for (int i = 1; i < N; ++i) {
      for (int j = 2; j < N; ++j) {
        F[i][j] = F[i][j - 1];
        if (i >= j) F[i][j] += F[i - j][j];
      }
    }
  }