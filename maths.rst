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

  //sieve mu
  int pri[N], mark[N], mu[N], cnt;
  void sieve_mu() {
    cnt = 0, mu[1] = mark[0] = mark[1] = 1;
    for (int i = 2; i < N; ++i) {
      if (!mark[i]) pri[cnt++] = mark[i] = i, mu[i] = -1;
      for (int j = 0; pri[j] * i < N; ++j) {
        mark[pri[j] * i] = pri[j];
        if (!(i % pri[j])) {
          mu[pri[j] * i] = 0;
          break;
        }
        mu[pri[j] * i] = -mu[i];
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

.. _square_free:

Squarefree
====================

.. code-block:: cpp

  // Number of (positive) squarefree numbers <= n.
  ll square_free_prefix_sum(ll n) {
     ll m = sqrtl(n), ret = 0;
     for (ll d = 1; d <= m; ++d) {
       ret += mu[d] * (n / (d * d));
     }
    return ret;
  }

.. _basis:

Basis
=================

.. code-block:: cpp

  int const B = 63;
  struct Basis {
    ll d[B]; int cnt = 0; bool zero;
    Basis() { rep(i, B) d[i] = 0; cnt = 0; zero = 0; }
    void rebuild() { for (int i = B - 1; i >= 0; --i) { for (int j = i - 1; j >= 0; --j) { if (d[i] >> j & 1) d[i] ^= d[j]; } } }
    bool ask(ll x) { for (int i = B - 1; i >= 0; --i) { if (x >> i & 1) x ^= d[i]; } return x == 0; }
    ll mx() { ll r = 0; for (int i = B - 1; i >= 0; --i) { if ((r ^ d[i]) > r) r ^= d[i]; } return r; }
    ll mi() { if (zero) return 0; rep(i, B) if (d[i]) { return d[i]; } return 0; }
    void debug() { rep(i, B) printf("%d ", d[i]); puts(""); }
    void add(ll x) {
      bool ins = 0;
      for (int i = B - 1; i >= 0; --i) { 
        if (x >> i & 1) { 
          if (d[i]) x ^= d[i]; 
          else { 
            d[i] = x; ins = 1; break; 
          } 
        } 
      }
      if (ins) ++cnt; else zero = 1;
    }
    ll kth(ll k) {
      if (zero) --k;
      if (k >= (1LL << cnt)) return -1;
      rebuild();
      ll r = 0; int top = 0; 
      rep(i, B) if (d[i]) {
        if (k >> top & 1) r ^= d[i];
        ++top;
      }
      return r;
    }
  } ba;

.. _mint:

mint
=================

.. code-block:: cpp

  int const mod = 998244353;
  struct mint {
    ll v;
    mint(ll x = 0) { v = x % mod; }
    mint& f(ll t) { v = t < 0 ? t + mod : t; return *this; }
    mint operator-() const { return mint() - *this; }
    mint &operator+=(const mint& rhs) { return f(v + rhs.v - mod); }
    mint &operator-=(const mint& rhs) { return f(v - rhs.v); }
    mint &operator*=(const mint& rhs) { v = v * rhs.v % mod; return *this; }
    mint &operator/=(const mint& rhs) { return *this *= rhs.inv(); }
    mint operator+(const mint& rhs) const { return mint(*this) += rhs; }
    mint operator-(const mint& rhs) const { return mint(*this) -= rhs; }
    mint operator*(const mint& rhs) const { return mint(*this) *= rhs; }
    mint operator/(const mint& rhs) const { return mint(*this) /= rhs; }
    friend mint operator+(ll x, const mint& rhs) { return mint(x) + rhs; }
    friend mint operator-(ll x, const mint& rhs) { return mint(x) - rhs; }
    friend mint operator*(ll x, const mint& rhs) { return mint(x) * rhs; }
    friend mint operator/(ll x, const mint& rhs) { return mint(x) / rhs; }
    bool operator<(const mint& rhs) const{ return v < rhs.v;}
    bool operator==(const mint& rhs) const{ return v == rhs.v;}
    bool operator!=(const mint& rhs) const{ return v != rhs.v;}
    operator bool() const { return v; }
    operator int() const { return v; }
    operator ll() const { return v; }
    mint inv() const { return pow(mod - 2); }
    mint pow(ll n) const {
      if (n < 0) return inv().pow(-n);
      mint r(1), x(*this);
      while (n) {
        if (n & 1) r *= x;
        x *= x;
        n >>= 1;
      }
      return r;
    }
    friend ostream& operator<<(ostream&os, const mint&t) { return os << t.v; }
    friend istream& operator>>(istream&is, mint&t) { ll x; is >> x; t = mint(x); return is; }
  };

.. _n_i enumerate:

n/i Enumerate
=================

.. code-block:: cpp

  for (int i = 1, j; i <= n; i = j + 1) {
    j = n / (n / i);
    // n / i : [i, j]
  }

.. _comb_mod:

Combination(mod)
=================

.. code-block:: cpp

  ll fac[N], inv[N];
  ll C(int n, int m) {
    if (n < m) return 0;
    return fac[n] * inv[m] % mod * inv[n - m] % mod;
  }
  void Cinit() {
    fac[0] = inv[0] = inv[1] = 1;
    for (int i = 1; i < N; ++i) fac[i] = fac[i - 1] * i % mod;
    for (int i = 2; i < N; ++i) inv[i] = inv[mod % i] * (mod - mod / i) % mod;
    for (int i = 2; i < N; ++i) inv[i] = inv[i - 1] * inv[i] % mod;
  }

.. _lucas

Lucas
=======

.. code-block:: cpp

  ll lucas(ll n, ll m) {
    if (m == 0) return 1;
    return (C(n % mod, m % mod) * lucas(n / mod, m / mod)) % mod;
  }

.. _lagrange_interpolation

Lagrange Interpolation
============================

.. code-block:: cpp

  // intern
  // a = {f(0), f(1), ... , f(n)}
  int const N = 10100;
  ll fac[N], inv[N];
  ll lagrange(ll a[N], int n, ll x, ll mod) {
    fac[0] = inv[0] = inv[1] = 1;
    for (int i = 1; i <= n; ++i) fac[i] = fac[i - 1] * i % mod;
    for (int i = 2; i <= n; ++i) inv[i] = inv[mod % i] * (mod - mod / i) % mod;
    for (int i = 2; i <= n; ++i) inv[i] = inv[i - 1] * inv[i] % mod;
    vector<ll> pre(n + 1), suf(n + 1); 
    x %= mod; pre[0] = x; Rep(i, n) pre[i] = pre[i - 1] * (x - i) % mod;
    suf[n] = (x - n) % mod; for (int i = n - 1; i >= 0; --i) suf[i] = suf[i + 1] * (x - i) % mod;
    ll ret = 0;
    rep(i, n + 1) {
      ll di = (i == 0 ? 1LL : pre[i - 1]) * (i == n ? 1LL : suf[i + 1]) % mod;
      ll t = di * inv[i] % mod * inv[n - i] % mod;
      if ((n - i) & 1) ret -= t * a[i];
      else ret += t * a[i]; 
      ret %= mod;
    }
    if (ret < 0) ret += mod;
    return ret;
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

  // intern
  // S2(n, m) = (1/m!) * Sum_{i=0..k} (-1)^(m-i)*binomial(m, i)*i^n.
  ll stirling2(int n, int m) {
    ll sum = 0;
    Rep(i, m) {
      ll t = C(m, i) * powMod(i, n, mod) % mod;
      if ((m - i) & 1) sum -= t;
      else sum += t;
      sum %= mod;
    }
    mint ret = sum; 
    ret /= fac[m];
    return ret;
  }


.. _fft:

FFT
=================

.. code-block:: cpp

  int const N = 100100;
  double const pi = atan2(0, -1);
  struct E {
    double x, y;
    E(double x = 0, double y = 0) : x(x), y(y) {}
    E operator-(const E &b) const { return E(x - b.x, y - b.y); }
    E operator+(const E &b) const { return E(x + b.x, y + b.y); }
    E operator*(const E &b) const { return E(x * b.x - y * b.y, x * b.y + y * b.x); }
  };

  struct FFT {
    int n, m, l, r[N*2]; ll s[N][2]; int re[N*2];
    E a[N*2], b[N*2];
    void fft(E *a, int sig) {
      rep(i, n) if (i < r[i]) swap(a[i], a[r[i]]);
      for (int i = 1; i < n; i <<= 1) {
        E wn(cos(pi / i), sig * sin(pi / i));
        for (int j = 0, p = i << 1; j < n; j += p) {
          E w(1, 0);
          for (int k = 0; k < i; ++k, w = w * wn) {
            E x = a[j + k], y = w * a[j + k + i];
            a[j + k] = x + y; a[j + k + i] = x - y;
          }
        }
      }
    }
    int solve(int _n, int _m) {
      n = _n - 1, m = _m - 1, l = 0;
      rep(i, n + 1) a[i] = E(s[i][0], 0);
      rep(i, m + 1) b[i] = E(s[i][1], 0);
      m += n; for (n = 1; n <= m; n <<= 1) ++l;
      rep(i, n) r[i] = (r[i >> 1] >> 1) | ((i & 1) << (l - 1));
      fft(a, 1); fft(b, 1);
      rep(i, n + 1) a[i] = a[i] * b[i];
      fft(a, -1);
      rep(i, m + 1) re[i] = (int)(a[i].x / n + 0.5);
      return m + 1;
    }
    void init() {
      clr(a, 0); clr(b, 0); clr(re, 0);
    }
  } fft;

.. _fwt:

FWT
=================
:math:`C_i = \sum_{j\oplus k=i}A_j\times B_k`

.. code-block:: cpp

  int const M = 17;
  int const mod = 998244353;
  int const inv2 = 499122177; // inv(2, 1, mod)
  ll a[1<<M], b[1<<M], ta[1<<M], tb[1<<M];

  void fwt_or(ll a[], ll n, int t) {
    for (int k = 1; k < n; k <<= 1) {
      for (int i = 0; i < n; i += (k << 1)) {
        rep(j, k) {
          a[i + j + k] = (a[i + j + k] + a[i + j] * t + mod) % mod;
        }
      }
    }
  }

  void fwt_xor(ll a[], ll n, int t) {
    for (int k = 1; k < n; k <<= 1) {
      for (int i = 0; i < n; i += (k << 1)) {
        rep(j, k) {
          ll x = a[i + j], y = a[i + j + k];
          a[i + j] = (x + y) * (t == 1 ? 1 : inv2) % mod;
          a[i + j + k] = (x - y + mod) * (t == 1 ? 1 : inv2) % mod;
        }
      }
    }
  }

  void fwt_and(ll a[], ll n, int t) {
    for (int k = 1; k < n; k <<= 1) {
      for (int i = 0; i < n; i += (k << 1)) {
        rep(j, k) {
          a[i + j] = (a[i + j] + a[i + j + k] * t + mod) % mod;
        }
      }
    }
  }

  void solve(void (*fwt)(ll *a, ll n, int t), ll a[], ll b[], int n) {
    rep(i, n) ta[i] = a[i];
    rep(i, n) tb[i] = b[i];
    fwt(ta, n, 1), fwt(tb, n, 1);
    rep(i, n) ta[i] = (ta[i] * tb[i]) % mod;
    fwt(ta, n, -1);
  }

.. _sum_specific_bit_count:

:math:`\sum_{i=1}^{n}(i\gg k\space\&\space 1)`

======================================================

.. code-block:: cpp

  ll sum_bit(ll n, int k) {
    ll x = (n + 1) / (1LL << (k + 1)) * (1LL << k);
    ll y = (n + 1) % (1LL << (k + 1)) - (1LL << k);
    return x + max(y, 0LL);
  }

.. _omega_d_table:

ω(n) and d(n) Table
======================================================

+------------------------+------------------------+-------------------+
|n ≤                     |:math:`max\{\omega(n)\}`|:math:`max\{d(n)\}`|
+========================+========================+===================+
|:math:`10^1`            | 2                      |4                  |
+------------------------+------------------------+-------------------+
|:math:`10^2`            | 3                      |12                 |
+------------------------+------------------------+-------------------+
|:math:`10^3`            | 4                      |32                 |
+------------------------+------------------------+-------------------+
|:math:`10^4`            | 5                      |64                 |
+------------------------+------------------------+-------------------+
|:math:`10^5`            | 6                      |128                |
+------------------------+------------------------+-------------------+
|:math:`10^6`            | 7                      |240                |
+------------------------+------------------------+-------------------+
|:math:`10^7`            | 8                      |448                |
+------------------------+------------------------+-------------------+
|:math:`10^8`            | 8                      |768                |
+------------------------+------------------------+-------------------+
|:math:`10^9`            | 9                      |1344               |
+------------------------+------------------------+-------------------+
|:math:`10^{10}`         | 10                     |2304               |
+------------------------+------------------------+-------------------+
|:math:`10^{11}`         | 10                     |4032               |
+------------------------+------------------------+-------------------+
|:math:`10^{12}`         | 11                     |6720               |
+------------------------+------------------------+-------------------+
|:math:`10^{13}`         | 12                     |10752              |
+------------------------+------------------------+-------------------+
|:math:`10^{14}`         | 12                     |17280              |
+------------------------+------------------------+-------------------+
|:math:`10^{15}`         | 13                     |26880              |
+------------------------+------------------------+-------------------+
|:math:`10^{16}`         | 13                     |41472              |
+------------------------+------------------------+-------------------+
|:math:`10^{17}`         | 14                     |64512              |
+------------------------+------------------------+-------------------+
|:math:`10^{18}`         | 15                     |103680             |
+------------------------+------------------------+-------------------+
