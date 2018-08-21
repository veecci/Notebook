.. _maths:

******
Maths
******

.. _gcdnlcm:

GCD & LCM
=========================

.. code-block:: cpp

  int gcd(int a, int b) { return b ? gcd(b, a % b) : a;}
  int lcm(int a, int b) { return a / gcd(a, b) * b; }
  db fgcd(db a, db b) {
    if(b > -eps && b < eps) return a;
    else return fgcd( b, fmod(a, b) );
  }

.. _powMod:

PowMod
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