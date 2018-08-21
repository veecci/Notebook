.. _number_theory:

*************
Number Theory
*************

.. _gcd:

Greatest Common Divisor
=========================

.. code-block:: cpp

  int gcd(int a, int b) {
    if (b == 0) {
      return a;
    } else {
      return gcd(b, a % b);
    }
  }

.. _sievePrimes:

Sieve Primes
=================

.. code-block:: cpp
  
  int pri[N], mark[N], cnt;
  void sieve() {
    mark[0] = mark[1] = 1;
    for (int i = 2; i < N; ++i) {
      if (!mark[i]) mark[i] = pri[cnt++] = i;
      for (int j = 0; pri[j] * i < N; ++j) {
        mark[i * pri[j]] = pri[j];
        if (i % pri[j] == 0) break;
      }
    }
  }
