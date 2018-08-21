.. notebook documentation master file, created by
   sphinx-quickstart on Tue Aug 21 12:07:37 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


NOXE
=====================================
 
DSU
------------
.. code-block:: cpp
    
  int parents[N];
  int find(int a) { return parents[a] < 0 ? a : parents[a] = find(parents[a]); }
  void uni(int a, int b) {
    a = find(a), b = find(b); if (a == b) return;
    if (parents[a] < parents[b]) { parents[a] += parents[b], parents[b] = a; }
    else { parents[b] += parents[a], parents[a] = b; }
  }
  void init() { clr(parents, 0xff); }


``rsimpson(0, r1, 1e-4, simpson(0, r1)``

.. note::
     ``rsimpson(0, r1, 1e-4, simpson(0, r1)``
    
:Date: 2001-08-16
:Version: 1
:Authors: - Me
          - Myself
          - I
:Indentation: Since the field marker may be quite long, the second
   and subsequent lines of the field body do not have to line up
   with the first line, but they must be indented relative to the
   field name marker, and they must line up with each other.
:Parameter i: integer

=====  =====  =======
A      B      A and B
=====  =====  =======
False  False  False
True   False  False
False  True   False
True   True   True
=====  =====  =======

This is a paragraph that contains `a link`_.

.. _a link: https://domain.invalid/

=================
This is a heading
=================

.. cpp:function:: bool myMethod(int, double)

   A function with unnamed parameters.

.. cpp:function:: const T &MyClass::operator[](std::size_t i) const

   An overload for the indexing operator.

.. cpp:function:: operator bool() const

   A casting operator.

.. cpp:function:: constexpr void foo(std::string &bar[2]) noexcept

   A constexpr function.

.. cpp:function:: MyClass::MyClass(const MyClass&) = default

   A copy constructor with default implementation.
        
 
API文档
-----------------
 
.. toctree::
   :maxdepth: 2
 
 
Indices and tables
==================
 
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
