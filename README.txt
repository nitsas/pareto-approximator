Author:  Christos Nitsas
Date:    2012-2013


Approximate the Pareto set of a multiobjective optimization problem using 
the Chord and PGEN algorithms. (both algorithms based on the weighted 
sum method)

The software presented here is general and can be applied to any 
multicriteria optimization problem, as long as that problem has a COMB 
routine that optimizes (minimizes) linear combinations of the objectives.

Please read [1] for more info on the Chord algorithm. (biobjective case)

Please read [2] for more info on the PGEN algorithm. (3+ objectives)


In-code comments:
------------------------------
In-code comments follow the Doxygen QT style so if you don't understand some
symbol in the comments (comments contain symbols/commands similar to Latex) or 
are just curious about Doxygen please check here:
http://www.stack.nl/~dimitri/doxygen/docblocks.html
and here:
http://www.stack.nl/~dimitri/doxygen/manual.html


Requirements:
------------------------------
1. The Armadillo C++ linear algebra library is required (chiefly) for some 
   linear algebra operations inside the Facet class (and some less important 
   operations elsewhere). The earliest version I've tried is 3.6.1 and it 
   worked. 
   (you can find it on http://arma.sourceforge.net/)
2. The Qhull program (specifically "qconvex") is required by PGEN to compute 
   the convex hull of sets of points. 
   (you can find it on http://www.qhull.org/)


Usage:
------------------------------
Let's say that we have a specific instance S of a multicriteria optimization 
problem P and we want to use Chord or PGEN to approximate its Pareto set. We 
will need to create a class that inherits from BaseProblem, let's call it 
MyProblem, (we can put any problem specific data (e.g. graphs e.t.c.) inside 
MyProblem) and implement MyProblem's comb() function. comb() is declared 
virtual (and null) in BaseProblem and used by the Chord and PGEN algorithms 
to access the problem's Pareto set (i.e. generate Pareto points).

The comb() function is the implementation of the theoretical COMB routine. 
(see below for more info on the COMB routine)

Please check the examples (./examples/) and experiments 
(./experiments/vs_namoa_star/) for examples of how to use the software.


What is the COMB routine:
------------------------------
The COMB routine is a routine that optimizes (minimizes) linear combinations 
of the objective functions. We give it a set of weights w_{i} and it 
optimizes (minimizes) the combined objective:
\f$ w_{1} f_{1} + w_{2} f_{2} + ... + w_{n} f_{n} \f$,
where f_{i} is the ith objective function. The COMB routine is specific to 
the underlying problem (in this case Multicriteria Shortest Paths). Chord 
and PGEN use the COMB routine as a black-box to generate Pareto optimal 
points (on the convex hull of the Pareto set).

We have assumed (w.l.o.g.) that all objectives are supposed to be minimized. 
If some objectives are maximization objectives they can easily be converted 
to equivalent minimization objectives.


References:
------------------------------
[1] C. Daskalakis, I. Diakonikolas and M. Yannakakis: "How good is the Chord 
algorithm?". Proceedings of the Twenty-First Annual ACM-SIAM Symposium on 
Discrete Algorithms. Society for Industrial and Applied Mathematics, 2010.

[2] D. L. Craft, T. F. Halabi, H. A. Shih and T. R. Bortfeld: "Approximating 
convex Pareto surfaces in multiobjective radiotherapy planning". Med. Phys., 
33, pp. 3399-3407, 2006.
