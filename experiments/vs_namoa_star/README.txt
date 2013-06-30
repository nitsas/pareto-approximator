This folder contains the code for the experiments we did to compare 
Chord and PGEN with NAMOA*. We used the PGL library's NAMOA* implementation.
We also used PGL's implementation of the A* and (single objective) Dijkstra 
algorithms for our COMB routine. 

The COMB routine is a routine that optimizes (minimizes) linear combinations 
of the objective functions. We give it a set of weights w_{i} and it 
optimizes (minimizes) the combined objective:
\f$ w_{1} f_{1} + w_{2} f_{2} + ... + w_{n} f_{n} \f$,
where f_{i} is the ith objective function. The COMB routine is specific to 
the underlying problem (in this case Multicriteria Shortest Paths). Chord 
and PGEN use the COMB routine as a black-box to generate Pareto optimal 
points (on the convex hull of the Pareto set).

In this folder we use the words "Chord" and "PGEN" interchangeably. Both 
algorithms essentially do the same thing but for different numbers of 
objective functions each (Chord for 2, PGEN for 3+. PGEN can also work 
for 2 objective functions but Chord is simpler).

Hardcoded stuff:
-------------------------------
The option to use 2 or 3 objectives is currently hardcoded. The places where 
changes have to be made are labeled with CHANGE-HERE.
