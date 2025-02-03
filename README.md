# How to use EIQP
EQIP solves the convex QP: $\min \frac{1}{2} z^\top Q z + z^\top c, \text{s.t.} Az\geq b, z\geq0$

## How to compile in Matlab:  mex -O EIQP.c -lblas -llapack
