# How to use EIQP
EQIP solve the convex QP: min 1/2 z^\top Q z + z^\top c, s.t. Az>=b, z>=0

## How to compile in Matlab:  mex -O EIQP.c -lblas -llapack
