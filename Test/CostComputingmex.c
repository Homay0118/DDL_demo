


#include<mex.h>
#include <stdio.h>
#include <math.h>
#include <float.h>


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

	int *dims, M, N, index1, index, numDisp;
    double *Alpha1,*Alpha2, *dispVol,*A2;
	int width,height,i,j,k,u,v,sz,numelements;


    if (nrhs < 1) {
        mexErrMsgTxt("At least one argument is required.") ;
    } else if(nrhs > 6) {
        mexErrMsgTxt("Too many input arguments.");
    }
	
    numelements   = mxGetNumberOfElements(prhs[0]) ;
    dims  = mxGetDimensions(prhs[0]) ;
    Alpha1 = (double*)mxGetData(prhs[0]) ;//mxGetData returns a void pointer, so cast it

    width = dims[1]; height = dims[0];//Note: first dimension provided is height and second is width
   
	Alpha2 = (double*)mxGetData(prhs[1]) ;

	numDisp = mxGetScalar(prhs[2]);
	M = mxGetScalar(prhs[3]);
	N = mxGetScalar(prhs[4]);
	sz = M*N;

	dispVol = mxMalloc(sizeof(double)*M*N*numDisp);

	plhs[0] = mxCreateNumericMatrix(M*N,numDisp,mxDOUBLE_CLASS, mxREAL);
	dispVol = (double *)mxGetData(plhs[0]);//gives a void*, cast it to int*
	
	for (i = 0; i < numDisp*M*N; i++)
		    dispVol[i] = 0;

	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++)
		{
			index = j*M + i;
			for (k = 0; k < numDisp; k++)
			{
				index1 = max(0,j-k-1)*M+i;

				for (u = 0; u < height; u++)

					dispVol[index + k*sz] = dispVol[index + k*sz] + fabs(Alpha1[index*height + u] - Alpha2[index1*height + u]);
			}
		}
}
