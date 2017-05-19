

#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <float.h>


void Path1(double *dispVol,double *I,double *AA,double *BB,double p1,double sigma,int width,int height,int sz,int numDisp,double *S)
{
	double temp, temp_r, temp_g, temp_b, *P2, d_prev, Tmin, Tmin1, Tmin2, Tmin3;
	double *aggr;
	double c = 0.05;
	int i,j,d,d1;
	P2 = mxMalloc(sizeof(double)*(width-1)*height);
	aggr = mxMalloc(sizeof(double)*height*width*numDisp);
	Tmin1 = 1000;Tmin2 = 1000;Tmin3 = 1000;
	for(i= 0 ;i<height;i++)
		for(j= 1;j<width;j++)
		{
			temp_r = fabs(I[(j-1)*height+i] - I[j*height+i]);
			temp_g = fabs(I[(j-1)*height+i+sz] - I[j*height+i+sz]);
			temp_b = fabs(I[(j-1)*height+i+2*sz] - I[j*height+i+2*sz]);
			temp = max(max(temp_r,temp_g),temp_b);
			temp = max(temp,c);
			temp = max(p1,sigma/temp);
			P2[(j-1)*height+i] = temp;
		}
	for(i=0;i<height;i++)
		for(d=0;d<numDisp;d++)
		{
			aggr[i+d*sz] = dispVol[i+d*sz];
			S[i + d*sz] = 0;
		}

for(j=1;j<width;j++)
	for (i = 0; i<height; i++)
		{	
			 d_prev = aggr[(j-1)*height+i+0*sz];
			 for(d=0;d<numDisp;d++)
				{
				d_prev = min(d_prev,aggr[(j-1)*height+i+d*sz]);
				}
				for(d = 0;d<numDisp;d++)
				{
					for(d1=0;d1<numDisp;d1++)
					{
						if(AA[d+d1*numDisp]==1)
						{
								Tmin1 = min(Tmin1,aggr[(j-1)*height+i+d1*sz]);
						}
						if(BB[d+d1*numDisp]==1)
						{
								Tmin2 = min(Tmin2,aggr[(j-1)*height+i+d1*sz]);
						}
					}
				Tmin1 = Tmin1+P2[(j-1)*height+i];
				Tmin2 = Tmin2+p1;
				Tmin3 = aggr[(j-1)*height+i+d*sz];
				Tmin = min(min(Tmin1,Tmin2),Tmin3);
				aggr[j*height+i+d*sz] = dispVol[j*height+i+d*sz]+Tmin-d_prev;
				S[j*height + i + d*sz] = aggr[j*height + i + d*sz];
				Tmin1 = 1000;Tmin2 = 1000;Tmin3 = 1000;
				}
		}
}

void Path2(double *dispVol,double *I,double *AA,double *BB,double p1,double sigma,int width,int height,int sz,int numDisp,double *S)
{
	double temp,temp_r,temp_g,temp_b,*P2,*temp_aggr,d_prev,Tmin,Tmin1,Tmin2,Tmin3;
	double c = 0.05;
	int i,j,d,d1;
	double *aggr;
	P2 = mxMalloc(sizeof(double)*(width-1)*height);
	aggr = mxMalloc(sizeof(double)*height*width*numDisp);
	Tmin1 = 1000;Tmin2 = 1000;Tmin3 = 1000;
	for(i= 0;i<height;i++)
		for(j=1;j<width;j++)
		{
			temp_r = fabs(I[(j-1)*height+i] - I[j*height+i]);
			temp_g = fabs(I[(j-1)*height+i+sz] - I[j*height+i+sz]);
			temp_b = fabs(I[(j-1)*height+i+2*sz] - I[j*height+i+2*sz]);
			temp = max(max(temp_r,temp_g),temp_b);
			temp = max(temp,c);
			temp = max(p1,sigma/temp);
			P2[(j-1)*height+i] = temp;
		}
	for(i=0;i<height;i++)
		for(d=0;d<numDisp;d++)
		{
			aggr[(width-1)*height+i+d*sz] = dispVol[(width-1)*height+i+d*sz];
			S[(width - 1)*height + i + d*sz] = 0;
		}
   for(j=width-2;j>=0;j--)
	for (i = 0; i<height; i++)
		{	
			d_prev = aggr[(j+1)*height+i+0*sz];
			for(d=0;d<numDisp;d++)
			{
				d_prev = min(d_prev,aggr[(j+1)*height+i+d*sz]);
			}
				for(d = 0;d<numDisp;d++)
				{
					for(d1=0;d1<numDisp;d1++)
					{
						if(AA[d+d1*numDisp]==1)
							{
						      Tmin1 = min(Tmin1,aggr[(j+1)*height+i+d1*sz]);
							}
					    if(BB[d+d1*numDisp]==1)
							{
							 Tmin2 = min(Tmin2,aggr[(j+1)*height+i+d1*sz]);
							}	
					}
				Tmin1 = Tmin1+P2[j*height+i];
				Tmin2 = Tmin2+p1;
				Tmin3 = aggr[(j+1)*height+i+d*sz];
				Tmin = min(min(Tmin1,Tmin2),Tmin3);
				aggr[j*height+i+d*sz] = dispVol[j*height+i+d*sz]+Tmin-d_prev;
				S[j*height + i + d*sz] = aggr[j*height + i + d*sz];
				Tmin1 = 1000;Tmin2 = 1000;Tmin3 = 1000;
				}
		}
}


void Path3(double *dispVol,double *I,double *AA,double *BB,double p1,double sigma,int width,int height,int sz,int numDisp,double *S)
{
	double temp,temp_r,temp_g,temp_b,*P2,*temp_aggr,d_prev,Tmin,Tmin1,Tmin2,Tmin3;
	double c = 0.05;
	int i,j,d,d1;
	double *aggr;
	P2 = mxMalloc(sizeof(double)*(height-1)*width);
	aggr = mxMalloc(sizeof(double)*height*width*numDisp);
	Tmin1 = 1000;Tmin2 = 1000;Tmin3 = 1000;

	for(i= 1;i<height;i++)
		for(j= 0;j<width;j++)
		{
			temp_r = fabs(I[j*height + i - 1] - I[j*height + i]);
			temp_g = fabs(I[j*height + i - 1 + sz] - I[j*height + i + sz]);
			temp_b = fabs(I[j*height + i - 1 + 2 * sz] - I[j*height + i + 2 * sz]);
			temp = max(max(temp_r,temp_g),temp_b);
			temp = max(temp,c);
			temp = max(p1,sigma/temp);
			P2[j*(height-1)+i-1] = temp;
		}
	for(j=0;j<width;j++)
		for(d=0;d<numDisp;d++)
		{
			aggr[j*height+d*sz] = dispVol[j*height+d*sz];
			S[j*height + d*sz] = 0;
		}
	
	for(i=1;i<height;i++)
		for (j = 0; j<width; j++)
	      {	
			d_prev = aggr[j*height+i-1+0*sz];
			for(d=0;d<numDisp;d++)
			{
				d_prev = min(d_prev,aggr[j*height+i-1+d*sz]);
			}
				for(d = 0;d<numDisp;d++)
				{
					for(d1=0;d1<numDisp;d1++)
					{
						if(AA[d+d1*numDisp]==1)
							 {
								Tmin1 = min(Tmin1,aggr[j*height+i-1+d1*sz]);
							 }
						if(BB[d+d1*numDisp]==1)
							{
								Tmin2 = min(Tmin2,aggr[j*height+i-1+d1*sz]);
							}	
					}
				Tmin1 = Tmin1+P2[j*(height-1)+i-1];
				Tmin2 = Tmin2+p1;
				Tmin3 = aggr[j*height+i-1+d*sz];
				Tmin = min(min(Tmin1,Tmin2),Tmin3);
				aggr[j*height+i+d*sz] = dispVol[j*height+i+d*sz]+Tmin-d_prev;
				S[j*height + i + d*sz] = aggr[j*height + i + d*sz];
				Tmin1 = 1000;Tmin2 = 1000;Tmin3 = 1000;
				}
	}
}


void Path4(double *dispVol,double *I,double *AA,double *BB,double p1,double sigma,int width,int height,int sz,int numDisp,double *S)
{
	double temp,temp_r,temp_g,temp_b,*P2,*temp_aggr,d_prev,Tmin,Tmin1,Tmin2,Tmin3;
	double c = 0.05;
	int i,j,d,d1;
	double *aggr;
	P2 = mxMalloc(sizeof(double)*(height-1)*width);
	aggr = mxMalloc(sizeof(double)*height*width*numDisp);
	Tmin1 = 1000;Tmin2 = 1000;Tmin3 = 1000;
	for(i= 1 ;i<height;i++)
		for(j= 0;j<width;j++)
		{
			temp_r = fabs(I[j*height+i-1] - I[j*height+i]);
			temp_g = fabs(I[j*height + i - 1 + sz] - I[j*height + i + sz]);
			temp_b = fabs(I[j*height + i - 1 + 2 * sz] - I[j*height + i + 2 * sz]);
			temp = max(max(temp_r,temp_g),temp_b);
			temp = max(temp,c);
			temp = max(p1,sigma/temp);
			P2[j*(height-1)+i-1] = temp;
		}

	for(j=0;j<width;j++)
		for(d=0;d<numDisp;d++)
		{
			aggr[j*height+height-1+d*sz] = dispVol[j*height+height-1+d*sz];
			S[j*height + height - 1 + d*sz] = 0;
		}

   for(i=height-2;i>=0;i--)
		 for (j = 0; j<width; j++)
		{	
			d_prev = aggr[j*height+i+1+0*sz];
			for(d=0;d<numDisp;d++)
			{
				d_prev = min(d_prev,aggr[j*height+i+1+d*sz]);
			}
			for(d = 0;d<numDisp;d++)
			{
				for(d1=0;d1<numDisp;d1++)
				{
					if(AA[d+d1*numDisp]==1)
					{
							Tmin1 = min(Tmin1,aggr[j*height+i+1+d1*sz]);
					}
					if(BB[d+d1*numDisp]==1)
					{
							Tmin2 = min(Tmin2,aggr[j*height+i+1+d1*sz]);
					}	
				}
				Tmin1 = Tmin1+P2[j*(height-1)+i];
				Tmin2 = Tmin2+p1;
				Tmin3 = aggr[j*height+i+1+d*sz];
				Tmin = min(min(Tmin1,Tmin2),Tmin3);
				aggr[j*height+i+d*sz] = dispVol[j*height+i+d*sz]+Tmin-d_prev;
				S[j*height + i + d*sz] = aggr[j*height + i + d*sz];
				Tmin1 = 1000;Tmin2 = 1000;Tmin3 = 1000;
			}
		}
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    int *dims;
    double *dispVol,*I;
	int width,height,i,j,numDisp,sz,numelements;
	double p1,sigma;
	double *DD,*AA,*BB,*S,*S1,*S2,*S3,*S4;
	double *P2;
    if (nrhs < 1) {
        mexErrMsgTxt("At least one argument is required.") ;
    } else if(nrhs > 5) {
        mexErrMsgTxt("Too many input arguments.");
    }


    numelements   = mxGetNumberOfElements(prhs[0]) ;
    dims  = mxGetDimensions(prhs[0]) ;
    dispVol = (double*)mxGetData(prhs[0]) ;//mxGetData returns a void pointer, so cast it
    width = dims[1]; height = dims[0];//Note: first dimension provided is height and second is width
    sz = width*height;
	numDisp = numelements/sz;

	I = (double*)mxGetData(prhs[1]) ;
    //---------------------------

    p1  =  mxGetScalar(prhs[2]);
    sigma  =  mxGetScalar(prhs[3]);
    
	DD = mxMalloc(sizeof(double)*numDisp*numDisp);
	AA = mxMalloc(sizeof(double)*numDisp*numDisp);
	BB = mxMalloc(sizeof(double)*numDisp*numDisp);
	S =  mxMalloc(sizeof(double)*width*height*numDisp);
	S1 = mxMalloc(sizeof(double)*width*height*numDisp);
	S2 = mxMalloc(sizeof(double)*width*height*numDisp);
	S3 = mxMalloc(sizeof(double)*width*height*numDisp);
	S4 = mxMalloc(sizeof(double)*width*height*numDisp);
	P2 = mxMalloc(sizeof(double)*height*(width-1));

	//for (i = 0; i < width*height*numDisp; i++)
	//	S[i] = 0;
	//S1 = S;
	//S2 = S;
	//S3 = S;
	//S4 = S;

	for(i=0;i<numDisp;i++)
		for(j=0;j<numDisp;j++)
		{
			DD[j*numDisp+i] = fabs(j-i);
			AA[j*numDisp+i] = 0;
			BB[j*numDisp+i] = 0;
		}
	for(i=0;i<numDisp;i++)
		for(j=0;j<numDisp;j++)
		{
			if(DD[j*numDisp+i]>1)
				AA[j*numDisp+i]=1;
			else if(DD[j*numDisp+i]==1)
				BB[j*numDisp+i]=1;
		}

    plhs[0] = mxCreateNumericMatrix(width*height,numDisp,mxDOUBLE_CLASS,mxREAL);
    S = (double*)mxGetData(plhs[0]);//gives a void*, cast it to int*
    Path1(dispVol,I,AA,BB,p1,sigma ,width, height,sz,numDisp,S1);
	Path2(dispVol,I,AA,BB,p1,sigma ,width, height,sz,numDisp,S2);
	Path3(dispVol,I,AA,BB,p1,sigma ,width, height,sz,numDisp,S3);
	Path4(dispVol,I,AA,BB,p1,sigma ,width, height,sz,numDisp,S4);
	for(i=0;i<width*height*numDisp;i++)
	{
		S[i] = S1[i]+S2[i]+S3[i]+S4[i];
	}

}
