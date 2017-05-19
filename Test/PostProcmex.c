

#include<mex.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>

typedef int int32_t;

int comp(const void *a,const void *b)
{
	return ((*(double*)a-*(double*)b>0)?1:-1);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    int *dims,index;
	double *I, *I_filter, *hist, *I_padding;
	int *disp,*disp_filled,*medianFiltered,left_d,right_d;
	int width,height,i,j,k,numDisp,sz,numelements,i_start,i_end,j_start,j_end;
	double gama_c = 0.3;
	double gama_d = 8;
	double r_c,g_c,b_c;
	double hist_sum =0, hist_cumsum =0;
	int  R = 8,r1,r2;
	int  flag_l, flag_r;
	double *temp_r;
	double *temp_g;
	double *temp_b;
	double temp_rgb,temp_dist;
    if (nrhs < 1) {
        mexErrMsgTxt("At least one argument is required.") ;
    } else if(nrhs > 6) {
        mexErrMsgTxt("Too many input arguments.");
    }


    numelements   = mxGetNumberOfElements(prhs[0]) ;
    dims  = mxGetDimensions(prhs[0]) ;
    disp = (int*)mxGetData(prhs[0]) ;//mxGetData returns a void pointer, so cast it

    width = dims[1]; height = dims[0];//Note: first dimension provided is height and second is width
    sz = width*height;
	
	I = (double*)mxGetData(prhs[1]) ;
	numDisp = mxGetScalar(prhs[2]);
	
	disp_filled = mxMalloc(sizeof(int)*height*width);
    medianFiltered = mxMalloc(sizeof(int)*height*width);
	I_padding = mxMalloc(sizeof(double)*(height + 2)*(width + 2) * 3);
	I_filter = mxMalloc(sizeof(double)*height*width* 3);
	hist = mxMalloc(sizeof(double)*numDisp);
	temp_r = mxMalloc(sizeof(double)*9);
	temp_g = mxMalloc(sizeof(double)*9);
	temp_b = mxMalloc(sizeof(double)*9);

/*******************Streak-based filling**************************/

	plhs[0] = mxCreateNumericMatrix(height,width, mxINT32_CLASS, mxREAL);
	medianFiltered = (int *)mxGetData(plhs[0]);//gives a void*, cast it to int*
	//plhs[0] = mxCreateNumericMatrix(height*width,3, mxDOUBLE_CLASS, mxREAL);
	//I_filter = (double *)mxGetData(plhs[0]);
	for(i=0;i<height;i++)
		for(j=0;j<width;j++)
		{   
		    disp_filled[j*height + i] = disp[j*height + i];
			if(disp[j*height+i]<0)
			{
				flag_r = 0;
				flag_l = 0;
					for(k=j+1;k<width;k++)
						{
						  if (disp[k*height+i]>0)
							{
							  right_d=disp[k*height+i];
							  flag_r=1;
							  break;
							}
						}
					for(k=j-1;k>=0;k--)
						{
						  if (disp[k*height+i]>0)
							{
								left_d=disp[k*height+i];
								flag_l=1;
								break;
							}
						}
				
					if (flag_l==1&&flag_r==1)
					{
						disp_filled[j*height+i]=min(left_d,right_d);
					}
					else if (flag_l==1&&flag_r==0)
					{
						disp_filled[j*height+i]=left_d;
					}
					else if (flag_l==0&&flag_r==1)
					{
							disp_filled[j*height+i]=right_d;
					}
					else
					{
						disp_filled[j*height+i]=0;
					}
			}
		}
/******************medfilt2**********************/
	for (i = 0; i < (height + 2)*(width + 2) * 3;i++) //padding
		I_padding[i] = 0;
	
	for (i = 1; i < height+1; i++)
		for (j = 1; j < width+1;j++) 
			{
			I_padding[j*(height + 2) + i] = I[(j - 1)*height + i - 1];
			I_padding[j*(height + 2) + i + (height + 2)*(width + 2)] = I[(j - 1)*height + i - 1 + sz];
			I_padding[j*(height + 2) + i + 2 * (height + 2)*(width + 2)] = I[(j - 1)*height + i - 1 + 2 * sz];
			}

		for(i=1;i<height+1;i++)
			for(j=1;j<width+1;j++)
			{
			temp_r[0] = I_padding[(j - 1)*(height+2) + i - 1];
			temp_r[1] = I_padding[j*(height + 2) + i - 1];
			temp_r[2] = I_padding[(j + 1)*(height + 2) + i - 1];
			temp_r[3] = I_padding[(j - 1)*(height + 2) + i];
			temp_r[4] = I_padding[j*(height + 2) + i];
			temp_r[5] = I_padding[(j + 1)*(height + 2) + i];
			temp_r[6] = I_padding[(j - 1)*(height + 2) + i + 1];
			temp_r[7] = I_padding[j*(height + 2) + i + 1];
			temp_r[8] = I_padding[(j + 1)*(height + 2) + i + 1];

			temp_g[0] = I_padding[(j - 1)*(height + 2) + i - 1 + (height + 2)*(width + 2)];
			temp_g[1] = I_padding[j*(height + 2) + i - 1 + (height + 2)*(width + 2)];
			temp_g[2] = I_padding[(j + 1)*(height + 2) + i - 1 + (height + 2)*(width + 2)];
			temp_g[3] = I_padding[(j - 1)*(height + 2) + i + (height + 2)*(width + 2)];
			temp_g[4] = I_padding[j*(height + 2) + i + (height + 2)*(width + 2)];
			temp_g[5] = I_padding[(j + 1)*(height + 2) + i + (height + 2)*(width + 2)];
			temp_g[6] = I_padding[(j - 1)*(height + 2) + i + 1 + (height + 2)*(width + 2)];
			temp_g[7] = I_padding[j*(height + 2) + i + 1 + (height + 2)*(width + 2)];
			temp_g[8] = I_padding[(j + 1)*(height + 2) + i + 1 + (height + 2)*(width + 2)];

			temp_b[0] = I_padding[(j - 1)*(height + 2) + i - 1 + 2 * (height + 2)*(width + 2)];
			temp_b[1] = I_padding[j*(height + 2) + i - 1 + 2 * (height + 2)*(width + 2)];
			temp_b[2] = I_padding[(j + 1)*(height + 2) + i - 1 + 2 * (height + 2)*(width + 2)];
			temp_b[3] = I_padding[(j - 1)*(height + 2) + i + 2 * (height + 2)*(width + 2)];
			temp_b[4] = I_padding[j*(height + 2) + i + 2 * (height + 2)*(width + 2)];
			temp_b[5] = I_padding[(j + 1)*(height + 2) + i + 2 * (height + 2)*(width + 2)];
			temp_b[6] = I_padding[(j - 1)*(height + 2) + i + 1 + 2 * (height + 2)*(width + 2)];
			temp_b[7] = I_padding[j*(height + 2) + i + 1 + 2* (height + 2)*(width + 2)];
			temp_b[8] = I_padding[(j + 1)*(height + 2) + i + 1 + 2 * (height + 2)*(width + 2)];

			qsort(temp_r, 9, sizeof(temp_r[0]), comp);
			qsort(temp_g, 9, sizeof(temp_g[0]), comp);
			qsort(temp_b, 9, sizeof(temp_b[0]), comp);

				I_filter[(j-1)*height + i-1] = temp_r[4];
				I_filter[(j-1)*height + i-1 + sz] = temp_g[4];
				I_filter[(j-1)*height + i-1 + 2 * sz] = temp_b[4];
			}
/***********weightMedian***********/
    for (k = 0; k < numDisp; k++)
			hist[k] = 0;

	for(i=0;i<height;i++)
		for (j = 0; j < width; j++)
		{
			medianFiltered[j*height + i] = disp[j*height + i];
			if (disp[j*height + i] < 0)
			{
				r_c = I_filter[j*height + i];
				g_c = I_filter[j*height + i + sz];
				b_c = I_filter[j*height + i + 2 * sz];
				i_start = max(0, i - R);
				i_end = min(height - 1, i + R);
				j_start = max(0, j - R);
				j_end = min(width - 1, j + R);
					for (r1 = i_start; r1 <= i_end; r1++)
						for (r2 = j_start; r2 <= j_end; r2++)
						{
						temp_rgb = (r_c - I_filter[r2*height + r1])*(r_c - I_filter[r2*height + r1]) + (g_c - I_filter[r2*height + r1 + sz])*(g_c - I_filter[r2*height + r1 + sz]) + (b_c - I_filter[r2*height + r1 + 2 * sz])*(b_c - I_filter[r2*height + r1 + 2 * sz]);
						temp_rgb = sqrt(temp_rgb) / (gama_c*gama_c);
						temp_dist = (r1 - i)*(r1 - i) + (r2 - j)*(r2 - j);
						temp_dist = sqrt(temp_dist) / (gama_d*gama_d);
						index = disp_filled[r2*height + r1];
						hist[index-1] = hist[index-1] + exp(-1*temp_rgb)* exp(-1*temp_dist);
						hist_sum = hist_sum + exp(-1 * temp_rgb)* exp(-1 * temp_dist);
						}
					for (k = 0; k<numDisp; k++)
					{
						hist_cumsum = hist[k] + hist_cumsum;
						if (hist_cumsum > hist_sum*0.5)
							{
								medianFiltered[j*height + i] = k+1;
								break;
							}
					}
			hist_sum = 0;
			hist_cumsum = 0;
			for (k = 0; k < numDisp; k++)
				    hist[k] = 0;
			}
		}
}
