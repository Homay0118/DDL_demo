function ErrorRate=mainfunction_2006(year,testimage)
tic
%----- Load input images and Groundtruth, and set disparity ranges and the scale of disparity map to write------%
if year==2001
I1=imread(['.\img_stereo\',testimage,'\im2.ppm']);
I2=imread(['.\img_stereo\',testimage,'\im6.ppm']);    
scale=8;
GT1=imread(['.\img_stereo\',testimage,'\disp2.pgm']);
GT2=imread(['.\img_stereo\',testimage,'\disp6.pgm']);
[GT,mask,numDisp]=nonoccmask(GT1,GT2,scale);
elseif year==20011
I1=imread(['.\img_stereo\',testimage,'\imL.png']);
I2=imread(['.\img_stereo\',testimage,'\imR.png']);
scale=16;
GT=fix(double(imread(['.\img_stereo\',testimage,'\truedisp.pgm']))./scale);
mask=double(imread(['.\img_stereo\',testimage,'\nonocc.png']))==1;
numDisp=max(GT(:));
numDisp=unique(numDisp);
elseif year==2003
I1=imread(['.\img_stereo\',testimage,'\imL.png']);
I2=imread(['.\img_stereo\',testimage,'\imR.png']);
scale=4;
GT1=imread(['.\img_stereo\',testimage,'\disp2.png']);
GT2=imread(['.\img_stereo\',testimage,'\disp6.png']);
[GT,mask,numDisp]=nonoccmask(GT1,GT2,scale);
elseif year==2005||year==2006
I1=imread(['.\img_stereo\',testimage,'\view1.png']);
I2=imread(['.\img_stereo\',testimage,'\view5.png']);
scale=3;
GT1=imread(['.\img_stereo\',testimage,'\disp1.png']);
GT2=imread(['.\img_stereo\',testimage,'\disp5.png']);
[GT,mask,numDisp]=nonoccmask(GT1,GT2,scale);
else
I1=imread(['.\img_stereo\',testimage,'\view1.png']);
I2=imread(['.\img_stereo\',testimage,'\view5.png']);
scale=3;
GT1=imread(['.\img_stereo\',testimage,'\disp1.png']);
GT2=imread(['.\img_stereo\',testimage,'\disp5.png']);
[GT,mask,numDisp]=nonoccmask(GT1,GT2,scale);
end

%------------- load dictionary and remove the zero atoms-------------------%
load ..\training\all15\D_Middlebury_2014_5.mat;
D(:,sum(D)==0)=[];

%------------------------------ Parameter settings ------------------------%
b=5;
K=size(D,2);
param.K=K; 
param.lambda=1.2/b;
param.numThreads=-1; 
param.batchsize=512;
param.verbose=false;
param.iter=1000;

%-------------------------------------- data preparation -----------------%
[Data0_0,Data1_0,Data0_1,Data1_1]=DataPreparing(I1,I2,b);
[M,N,~]=size(I1); 

 %---------------------------------- sparse coding ------------------------%
 Alpha1 = mexLasso(Data0_0, D, param);
 Alpha2 = mexLasso(Data1_0, D, param);
 Alpha1_1 = mexLasso(Data0_1, D, param);
 Alpha2_1 = mexLasso(Data1_1, D, param);
     
 %----------------------------- stereo matching ---------------------------%
 dispVol00=CostComputingmex(full(Alpha1),full(Alpha2),numDisp,M,N);
 dispVol100=CostComputingmex(full(Alpha1_1),full(Alpha2_1),numDisp,M,N);
 dispVol00=reshape(dispVol00,[M,N,numDisp]);
 dispVol100=flip(reshape(dispVol100,[M,N,numDisp]),2);    
 %---------------------------- smooth constraint --------------------------%
 I1=double(I1)/255;
 I2=double(I2)/255;
 P1=1.2;
 sigma=0.4;
 dispVol0=SGMmex(dispVol00,I1,P1,sigma);
 dispVol10=SGMmex(dispVol100,I2,P1,sigma);
 dispVol0=reshape(dispVol0,[M,N,numDisp]);
 dispVol10=reshape(dispVol10,[M,N,numDisp]);
 [~,DisparityMap1]=min(dispVol0,[],3);
 [~,DisparityMap2]=min(dispVol10,[],3);     
%--------------------- Left-right consistency check -----------------------%
 [M,N,~]=size(I1);
 Y = repmat((1:M)', [1 N]);
 X = repmat(1:N, [M 1]) - DisparityMap1;
 X(X<1) = 1;
 indices = sub2ind([M,N],Y,X);
 final_labels = DisparityMap1;
 final_labels(abs(DisparityMap1 - DisparityMap2(indices))>=1) = -1;   
%------------ Fill and filter (post-process) pixels that fail the consistency check ------------%
 inputLabels = int32(final_labels);
 final_labels = PostProcmex(inputLabels,I1,numDisp);
 final_labels = double(final_labels);
 time_taken=toc;
%------------------------Compute the error rate ---------------------------%
 Error = abs(final_labels- GT) > 1;
 Error(~mask) = 0;
 ErrorRate = sum(Error(:))/sum(mask(:));
 fprintf('%s = %f\n',testimage, ErrorRate);
%--------------------------------.png format-------------------------------%
 %  Output disparity maps and timing results for Middlebury evaluation
 mkdir(['MiddEval2results\',testimage]);
 imwrite(uint8(final_labels*scale), ['MiddEval2results\',testimage,'\disp0.png']);
 fid = fopen(['MiddEval2results\',testimage,'\time.txt'],'w');   
 fprintf(fid,'%f',time_taken);
 fclose(fid);
end