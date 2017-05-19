function [GT,mask,numDisp]=nonoccmask(disp1,disp2,scale)
   GT=double(disp1)./scale;
   disp1=fix(double(disp1)./scale);
   disp2=fix(double(disp2)./scale);
   [M,N]=size(disp1);
   [X,Y]=meshgrid(1:N,1:M);
   X2 = X-disp1;
   X2(X2<1) = 1;
   indices = sub2ind([M,N],Y,X2);
   mask=abs(disp1-disp2(indices))==0;
   numDisp=max(disp1(:));
   numDisp=unique(numDisp);
   numDisp=ceil(numDisp/16)*16;
end