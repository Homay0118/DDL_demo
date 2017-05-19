function [Data0_0,Data1_0,Data0_1,Data1_1]=DataPreparing(I1,I2,b)
   %% converte image to patches
    I1_g = double(rgb2gray(I1))./255;
    I2_g = double(rgb2gray(I2))./255;
    Data0=img2patch(I1_g,b);
    Data1=img2patch(I2_g,b);
   %% normalization
    Data0=Data0-repmat(mean(Data0),[size(Data0,1) 1]);
    Data0_0=Data0 ./ repmat(sqrt(sum(Data0.^2)),[size(Data0,1) 1]);
    Data1=Data1-repmat(mean(Data1),[size(Data1,1) 1]);
    Data1_0=Data1 ./ repmat(sqrt(sum(Data1.^2)),[size(Data1,1) 1]);
    ind=isnan(sum(Data0_0));
    ind=find(ind~=0);
    Data0_0(:,ind)=Data0(:,ind);
    ind=isnan(sum(Data1_0));
    ind=find(ind~=0);
    Data1_0(:,ind)=Data1(:,ind);
   %% another view
    I1_g_1 = flip(double(rgb2gray(I2))/255,2);
    I2_g_1 = flip(double(rgb2gray(I1))/255,2);
    Data0=img2patch(I1_g_1,b);
    Data1=img2patch(I2_g_1,b);
    Data0=Data0-repmat(mean(Data0),[size(Data0,1) 1]);
    Data0_1=Data0 ./ repmat(sqrt(sum(Data0.^2)),[size(Data0,1) 1]);
    Data1=Data1-repmat(mean(Data1),[size(Data1,1) 1]);
    Data1_1=Data1 ./ repmat(sqrt(sum(Data1.^2)),[size(Data1,1) 1]);
    ind=isnan(sum(Data0_1));
    ind=find(ind~=0);
    Data0_1(:,ind)=Data0(:,ind);
    ind=isnan(sum(Data1_1));
    ind=find(ind~=0);
    Data1_1(:,ind)=Data1(:,ind);
end