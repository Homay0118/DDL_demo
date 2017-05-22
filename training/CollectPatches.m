% This script is used to collect training samples
clc
clear
close all
b=5;
filenames=dir('../MiddEval3GT/trainingQ/');
filenames(1:2)=[];
Left=[];
Right=[];
Left_ng=[];
Right_ng=[];
for im=1:15   %or 1:8, change the save path to first8\ or 9:15,  change the save path to last7\
%% left view
names=filenames(im).name;
I=imread(strcat('../MiddEval3/trainingQ/',names,'/im0.png'));
I0=rgb2gray(I);
%% right view
names=filenames(im).name;
I=imread(strcat('../MiddEval3/trainingQ/',names,'/im1.png'));
I1=rgb2gray(I);
%% image to columns
Data00=im2col(I0,[b,b],'sliding');
Data10=im2col(I1,[b,b],'sliding');
Data0=Data00;
Data1=Data10;
%% disp
names=filenames(im).name;
Disp=round(double(readpfm(strcat('../MiddEval3GT/trainingQ/',names,'/disp0GT.pfm'))));
mask0nocc=double(imread(strcat('../MiddEval3GT/trainingQ/',names,'/mask0nocc.png')))==255;
Mask=im2col(~mask0nocc,[b,b],'sliding');
ind_remove=find(sum(Mask)~=0); % remove occluded regions
%% matching pairs
[M,N]=size(Disp);
[X,Y]=meshgrid(1:N,1:M);
Disp=Disp((b+1)/2:M-(b-1)/2,(b+1)/2:N-(b-1)/2);
X=X((b+1)/2:M-(b-1)/2,(b+1)/2:N-(b-1)/2);
Y=Y((b+1)/2:M-(b-1)/2,(b+1)/2:N-(b-1)/2);
M1=M-b+1;
N1=N-b+1;
X1=min(max((b+1)/2,X-Disp+randi([-1,1],M1,N1)),N-(b-1)/2);
X1=X1-(b-1)/2;
Y1=Y-(b-1)/2;
IND=sub2ind([M1,N1],Y1,X1);
Data1=Data1(:,IND);
Data0(:,ind_remove)=[];
Data1(:,ind_remove)=[];
%% selecting patches with large variance
Mean_left=repmat(mean(double(Data0)),[size(Data0,1),1]);
Var_left=sum((double(Data0)-Mean_left).^2);
[Var_left,ind1]=sort(Var_left,'descend');
K0=size(Data0,2)/2;
Data0=Data0(:,ind1(1:K0));
Data1=Data1(:,ind1(1:K0));
%% randomly select 20000 patch pairs
sampling=randperm(length(Data0),min(length(Data0),20000));
Left=[Left,Data0(:,sampling)];
Right=[Right,Data1(:,sampling)];
%% non-matching pairs
X2=min(max((b+1)/2,X-Disp+randi([2 8],M1,N1).*(-1).^randi([0,1],M1,N1)),N-(b-1)/2);
X2=X2-(b-1)/2;
Y2=Y-(b-1)/2;
IND=sub2ind([M1,N1],Y2,X2);
Data11=Data10(:,IND);
Data11(:,ind_remove)=[];
Data11=Data11(:,ind1(1:K0));
Left_ng=[Left_ng,Data0(:,sampling)];
Right_ng=[Right_ng,Data11(:,sampling)];
end

%% save patches
MatchPairs.Left=Left;
MatchPairs.Right=Right;
save(['all15/MatchPairs_2014_',num2str(b),'_15.mat'],'MatchPairs');
MatchPairs_ng.Left=Left_ng;
MatchPairs_ng.Right=Right_ng;
save(['all15/MatchPairs_ng_2014_',num2str(b),'_15.mat'],'MatchPairs_ng');
disp('Done!')