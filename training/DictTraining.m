clc
clear
tic
%% Data loading
load all15/MatchPairs_2014_5_15.mat;
load all15/MatchPairs_ng_2014_5_15.mat;
K=256;
b=5;
param.K=K;
param.lambda=1.2/b;
param.numThreads=-1;
param.batchsize=512;
param.verbose=false;
param.iter=1000;
param.L=b^2;
%% Data preprocessing
L=double(MatchPairs.Left)./255;
R=double(MatchPairs.Right)./255;
R_ng=double(MatchPairs_ng.Right)./255;
L=L-repmat(mean(L),[size(L,1) 1]);
L=L ./ repmat(sqrt(sum(L.^2)),[size(L,1) 1]);
R=R-repmat(mean(R),[size(R,1) 1]);
R=R ./ repmat(sqrt(sum(R.^2)),[size(R,1) 1]);
R_ng=R_ng-repmat(mean(R_ng),[size(R_ng,1) 1]);
R_ng=R_ng ./ repmat(sqrt(sum(R_ng.^2)),[size(R_ng,1) 1]);
%% Remove the homogeneous columns
Ind=sum(isnan(L));
L(:,Ind~=0)=[];
R(:,Ind~=0)=[];
R_ng(:,Ind~=0)=[];
Ind=sum(isnan(R));
L(:,Ind~=0)=[];
R(:,Ind~=0)=[];
R_ng(:,Ind~=0)=[];
Ind=sum(isnan(R_ng));
L(:,Ind~=0)=[];
R(:,Ind~=0)=[];
R_ng(:,Ind~=0)=[];
%% Initial dictionary learning
X=[L,R,R_ng];
D = mexTrainDL(X, param);
save('all15/D0_Middlebury_2014_5_new.mat','D');
W=ones(K,size(L,2));
V=ones(K,1);
%% Iterating
for iter=1:200
 % Update Alpha
     Alphal = mexLassoWeighted(L,D,W,param);
     Alphar = mexLassoWeighted(R,D,W,param);
     Alphar_ng = mexLassoWeighted(R_ng,D,W,param);
 % Update D
     Alpha=[Alphal Alphar Alphar_ng];
     for i=1:K
         ai        =    Alpha(i,:);
         Y         =    X-D*Alpha+D(:,i)*ai;
         di        =    Y*ai';
         di        =    di./(norm(di,2) + eps);
         D(:,i)    =    di;
     end
 % Update W
    A=(Alphal-Alphar_ng)*(Alphal-Alphar_ng)';
    B=(Alphal-Alphar)*(Alphal-Alphar)';
    idx1=find(sum(A)==0);
    idx2=find(sum(B)==0);
    idx=unique([idx1,idx2]);
    A(idx,:)=[];
    A(:,idx)=[];
    B(idx,:)=[];
    B(:,idx)=[];
    [v,Eig]=eigs(A,B,1);
    idx1=find(sum(A)~=0);
    idx2=find(sum(B)~=0);
    idx_1=intersect(idx1,idx2);
    v=v./norm(v);
    v(abs(v)<10^-3)=0;
    V(idx)=0;
    V(idx_1)=v;
    W=1./(log(1+abs(V)/0.01)+0.0001);
    W=W./norm(W);
    W=repmat(W,1,size(L,2));
%  Compute the energy value
    P1 = X - D * Alpha;
    P1 = P1(:)'*P1(:) / 2;
    P2 = norm(Alpha, 1);
    f1=sum(abs(Alphal-Alphar));
    f2=sum(abs(Alphal-Alphar_ng));
    P3=full(sum(f2)/sum(f1));
    f(iter) =P1 + param.lambda*P2 - P3;
    disp('||X-D*a||_2    ||a||_2   ||V(a_r-a_l)||_2:   Total:');
    fprintf('%6.4f,     %6.4f,     %6.4f,      %6.4f\n',P1,param.lambda*P2,P3,f(iter));
%  Termination condition
    if iter==1
        continue;
    elseif abs(f(iter)-f(iter-1))/f(iter-1)<0.001
         save('all15/D_Middlebury_2014_5_new.mat','D','V','param','f');
         break;
    end
end

