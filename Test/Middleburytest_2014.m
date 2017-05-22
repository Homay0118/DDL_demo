clear variables;
warning('off','all')
close all;
% UseParallelToolbox = true; % Set true if you want to take advantage of the Matlab parallel computing toolbox
% ParallelWorkers = 4; % How many workers should be used by the parallel computing toolbox (should be equal or less the number of available CPU cores)
% 
% % Set up parallel computing toolbox
% if (UseParallelToolbox)
%  if isempty(gcp('nocreate'))
%  parpool(ParallelWorkers);
%  end
% end

% Are you going to use the training or test set?
%   imgset = 'training';
imgset = 'test';

% Specify which resolution you are using for the stereo image set (F, H, or Q?)
  imgsize = 'Q';

% What are you calling your method?
  methodname = 'DDL';

mkdir('MiddEval3results');
mkdir(['MiddEval3results/',imgset,imgsize]);
if strcmp(imgset,'training')
    image_names{1} = 'Adirondack';
    image_names{2} = 'ArtL';
    image_names{3} = 'Jadeplant';
    image_names{4} = 'Motorcycle';
    image_names{5} = 'MotorcycleE';
    image_names{6} = 'Piano';
    image_names{7} = 'PianoL';
    image_names{8} = 'Pipes';
    image_names{9} = 'Playroom';
    image_names{10} = 'Playtable';
    image_names{11} = 'PlaytableP';
    image_names{12} = 'Recycle';
    image_names{13} = 'Shelves';
    image_names{14} = 'Teddy';
    image_names{15} = 'Vintage';
    ndisp = [290, 256, 640, 280, 280, 260, 260, 300, 330, 290, 290, 260, 240, 256, 760];
else
    image_names{1} = 'Australia';
    image_names{2} = 'AustraliaP';
    image_names{3} = 'Bicycle2';
    image_names{4} = 'Classroom2';
    image_names{5} = 'Classroom2E';
    image_names{6} = 'Computer';
    image_names{7} = 'Crusade';
    image_names{8} = 'CrusadeP';
    image_names{9} = 'Djembe';
    image_names{10} = 'DjembeL';
    image_names{11} = 'Hoops';
    image_names{12} = 'Livingroom';
    image_names{13} = 'Newkuba';
    image_names{14} = 'Plants';
    image_names{15} = 'Staircase';
    ndisp = [290, 290, 250, 610, 610, 256, 800, 800, 320, 320, 410, 320, 570, 320, 450];
end
%------------------------------ load dictionary ---------------------------%
if strcmp(imgset,'training')
%   load ..\training\last7\D_Middlebury_2014_5.mat
%   test_start = 1;
%   test_end = 8;
  load ..\training\first8\D_Middlebury_2014_5.mat
  test_start = 9;
  test_end = 15;
else
  load ..\training\all15\D_Middlebury_2014_5.mat
  test_start = 1;
  test_end = 15;
end
D(:,sum(D)==0)=[];
%------------------- Sparse Coding Parameters  ----------------------------%
b=5;
K=size(D,2);
param.K=K;
param.lambda=1.2/b;
param.numThreads=-1;
param.batchsize=512;
param.verbose=false;
param.iter=1000;
%--------------------- Cost Aggregation Parameters  -----------------------%
P1=1.2;
sigma=0.4;
ErrorRate = zeros(1,15);
%--------------------------------- Executing ------------------------------%
% for im_num = test_start: test_end   % if all the datasets are downloaded.
for im_num =2:2
    %----------- Adjust the range of disparities to the chosen resolution -------------------------------%
    if imgsize == 'Q'
        DisparityRange = [1,round(ndisp(im_num)/4)];
    elseif imgsize == 'H'
        DisparityRange = [1,round(ndisp(im_num)/2)];
    else
        DisparityRange = [1,round(ndisp(im_num))];
    end
    
    tic
    %------------------------------- input color images  ---------------------------------------------%
    I1= imread(['../MiddEval3/',imgset,imgsize,'/',image_names{im_num},'/im0.png']);
    I2= imread(['../MiddEval3/',imgset,imgsize,'/',image_names{im_num},'/im1.png']);
    
    %-------------------------------------- data preparation -------------------------%
    [Data0_0,Data1_0,Data0_1,Data1_1]=DataPreparing(I1,I2,b);
    [M,N,~]=size(I1);
    
    %---------------------------------- sparse coding -------------------------------%
    Alpha1 = mexLasso(Data0_0, D, param);
    Alpha2 = mexLasso(Data1_0, D, param);
    Alpha1_1 = mexLasso(Data0_1, D, param);
    Alpha2_1 = mexLasso(Data1_1, D, param);
    
   %----------------------------- stereo matching ---------------------------------%
    dispVol00=CostComputingmex(full(Alpha1),full(Alpha2),DisparityRange(2),M,N);
    dispVol100=CostComputingmex(full(Alpha1_1),full(Alpha2_1),DisparityRange(2),M,N);
    dispVol00=reshape(dispVol00,[M,N,DisparityRange(2)]);
    dispVol100=flip(reshape(dispVol100,[M,N,DisparityRange(2)]),2);
    
    %-------------------------- smooth constraint--------------------------%
    I1=double(I1)/255;
    I2=double(I2)/255;
    dispVol0=SGMmex(dispVol00,I1,P1,sigma);
    dispVol10=SGMmex(dispVol100,I2,P1,sigma);
    dispVol0=reshape(dispVol0,[M,N,DisparityRange(2)]);
    dispVol10=reshape(dispVol10,[M,N,DisparityRange(2)]);
    %----------------------------Disparity Computing-----------------------%
    [~,DisparityMap1]=min(dispVol0,[],3);
    [~,DisparityMap2]=min(dispVol10,[],3);
    
    %----------------- Left-right consistency check -----------------------%
    [M,N,c]=size(I1);
    Y = repmat((1:M)', [1 N]);
    X = repmat(1:N, [M 1]) - DisparityMap1;
    X(X<1) = 1;
    indices = sub2ind([M,N],Y,X);
    final_labels = DisparityMap1;
    final_labels(abs(DisparityMap1 - DisparityMap2(indices))>=1) = -1;
    
    %------------ Fill and filter (post-process) pixels that fail the consistency check ------------%
    inputLabels = int32(final_labels);
    final_labels = PostProcmex(inputLabels,I1,DisparityRange(2));
    final_labels = double(final_labels);
    time_taken=toc;
    
    %-------------------------------If possible, compute the error rate------------------------------%
    if strcmp(imgset,'training')
        GT = readpfm(['../MiddEval3GT/training',imgsize,'/',image_names{im_num},'/disp0GT.pfm']);
        mask = imread(['../MiddEval3GT/training',imgsize,'/',image_names{im_num},'/mask0nocc.png']);
        mask = mask == 255;
        Error = abs(final_labels- GT) > 1;
        Error(~mask) = 0;
        ErrorRate(im_num) = sum(Error(:))/sum(mask(:));
%         fprintf('%s = %f\n', image_names{im_num}, ErrorRate(im_num));
    end
    
    %--------------------------.png format-------------------------------------------%
    % Output disparity maps and timing results for Middlebury evaluation
    mkdir(['MiddEval3results/',imgset,imgsize,'/',image_names{im_num}]);
    pfmwrite(single(final_labels), ['MiddEval3results/',imgset,imgsize,'/',image_names{im_num},'/disp0',methodname,'.pfm']);
    fid = fopen(['MiddEval3results/',imgset,imgsize,'/',image_names{im_num},'/time',methodname,'.txt'],'w');
    fprintf(fid,'%f',time_taken);
    fclose(fid);
end
if strcmp(imgset,'training')
    ErrorRateMean = (sum(ErrorRate([1,2,3,4,5,6,8,11,12,14])) + 0.5*sum(ErrorRate([7,9,10,13,15])))/12.5;
    fprintf('Overall = %f%s\n', ErrorRateMean*100,'%');
end

                 
          
                      
            

