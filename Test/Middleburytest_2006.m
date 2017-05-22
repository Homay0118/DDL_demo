clc
clear
close all
% UseParallelToolbox = true; % Set true if you want to take advantage of the Matlab parallel computing toolbox
% ParallelWorkers =4; % How many workers should be used by the parallel computing toolbox (should be equal or less the number of available CPU cores)
% 
% % Set up parallel computing toolbox
% if (UseParallelToolbox)
%     if isempty(gcp('nocreate'))
%           parpool(ParallelWorkers);
%     end
% end
% load images
imgNum=0;
Datasets={2001,20011,2003,2005,2006};
mkdir('MiddEval2results');
for Y=5:5   %1:5 if all the datasets are downloaded.
    setNum=Datasets{Y};
    if setNum==2001
        datasets={'barn1','barn2','bull','poster','sawtooth','Venus'};
        maxNum=6;
    elseif setNum==20011
        datasets={'Tsukuba'};
        maxNum=1;
    elseif setNum==2003
        datasets={'Cones','Teddy'};
        maxNum=2;
    elseif setNum==2005
        datasets={'Art','Books','Dolls','Laundry','Reindeer','Moebius'};
        maxNum=6;
    elseif setNum==2006
    datasets={'Bowling1','Aloe','Baby1','Baby2','Baby3','Bowling2','Cloth1','Cloth2','Cloth3','Cloth4', 'Flowerpots','Lampshade1',...
              'Lampshade2','Rocks1','Rocks2'};
%     maxNum=15; %If all the datasets are downloaded.
     maxNum =1;
    end
    for testimage =1: maxNum
        data=datasets{testimage};
        imgNum=imgNum+1;
        ErrorRate(imgNum)=mainfunction_2006(setNum,data);
    end
end
fprintf('%s = %.2f%s\n','Overall',mean(ErrorRate)*100,'%');
 