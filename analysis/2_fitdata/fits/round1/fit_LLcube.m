% function fit_LLcube(subjid, model)

subjid = 'PM';
model = 1;
filetail = [];
nSamples = 1e5;

% loading real data
X = aldkfjl; % data struct/cell

% loading LLcube
load(sprintf('LLcube_model%d_subjid%s%s.mat',...
    model,subjid,filetail));
v2struct(data);

% getting all parameter combinations
Jbar1Vec = exp(linspace(gridMat(1,1),gridMat(1,2),gridMat(1,3)));
Jbar2Vec = exp(linspace(gridMat(2,1),gridMat(2,2),gridMat(2,3)));
tauVec = exp(linspace(gridMat(3,1),gridMat(3,2),gridMat(3,3)));
[x,y,z] = meshgrid(Jbar1Vec,Jbar2Vec,tauVec);
allParamCombo = [x(:) y(:) z(:)];

% find best fitting parameter
idx = find(LLMat(:) == max(LLMat(:)));
bestFitParam = allParamCombo(idx,:);


[~,prmat] = AhyBCL_datalikeall(bestFitParam,X,model,nSamples);