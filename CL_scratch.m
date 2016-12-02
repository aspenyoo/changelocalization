%% % % % % % % % % % % % % % % % % % % % % % % % % 
%    BETA DISTRIBUTION
% % % % % % % % % % % % % % % % % % % % % % % % % 
xx = linspace(0,1,100);
plot(xx,betapdf(xx,1.3,4))


%% logistic regression: pc ~ condition + delta
clear all
clc

subjid = 'AHY';
load(sprintf('data_ChangeLocalization_subj%s.mat',subjid));

deltaVec = 2.5:5:87.5;              % possible change values
respVec = [3 0 4 0 0 0 2 0 1];      % go from response key to location index

designMat = designMat(designMat(:,2) ~= 2,:);

% set X and y matrices
y = designMat(:,7); % correct/not correct
X = [designMat(:,2) abs(designMat(:,1))]; % delay, delta
X(X(:,1) == 0,2) = 1; % delay for 0 condition is 1
X(X(:,2) == 4,2) = 3; % delay for 4 condition is 3

% logistic regression
mdl = fitglm(X,y,'Distribution','binomial')
ci = coefCI(mdl)

%% linear regression: confidence ~ condition + delta

clear all

subjid = 'AHY';
load(sprintf('data_ChangeLocalization_subj%s.mat',subjid));

deltaVec = 2.5:5:87.5;              % possible change values
respVec = [3 0 4 0 0 0 2 0 1];      % go from response key to location index

designMat = designMat(designMat(:,2) ~= 2,:);

% set X and y matrices
y = designMat(:,9); % confidence rating
X = [designMat(:,2) abs(designMat(:,1))]; % delay, delta
X(X(:,1) == 0,2) = 1; % delay for 0 condition is 1
X(X(:,2) == 4,2) = 3; % delay for 4 condition is 3

% linear regression
X = [ones(size(designMat,1),1) X];
[b,bint,~,~,stats] = regress(y,X) % stats: R^2, F, p, and estimate of error variance

% % ordinal multinomial regression
% [B,dev,stats] = mnrfit(X,y,'model','ordinal')
% LL = stats.beta - 1.96*stats.se
% UL = stats.beta + 1.96*stats.se
%
% % plot regression fits
% xx = linspace(5,85,50);
% XX = [ones(50,1) 3*ones(50,1) xx'];
% BMat = [B(1:3)'; repmat(B(4:end),1,3)];
% logoddsMat = XX*BMat; % log odds for each of the 3 equations
% logoddssign = sign(logoddsMat);
% difflogoddssign = diff(logoddssign,1,2); % difference between the signs of log odds
% [row,col] = find(difflogoddssign == 2);
% predictedConf = nan(50,1);
% predictedConf(row) = col+1;
% predictedConf(logoddssign(:,1) == 1 ) = 1;
% predictedConf(logoddssign(:,3) == -1) = 4;
% hold on; plot(xx,predictedConf);


%% create model fit datat psychometric fn
clear all

load('modelfits.mat')
model = 2;
subjids = {'AHY','LA'};
nSubj = length(subjids);
deltaVec = 2.5:5:87.5;

for isubj = 1:nSubj;
    subjid = subjids{isubj};
    bfp = squeeze(xbest(isubj,model,:));
    
    X = []; % data in luigiform
    Nsamples = 1e5; % default if empty
    
    [~,prmat,X] = AhyBCL_datalikeall(bfp,X,model,Nsamples);
    data = X;
    save(sprintf('experiment_data/grant_data/luigidata_ChangeLocalization_subj%s_model%d.mat',subjid,model),'data','prmat');
end

for isubj = 1:nSubj;
    subjid = subjids{isubj};
    load(sprintf('experiment_data/grant_data/luigidata_ChangeLocalization_subj%s_model%d.mat',subjid,model));
    
    % plot psychometric function
    figure;
    hold on;
    colors = aspencolors(4,'pastel');%['b'; 'y'; 'g'; 'r'];
    condlegendd = {'0/2', '2/1','2/2','4/1'};
    for icond = 1:4;
        subplot(2,2,icond)
        plot_psychometricfunction({subjid},icond); % plot real data
        pc = prmat{icond}(:,1)';
        pc_std = pc.*(1-pc)./24;
        plot_summaryfit(deltaVec,[],[],pc,pc_std,colors(icond,:),colors(icond,:));
        title(condlegendd{icond})
        axis([0 90 0 1])
        ax = gca;
        ax.XTick = [0 30 60 90];
        ax.YTick = 0:.5:1;
        %     fill([deltaVec fliplr(deltaVec)],[pc-pc_std fliplr(pc + pc_std)],colors(icond,:))
    end
    
end

%% LOOK AT SOME STIMULI ON PSYCHTOOLBOX

% line stimuli
numDegrees = 180;
% prefs.lineArea = prefs.stimArea; % so line is same area as ellipse
lineWidth = 3; % pixels
lineLength = 100;
lineCoordinates = nan(2,numDegrees);
for j = 1:numDegrees;
    [lineCoordinates(1,j), lineCoordinates(2,j)] = lineCoord(lineLength,j);
end

% psychtoolbox
screenNumber = max(Screen('Screens'));       % use external screen if exists
[w, h] = Screen('WindowSize', screenNumber);  % screen resolution of smaller display
screenCenter = [w h]/2;
windowPtr = Screen('OpenWindow',screenNumber,128*ones(1,3),[],32,2);
Screen(windowPtr,'BlendFunction',GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

% make gabor patch with full contrast
stimLength = 100;
StimSizes = [stimLength stimLength];
reliability = 1; % contrast of 1
im = makeGabor(stimLength, reliability);
StimPatch= Screen('MakeTexture',windowPtr,im);


pos = screenCenter + [100 100];
lineVec = [0 30];
nLines = length(lineVec);
for iline = 1:nLines;
    deg = mod(90-lineVec(iline),180);
    if deg == 0; deg = 180; end
    
    % gabor
    srcrect = [0 0 StimSizes];
    destrect = CenterRectOnPoint(srcrect,pos(1),pos(2));
    Screen('DrawTexture', windowPtr, StimPatch, srcrect,destrect, 180-90+deg);
    
    % line
    lineStim = lineCoordinates(:,deg);
    xy = [-lineStim lineStim];
    Screen('DrawLines',windowPtr, xy, lineWidth,256,pos,1);
    
    Screen('flip',windowPtr); % tic;
    pause;
end

sca;
ShowCursor;


%% FIT MOC psychometric function
% 08.16.2016

clear all

mixedcontrast = 1; 
subjid = 'AMM';
psychofun = @(x) lambda/2 + (1-lambda).*0.5*(1+erf((x-mu)./(sqrt(2)*sigma)));

if (mixedcontrast)
    concatvars = concatcode('experiment_data/output_mat/',['ContrastFn_MOC_' subjid '_'],{'designMat','stimuliMat'});
    v2struct(concatvars);
%     load('ContrastFn_MOC_MIX_20160804T220952.mat')
    designMat = designMat(stimuliMat(:,5) >= 3,:);
    titlee = 'mixed contrast';
else
    concatvars = concatcode('experiment_data/output_mat/',['ContrastFn_MOC_' subjid '_'],{'designMat','stimuliMat'});
    v2struct(concatvars);
%     load('ContrastFn_MOC_AHY_20160804T215847.mat')
titlee = 'same contrast';
end

contrastVec = unique(designMat(:,1));
nContrasts = length(contrastVec);

PC = nan(1,nContrasts);
for icontrast = 1:nContrasts;
    contrast = contrastVec(icontrast);
    
    PC(icontrast) = nanmean(designMat(designMat(:,1) == contrast,3));
end

figure
plot((contrastVec),PC,'k*')
defaultplot
ylabel('percent correct')
xlabel('contrast')
ylim([0.5 1])
title(titlee)


%% COMPILE MODEL FITS FOR ALL SUBJECTS

modelVec = [1 2];
subjids = {'ALM','DR','EN','MR'};
exptypeVec = {'Contrast','Delay'};

nSubjs = length(subjids);
nModels = length(modelVec);
nExptype = length(exptypeVec);
nParams = 3;


for imodel = 1:nModels
    model = modelVec(imodel);
    modelstr = ['model' num2str(model)];
    
    for iexptype = 1:nExptype
        exptype = exptypeVec{iexptype};
        
        
        for isubj = 1:nSubjs
            subjid = subjids{isubj};
            filename = sprintf('ChangeLocalization_%s_model%d_subj%s.txt',exptype,model,subjid);
            alldata = dlmread(filename);
            
            datasorted = sortrows(alldata,nParams+1);
            fits.(modelstr).(exptype).bfpMat(isubj,:)= datasorted(1,1:nParams);
            fits.(modelstr).(exptype).LLVec(isubj)= datasorted(1,nParams+1);
        end
    end
    
end


%% MODEL COMPARISON

modcomp = [fits.model1.Contrast.LLVec; fits.model2.Contrast.LLVec];
modcomp = modcomp + [fits.model1.Delay.LLVec; fits.model2.Delay.LLVec]

% positive number means 
bsxfun(@minus,modcomp(1,:),modcomp)