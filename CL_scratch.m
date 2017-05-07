%% BETA DISTRIBUTION

xx = linspace(0,1,100);
plot(xx,betapdf(xx,1.3,4))

%% gamma distribution

Jbar = 35;
tau = 95;
samples = gamrnd(Jbar/tau,tau,[1,1e6]);
% figure
subplot(1,2,2)
hist(samples,100)
axis([0 200, 0 1e5]);

%% check percent correct as a function of session
subjid = '6';
exptype = 'Delay';

filepath = 'experiment_data/output_mat/round4/';
filename = ['Exp_ChangeLocalization_' exptype '_subj' subjid '.mat'];

load([filepath filename],'progPC')
plot(progPC,'ko'); hold on
plot([0 60],[25 25],'Color',0.7*ones(1,3)); % chance
plot([10 10;20 20;30 30;40 40;50 50],[0 100],'Color',0.7*ones(1,3))
ylim([0 100])
defaultplot;
hold off
xlabel('block number')
ylabel('PC')


%% check percent correct as a function of delta and confidence
clear all

subjid = 'JP';
exptype = 'Contrast';

filepath = 'experiment_data/output_mat/round4/';
filename = ['Exp_ChangeLocalization_' exptype '_subj' subjid '.mat'];
load([filepath filename],'designMat','stimuliMat','names')

% plotting as a function of delta
deltaVec = unique(abs(designMat(:,1)));
nDeltas = length(deltaVec);
PCVec = nan(1,nDeltas);
for idelta = 1:nDeltas
    delta = deltaVec(idelta);
    idx = (abs(designMat(:,1)) == delta);
    
    PCVec(idelta) = nanmean(designMat(idx,10));
end
figure;
plot(deltaVec,PCVec,'ok')
xlabel('delta')
ylabel('PC')
ylim([0 1])
defaultplot

% plotting as a function of confidence
nConf = 4;
PCVec = nan(1,nConf);
for iconf = 1:nConf
    idx = (designMat(:,12) == iconf);
    
    PCVec(iconf) = mean(designMat(idx,10));
end
figure;
plot(PCVec,'ok')
xlabel('confidence')
ylabel('PC')
ylim([0 1])
defaultplot


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




%% ============================================
%       LOOKING AT PARAMETER FITS
% ============================================

%% COMPILE MODEL FITS FOR ALL SUBJECTS
clear all

modelVec = [1 2];
subjids = {'1','2','3','4','5','6'};
exptypeVec = {'Contrast','Delay'};
filepath = 'analysis/2_fitdata/fits/';

nSubjs = length(subjids);
nModels = length(modelVec);
nExptype = length(exptypeVec);
nParams = 4;

for imodel = 1:nModels
    model = modelVec(imodel)
    modelstr = ['model' num2str(model)];
    
    for iexptype = 1:nExptype
        exptype = exptypeVec{iexptype}
        
        for isubj = 1:nSubjs
            isubj
            subjid = subjids{isubj};
            filename = sprintf('%sChangeLocalization_%s_model%d_subj%s.txt',filepath,exptype,model,subjid);
            alldata = dlmread(filename);
            
            datasorted = sortrows(alldata,nParams+1);
            fits.(modelstr).(exptype).bfpMat(isubj,:)= datasorted(end,1:nParams);
            fits.(modelstr).(exptype).LLVec(isubj)= datasorted(end,nParams+1);
        end
    end
    
end

save([filepath 'modelfits.mat'],'fits')

%% MODEL COMPARISON

nParams = 4;
Contrast_AICMat = -2.*[fits.model1.Contrast.LLVec; fits.model2.Contrast.LLVec]+2*nParams;
Delay_AICMat = -2.* [fits.model1.Delay.LLVec; fits.model2.Delay.LLVec]+2*nParams;
modcomp = Contrast_AICMat + Delay_AICMat;
nSubj = size(modcomp,2);

comparetype = 2;
switch comparetype
    case 1 % just contrast
        comparevalue = Contrast_AICMat;
    case 2 % just delay
        comparevalue = Delay_AICMat;
    case 3 % both
        comparevalue = modcomp;
end

% positive number means 
moddiff = bsxfun(@minus,comparevalue,comparevalue(1,:))
mean_moddiff = mean(moddiff,2)
sem_moddiff = std(moddiff,[],2)./sqrt(nSubj)
figure;
bar(moddiff')

%% =======================
%  PLOTS
% ===================

clear all

filepath = 'analysis/2_fitdata/fits/';
load([filepath 'modelfits.mat'])
% subjids = {'1','2','3','4','5','6'};
subjids = {'6'};
exptype = 'Contrast';
model = 1;
modelstr = sprintf('model%d',model);

nSubj = length(subjids);
deltaVec = 2.5:5:87.5;

for isubj = 1:nSubj
subjid = subjids{isubj};
bfp = fits.(modelstr).(exptype).bfpMat(str2double(subjid),:);
try
if bfp(4) == 0
    bfp(4) = 1e-6;
end
end

X = []; % data in luigiform
Nsamples = 1000; % default if empty

[~,prmat,X] = AhyBCL_datalikeall(bfp,X,model,Nsamples);

% plot psychometric function
figure;
hold on;
colors = aspencolors(4,'pastel'); %['b'; 'y'; 'g'; 'r'];
condlegendd = {'0/2', '2/1','2/2','4/1'};
for icond = 1:4
    subplot(2,2,icond)
    plot_psychometricfunction({subjid},exptype,icond); % plot real data
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

%% average plot

clear all

filepath = 'analysis/2_fitdata/fits/';
load([filepath 'modelfits.mat'])
subjids = {'1','2','3','4','5','6'};
% subjids = {'ALM','DR','EN','MR'};
exptype = 'Delay';
model = 1;

modelstr = sprintf('model%d',model);
nSubj = length(subjids);
deltaVec = 2.5:5:87.5;

PRMat = cell(1,4);
for isubj = 1:nSubj
    subjid = subjids{isubj};
    bfp = fits.(modelstr).(exptype).bfpMat(isubj,:);
    try
    if bfp(4) == 0
        bfp(4) = 1e-6;
    end
    end
    
    X = []; % data in luigiform
    Nsamples = 1000; % default if empty
    
    [~,prmat,X] = AhyBCL_datalikeall(bfp,X,model,Nsamples);
    
    PRMat = cellfun(@(x,y) [x;y(:,1)'],PRMat,prmat,'UniformOutput',false);
    
end
mean_mod = cellfun(@mean,PRMat,'UniformOutput',false);
sem_mod = cellfun(@(x) std(x)/sqrt(size(x,1)),PRMat,'UniformOutput',false);
mean_mod= reshape(cell2mat(mean_mod),[18 4])';
sem_mod = reshape(cell2mat(sem_mod),[18 4])';

% plot psychometric function
figure;
hold on;
colors = aspencolors(4,'pastel');%['b'; 'y'; 'g'; 'r'];
condlegendd = {'0/2', '2/1','2/2','4/1'};
for icond = 1:4
    
    subplot(2,2,icond)
    plot_psychometricfunction(subjids,exptype,icond); % plot real data
    
    plot_summaryfit(deltaVec,[],[],mean_mod(icond,:),sem_mod(icond,:),colors(icond,:),colors(icond,:));
    title(condlegendd{icond})
    axis([0 90 0 1])
    ax = gca;
    ax.XTick = [0 30 60 90];
    ax.YTick = 0:.5:1;
    %     fill([deltaVec fliplr(deltaVec)],[pc-pc_std fliplr(pc + pc_std)],colors(icond,:))
end
