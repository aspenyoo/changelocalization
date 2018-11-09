%% BETA DISTRIBUTION

xx = linspace(0,1,100);
plot(xx,betapdf(xx,1.3,4))

%% gamma distribution

Jbar = 0.64;
tau = 262;
samples = gamrnd(Jbar/tau,tau,[1,1e6]);
% xmax = .001;
% samples(samples > xmax) = [];
% figure
hist(samples,100)
axis([0 200, 0 1e5]);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%           DATA RELATED
% % % % % % % % % % % % % % % % % % % % % % % % % % %

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

%% % % % % % % % % % % % % % % % % % % % % % %
%           STATS
% % % % % % % % % % % % % % % % % % % % % % %

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

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%        MODEL RELATED
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

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
%       REAL DATA STUFF
% ============================================

%% COMPILE MODEL FITS FOR ALL SUBJECTS
clear all

modelVec = [1 2];
subjids = {'1','2','3','4','5','6','7','8','9'};
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
clear all

comparetype = 'Contrast';

filepath = 'analysis/2_fitdata/fits/';
load([filepath 'modelfits.mat'],'fits')

nParams = 4;
Contrast_AICMat = -2.*[fits.model1.Contrast.LLVec; fits.model2.Contrast.LLVec]+2*nParams;
Delay_AICMat = -2.* [fits.model1.Delay.LLVec; fits.model2.Delay.LLVec]+2*nParams;
modcompp = Contrast_AICMat + Delay_AICMat;
nSubj = size(modcompp,2);


switch comparetype
    case 'Contrast' % just contrast
        comparevalue = Contrast_AICMat;
    case 'Delay' % just delay
        comparevalue = Delay_AICMat;
    case 'both' % both
        comparevalue = modcompp;
end

% positive number means
moddiff = bsxfun(@minus,comparevalue,comparevalue(1,:))
mean_moddiff = mean(moddiff,2)
sem_moddiff = std(moddiff,[],2)./sqrt(nSubj)
figure;
bar(moddiff')


%% mdoel comparison

% clear all
exptype = 'Delay';
comptype = 'BIC';

filepath = 'analysis/2_fitdata/fits/';
load([filepath 'modelfits.mat'],'fits')

nLLMat = [-fits.model1.(exptype).LLVec'  -fits.model2.(exptype).LLVec'];
nSubj = size(nLLMat,1);
nParamsVec = [4 4]; % both models have 4 parameters
nTrialsVec = 1440*ones(1,nSubj);

[AIC, BIC, AICc] = modcomp(nLLMat,nParamsVec,nTrialsVec);

switch comptype
    case 'nLL'
        val1 = nLLMat;
    case 'AIC'
        val1 = AIC;
    case 'AICc'
        val1 = AICc;
    case 'BIC'
        val1 = BIC;
end

moddiff = bsxfun(@minus,val1,val1(:,1))
mean_moddiff = mean(moddiff)
sem_moddiff = std(moddiff)./sqrt(nSubj)

% bar(moddiff(:,2))


% colors
terracotta = 'k';% aspencolors('terracotta');
green = 'k';%aspencolors('leafgreen');

figure;
bar(sort(moddiff(:,2)),'k');
ylim([-25 50])
defaultplot
set(gca,'YTick',-25:25:50)
% blah = moddiff(:,2);%sort(moddiff(:,2));
% idx1 = blah < 0;
% listt = 1:nSubj;
% 
% bar(listt(idx1),blah(idx1),'EdgeColor','none','FaceColor',terracotta);
% hold on;
% bar(listt(~idx1),blah(~idx1),'EdgeColor','none','FaceColor',green);
% set(gca,'XTick',[])
% defaultplot

%% bayesian model selection

BMSMat = -BIC./2;

[alpha2, expr2, xp2, pxp2, bor2] = spm_BMS(BMSMat)
% pmodel = expr*nums; % number of subjects from each model


%% expr barplot

a01 = sum(alpha1);
exprSD1 = sqrt(alpha1.*(a01-alpha1)./(a01^2 .*(a01+1)));

a02 = sum(alpha2);
exprSD2 = sqrt(alpha2.*(a02-alpha2)./(a02^2 .*(a02+1)));

figure
axis([0.5 2.5 0 1])
bar([2 1],expr1,'FaceColor','none');
hold on; 
errorb([2 1],expr1,exprSD1)
set(gca,'XTick', 1:2, 'XTickLabel',{'max','optimal'},'YTick',0:0.5:1)
axis([0.5 2.5 0 1])
defaultplot

figure
axis([0.5 2.5 0 1])
bar([2 1],expr2,'FaceColor','none');
hold on; 
errorb([2 1],expr2,exprSD2)
set(gca,'XTick', 1:2, 'XTickLabel',{'max','optimal'},'YTick',0:0.5:1)
axis([0.5 2.5 0 1])
defaultplot


%% pxp
figure
bar([2 1],pxp1,'FaceColor','none');%,'.','MarkerSize',18')
axis([0.5 2.5 0 1])
set(gca,'XTick', 1:2, 'XTickLabel',{'max','optimal'},'YTick',0:0.5:1)
defaultplot

figure
bar([2 1],pxp2,'FaceColor','none');%,'.','MarkerSize',18')
axis([0.5 2.5 0 1])
set(gca,'XTick', 1:2, 'XTickLabel',{'max','optimal'},'YTick',0:0.5:1)
defaultplot

%% =======================
%  PLOTS
% =======================


%% indvl subject psychometric function

clear all
subjids = {'1','2','3','4','5','6','7','8','9'};
exptype = 'Delay';
condVec = [1 4];
nSubj = length(subjids);

for isubj = 1:nSubj;
    
    figure;
    plot_psychometricfunction(subjids(isubj),exptype,condVec); % plot real data
    pause
end

%% VSS MODEL FIT PLOTS

clear all

filepath = 'analysis/2_fitdata/fits/';
load([filepath 'modelfits.mat'])
subjid = {'7'};
exptype = 'Contrast';

deltaVec = 2.5:5:87.5;
colorMat = [aspencolors('terracotta');aspencolors('leafgreen')];

% condlegendd = {'0/2', '2/1','2/2','4/1'};
subplotVec = [3 1 2 4];
figure;
[ha, pos] = tight_subplot(2,2,.05,[.1 .03],[.1 .03]);
for imodel = 1:2;
    bfp = fits.(['model' num2str(imodel)]).(exptype).bfpMat(str2double(subjid),:);
    try
        if bfp(4) == 0
            bfp(4) = 5e-3;
        end
    end
    
    X = []; % data in luigiform
    Nsamples = 1000; % default if empty
    
    [~,prmat,X] = AhyBCL_datalikeall(bfp,X,imodel,Nsamples);
    
    % plot psychometric function
    
    hold on;
    colors = colorMat(imodel,:);
    for icond = 1:4
        %         subplot(2,2,subplotVec(icond))
        axes(ha(subplotVec(icond)))
        pc = prmat{icond}(:,1)';
        pc_std = pc.*(1-pc)./24;
        plot_summaryfit(deltaVec,[],[],pc,pc_std,colors,colors);
        axis([0 90 0 1]);
        ax = ha(subplotVec(icond));
        ax.XTick = [0 30 60 90];
        ax.YTick = 0:.5:1;
        if any(subplotVec(icond) == [2 4]); ax.YTickLabel = []; else ax.YTickLabel = 0:.5:1; end
        if any(subplotVec(icond) == [1 2]); ax.XTickLabel = []; else ax.XTickLabel = [0 30 60 90]; end
        ax.FontSize = 14;
        %     fill([deltaVec fliplr(deltaVec)],[pc-pc_std fliplr(pc + pc_std)],colors(icond,:))
    end
end

for icond = 1:4
    axes(ha(subplotVec(icond))); hold on
    %     subplot(2,2,icond)
    plot_psychometricfunction(subjid,exptype,icond); % plot real data
end

%% indvl subjects model fits and data plots

clear all

filepath = 'analysis/2_fitdata/fits/';
load([filepath 'modelfits.mat'])
% subjids = {'1','2','3','4','5','6'};
subjids = {'3'};
exptype = 'Delay';
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
    colors = aspencolors(4,'passport'); %['b'; 'y'; 'g'; 'r'];
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

%% average plot subplot

clear all

filepath = 'analysis/2_fitdata/fits/';
load([filepath 'modelfits.mat'])
subjids = {'1','2','3','4','5','6','7','8','9'};
% subjids = {'ALM','DR','EN','MR'};
exptype = 'Contrast';
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

%% indvl subjects model fits and data plots

clear all

filepath = 'analysis/2_fitdata/fits/';
load([filepath 'modelfits.mat'])
subjids = {'3'};
exptype = 'Delay';
model = 2;
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
    colors = aspencolors(4,'passport'); %['b'; 'y'; 'g'; 'r'];
    condlegendd = {'0/2', '2/1','2/2','4/1'};
    for icond = 1:4
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

%% average plot on two plots

clear all

filepath = 'analysis/2_fitdata/fits/';
load([filepath 'modelfits.mat'])
subjids = {'1','2','3','4','5','6','7','8','9'};
% subjids = {'ALM','DR','EN','MR'};
exptype = 'Contrast';
model = 2;
condVec = [3];
% condlegendd = {'0/2', '2/1','2/2','4/1'};

modelstr = sprintf('model%d',model);
nSubj = length(subjids);
deltaVec = 2.5:5:87.5;
nConds = length(condVec);

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
% figure;
% hold on;
colors = aspencolors(4,'passport');%['b'; 'y'; 'g'; 'r'];


[ha, pos] = tight_subplot(1,2,.05,[.1 .03],[.1 .03]);
for icond = 1:nConds
    cond = condVec(icond);
    
    if mod(cond,3) == 1; axes(ha(1)); ax = ha(1); else axes(ha(2));ax = ha(2); end
    
    plot_psychometricfunction(subjids,exptype,cond); % plot real data
    plot_summaryfit(deltaVec,[],[],mean_mod(cond,:),sem_mod(cond,:),colors(cond,:),colors(cond,:));
    %     fill([deltaVec fliplr(deltaVec)],[pc-pc_std fliplr(pc + pc_std)],colors(icond,:))
    axis([0 90 0 1])

ax.XTick = [0 30 60 90];
ax.XTickLabel = [0 30 60 90];
ax.YTick = 0:.5:1;
ax.FontSize = 14;
end
ha(1).YTickLabel = 0:.5:1;

%% lines of PC w/in and across participants

clear all
subjidVec = {'1','2','3','4','5','6','7','8','9'};
exptype = 'Delay';
nSubj = length(subjidVec);

clf; hold on
axis([0.5 4.5 .3 .8])
condVec = [4 2 3 1];
colorMat = aspencolors(4,'passport');

for isubj = 1:nSubj;
    subj = subjidVec{isubj};
    
    filepath = 'experiment_data/output_mat/';
    load(sprintf('%sprocesseddata_ChangeLocalization_%s_subj%s.mat',filepath,exptype,subj));
    
    for icond = 1:4;
        PC(isubj,icond) = sum(data{icond}.Rmat(:,1))/sum(data{icond}.Rmat(:));
    end
    
    plot(1:4,PC(isubj,condVec),'Color',0.7*ones(1,3))
    for icond = 1:4;
        cond = condVec(icond);
        plot(icond,PC(isubj,cond),'Color',colorMat(cond,:),'Marker','.','MarkerSize',18)
    end
    
end

meanPC = mean(PC);

plot(1:4,meanPC(condVec),'k','LineWidth',3);
for icond = 1:4;
    cond = condVec(icond);
    plot(icond,meanPC(cond),'Color',colorMat(cond,:),'Marker','.','MarkerSize',36)
end

defaultplot;
set(gca,'XTick',[],'YTick',0.3:.1:.8)

%% errorbar of PC across participants

clear all
subjidVec = {'1','2','3','4','5','6','7','8','9'};
exptype = 'Delay';
nSubj = length(subjidVec);

for isubj = 1:nSubj;
    subj = subjidVec{isubj};
    
    filepath = 'experiment_data/output_mat/';
    load(sprintf('%sprocesseddata_ChangeLocalization_%s_subj%s.mat',filepath,exptype,subj));
    
    for icond = 1:4;
        PC(isubj,icond) = sum(data{icond}.Rmat(:,1))/sum(data{icond}.Rmat(:));
    end
end

meanPC = mean(PC);
semPC = std(PC)./sqrt(nSubj);

figure;
axis([0.5 4.5 .4 .7])
condVec = [1 4 3 2];
colorMat = aspencolors(4,'passport');
for icond = 1:4;
    cond = condVec(icond);
    hold on; 
    %     bar(icond,meanPC(cond))
    errorb(icond,meanPC(cond),semPC(cond),'Color',colorMat(cond,:))
end
defaultplot;
set(gca,'XTick',[],'YTick',0.4:.1:.7)

%% plot just high and low condition data
clear all
subjidVec = {'1','2','3','4','5','6','7','8','9'};
plot_psychometricfunction(subjidVec,'Delay',[1 4]);
axis([0 90 0 1])

ax.XTick = [0 30 60 90];
ax.XTickLabel = [0 30 60 90];
ax.YTick = 0:.5:1;
ax.FontSize = 14;
%% plot psychometric functions of subjects according to getbestcontrast!

clear all

% subjnameVec = {'ALM','DR','MR','EN','EK','DC','JP','MP','TC'};
subjnameVec = {'EK','DC','JP','MP','TC'};
subjidVec = {'5','6','7','8','9'};
nSubj = length(subjidVec);

for isubj = 1:nSubj;
    subjname = subjnameVec{isubj};
    contrastVec(isubj) = calculate_bestcontrast(subjname);

end

figure;
[B,I] = sort(contrastVec);
for isubj = 1:nSubj;
    subplot(1,5,isubj);
    
    subjid = subjidVec(I(isubj));
    plot_psychometricfunction(subjid,'Delay',[1 4]);
end

%% how many trials to get reasonable estimates on homogeneous conditions

clear all
% ====== GENERATE TRUE DATA ======
% parameters
model = 1;
theta = [3 6 1 1e-5];
nSamples = 1e4;
nAngles = 18;

% generate fake data just homogeneous condition
X{1}.Nset1 = 0; X{1}.Nset2 = 4; X{1}.Cset = 2;
X{2}.Nset1 = 4; X{2}.Nset2 = 0; X{2}.Cset = 1;
for iCnd = 1:numel(X); X{iCnd}.Rmat = []; end
[loglike,prmat,X] = AhyBCL_datalikeall(theta,X,model,nSamples);

% generate noisy data based on this
nTrials = 20;
% for iX = 1:numel(X)
iX = 1;
    mat = binornd(nTrials,prmat{iX}(:,1));
    mat = [mat nTrials-mat zeros(nAngles,1)];
    X{iX}.Rmat = mat;
    X{2}.Rmat = mat;
% end

% ====== ESTIMATE PARAMETERS =======
runmax = 10;
runlist = 1:10;
[bfp, LLVec, completedruns] = fit_parameters(X,model,runlist,runmax,nSamples);
bfp

% get model prediction
bfpbest = bfp(LLVec == max(LLVec),:);
bfpbest = [bfpbest(1:2) log(bfpbest(3)) bfpbest(4)];
[ll,prmat_predict,~] = AhyBCL_datalikeall(bfpbest,X,model,nSamples);

% plot
figure; hold on;
plot(prmat{1}(:,1),'Color',[0 0 1])
plot(prmat{2}(:,1), 'r')
plot(prmat_predict{1}(:,1), 'Color',[.5 .5 1])
plot(prmat_predict{2}(:,1), 'Color', [1 .5 .5])

% this doesnt look like a good fit, so seeing if the predicted parameters
% actually give a reasonalbe LL compared to true value
% -- the results of this show that the actual log likelihood is much lower
% for the true parameter value, which is comforting! now have to see why
% parameters aren't being estimated correctly
iX = 1;
[loglike,~] = AhyBCL_datalike1(theta,X{iX}.Nset1,X{iX}.Nset2,X{iX}.Cset,X{iX}.Rmat,model,nSamples)
[loglike,~] = AhyBCL_datalike1(bfpbest,X{iX}.Nset1,X{iX}.Nset2,X{iX}.Cset,X{iX}.Rmat,model,nSamples)
