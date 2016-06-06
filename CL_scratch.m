%% simulate data

model = 1; % 1: optimal. 2: fixed
X = []; % data in luigiform
Nsamples = []; % default if empty

theta = [3 10 1 0.01]; % [jbar1 jbar2 tau lapse]

[~,prmat,X] = AhyBCL_datalikeall(theta,X,model,Nsamples);

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
% X = [ones(size(designMat,1),1) X];
% [b,bint,~,~,stats] = regress(y,X) % stats: R^2, F, p, and estimate of error variance

% ordinal multinomial regression
[B,dev,stats] = mnrfit(X,y,'model','ordinal')
LL = stats.beta - 1.96*stats.se
UL = stats.beta + 1.96*stats.se

% plot regression fits
xx = linspace(5,85,50);
XX = [ones(50,1) 3*ones(50,1) xx'];
BMat = [B(1:3)'; repmat(B(4:end),1,3)];
logoddsMat = XX*BMat; % log odds for each of the 3 equations
logoddssign = sign(logoddsMat);
difflogoddssign = diff(logoddssign,1,2); % difference between the signs of log odds
[row,col] = find(difflogoddssign == 2);
predictedConf = nan(50,1);
predictedConf(row) = col+1;
predictedConf(logoddssign(:,1) == 1 ) = 1;
predictedConf(logoddssign(:,3) == -1) = 4;
hold on; plot(xx,predictedConf);


%% create model fit datat psychometric fn
clear all

load('modelfits.mat')
model = 1;
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


