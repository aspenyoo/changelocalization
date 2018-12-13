%% code related to Experiment 4: Orientation change Localizaton with Delay 
% Reliability manipulation

% Models: optimal and fixed
% model parameters: Jbar1 (low reliability condition), Jbar2, tau, lapse


bleh = '/Users/blobface/Google Drive File Stream/My Drive/Research/VSTM/Aspen Luigi - Reliability in VWM/Exp 4 Change localization delay reliability/data';
%% get data in a more intuitive organization
% a temporary fix, until the experiment outputs data as this is structured

% ideal designMat should be
% delta
% Nset1 (number of items in first primary presentation interval)
% Cset (primary presenation interval in which target lives)
% correct
% RT

% ideal stimuliMat should be
% item 1-4 interval interval (4 columns) (organized by visual quadrant)
% item orientations (4 columns) (organized by vis quadrant)
% item x and y coordinates (8 columns) (organized by vis quadrant)
% quadrant of target (1 column)

clear all
subjid = 'AA';
nTrials = 360;

load(sprintf('Exp_ChangeLocalization_Delay_subj%s.mat',subjid),...
    'designMat','stimuliMat','names')

% ============ EDITING STIMULIMAT =======

% indices of important information
idx_S_xcoords = 9:12; % x coordinates of items
idx_S_ycoords = 13:16; % y coordinates of items

% other useful calculations
screencenter_x = mean(mean(stimuliMat(:,idx_S_xcoords)));
screencenter_y = mean(mean(stimuliMat(:,idx_S_ycoords)));
stimuliMat(:,idx_S_xcoords) = stimuliMat(:,idx_S_xcoords) - screencenter_x;
stimuliMat(:,idx_S_ycoords) = stimuliMat(:,idx_S_ycoords) - screencenter_y;

% [sMat_interval, sMat_orientation, sMat_xcoord, sMat_ycoord] = deal(nan(nTrials,4));
sMat_interval = stimuliMat(:,1:4)';
sMat_orientation = stimuliMat(:,5:8)';
sMat_xcoord = stimuliMat(:,idx_S_xcoords)';
sMat_ycoord = stimuliMat(:,idx_S_ycoords)';

% creating sMat, which is stimuliMat reorganized so items are organized by
% quadrant, insead of how it was
sMat = nan(nTrials,17);
idxx = zeros(4,nTrials);
for iquad = 1:4
    switch iquad
        case 1
            idx = (sMat_xcoord > 0) & (sMat_ycoord > 0);
        case 2
            idx = (sMat_xcoord < 0) & (sMat_ycoord > 0);
        case 3
            idx = (sMat_xcoord < 0) & (sMat_ycoord < 0);
        case 4
            idx = (sMat_xcoord > 0) & (sMat_ycoord < 0);
    end
    idxx = idxx + iquad.*idx;
    
    sMat(:,iquad) = sMat_interval(idx);
    sMat(:,4+iquad) = sMat_orientation(idx);
    sMat(:,8+iquad) = sMat_xcoord(idx);
    sMat(:,12+iquad) = sMat_ycoord(idx);
end
sMat(:,17) = stimuliMat(:,17);


% =============== EDITING DESIGNMAT =========
li = sub2ind([nTrials 4],1:nTrials,sMat(:,17)');
sMat_interval = sMat(:,1:4);

dMat = nan(nTrials,5);
dMat(:,1) = designMat(:,1);                 % delta
dMat(:,2) = sum(sMat(:,1:4) == 1,2);        % Nset1 (number of items in first primary presentation interval)
dMat(:,3) = sMat_interval(li);              % Cset (primary presenation interval in which target lives)
dMat(:,end-1:end) = designMat(:,end-1:end); % correct & RT

% ===== names ======
names.designMat = {'delta','Nset1','Cset','correct?','RT'};
names.stimuliMat = {'col 1-4: presentation interval','col 4-8: orientation presentation',...
    'col 9-12: item x coordinates', 'col 13-16: item y coordinates', 'col 17: target quadrant number'};

designMat = dMat;
stimuliMat = sMat;

save(sprintf('%s/Exp_ChangeLocalization_Delay_subj%s_new.mat',bleh,subjid),...
    'designMat','stimuliMat','names')

%% get data in simple, analyzable organization

clear all

subjid = 'AA';
load(sprintf('Exp_ChangeLocalization_Delay_subj%s_new.mat',subjid),...
    'designMat','stimuliMat','names'); % load data file

% conditions
X{1}.Nset1 = 0; X{1}.Nset2 = 4; X{1}.Cset = 2; % high-reliability
X{2}.Nset1 = 4; X{2}.Nset2 = 0; X{2}.Cset = 1; % low-reliability
X{3}.Nset1 = 2; X{3}.Nset2 = 2; X{3}.Cset = 2; % mixed-rel. change in high
X{4}.Nset1 = 2; X{4}.Nset2 = 2; X{4}.Cset = 1; % mixed-rel. change in low
nCond = length(X);

% indices for information stuff
idx_D_delta = 1; % column index of amount of change in degrees
idx_D_correct = 10; % column index of whetheer correct or not
idx_S_interval = 1:4; % columns showing intervals each item is presented in
idx_S_quadrant = 17; % quadrant of target
idx_S_xcoords = 9:12; % x coordinates of items
idx_S_ycoords = 13:16; % y coordinates of items

% get data in simpler format
data = cell(1,nCond);
for icond = 1:nCond
    
    idx = (X{icond}.Nset1 == designMat(:,2));       % correct number of low-reliability trials
    idx = idx & (X{icond}.Cset == designMat(:,3));  % target in correct interval
    
    data{icond} = designMat(idx,[1 4 5]);
end


%% plot cumulative correct responses (NOT SAME AS PSYCHOMETRIC FUNCTION)

close all
for icond = 1:nCond
    d = data{icond}(:,1:2);
    d(:,1) = abs(d(:,1));
    d = sortrows(d,1);
    d(:,2) = cumsum(d(:,2));
    d(:,2) = d(:,2)./(d(end,2));
    
    plot(d(:,1),d(:,2),'o-')
    hold on;
    
    labels{icond} = sprintf('Nset %d, Cset, %d', X{icond}.Nset1, X{icond}.Cset);
    
end

legend(labels)
defaultplot

%% psychometric function

nQuants = 6;
close all
for icond = 1:nCond
    d = data{icond}(:,1:2);
    d(:,1) = abs(d(:,1));
    d = sortrows(d,1);
    
    ntrials = size(d,1);
    binedges = round(linspace(1,ntrials,nQuants+1));
    
    dd = nan(nQuants,2);
    for iquant = 1:nQuants
        dd(iquant,:) = mean(d(binedges(iquant):binedges(iquant+1),:));
    end
    
    plot(dd(:,1),dd(:,2),'o-')
    hold on;
    
    labels{icond} = sprintf('Nset %d, Cset, %d', X{icond}.Nset1, X{icond}.Cset);
end

legend(labels)
defaultplot

%% parameter recovery simulations

clear all;
vmprior = 8.742;
nTrials = 60;

delta = circ_vmrnd(0,vmprior, [nTrials,1])*180/pi;

%%

clear all
% ====== GENERATE TRUE DATA ======
% parameters
model = 1;
theta = [3 5 2 0.1];
nSamples = 1e4;
ds = (5:10:85)';
% ds = (2.5:5:87.5)';
nAngles = length(ds);

% generate fake data just homogeneous condition
X{1}.Nset1 = 0; X{1}.Nset2 = 4; X{1}.Cset = 2;
X{2}.Nset1 = 4; X{2}.Nset2 = 0; X{2}.Cset = 1;
for iCnd = 1:numel(X); X{iCnd}.Rmat = []; end
[loglike,prmat,X] = AhyBCL_datalikeall(theta,X,model,nSamples,ds);

% generate noisy data based on this
nTrials = 20;
for iX = 1:numel(X)
% iX = 1;
    mat = binornd(nTrials,prmat{iX}(:,1));
    mat = [mat nTrials-mat zeros(nAngles,1)];
    X{iX}.Rmat = mat;
    X{2}.Rmat = mat;
end

% ====== ESTIMATE PARAMETERS =======
runmax = 50;
runlist = 1:50;
[bfp, LLVec, completedruns] = fit_parameters(X,model,runlist,runmax,nSamples,ds);

% get model prediction
bfpbest = bfp(LLVec == max(LLVec),:);
bfpbest = [exp(bfpbest(1:3)) bfpbest(4)];
bfpbest
[ll,prmat_predict,~] = AhyBCL_datalikeall(bfpbest,X,model,nSamples,ds);

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
[loglike,~] = AhyBCL_datalike1(theta,X{iX}.Nset1,X{iX}.Nset2,X{iX}.Cset,X{iX}.Rmat,model,nSamples,ds)
[loglike,~] = AhyBCL_datalike1(bfpbest,X{iX}.Nset1,X{iX}.Nset2,X{iX}.Cset,X{iX}.Rmat,model,nSamples,ds)
