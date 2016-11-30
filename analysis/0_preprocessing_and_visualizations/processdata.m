function processdata(subjid,exptype)

subjid = upper(subjid);
filepath = 'experiment_data/output_mat/';
filename = ['Exp_ChangeLocalization_' exptype '_subj' subjid '.mat'];

load([filepath filename],'designMat','stimuliMat','names')

% deleting incompleted trials
stimuliMat(isnan(designMat(:,end-2)),:) = [];
designMat(isnan(designMat(:,end-2)),:) = [];

% getting column indices
targetlocationidx = 17;

deltaidx = 1;
correctidx = 10;
responseidx = 9;


switch exptype
    case 'Delay'
        % Nset1,Nset2,Cset,Rmat
        nLowRelVec = [0 2 2 4];            % number of stimuli in low rel condition
        targetCondVec = [2 1 2 1];            % which reliability condition the target was in (1: low rel. 2: high rel)
        respVec = [3 0 4 0 0 0 2 0 1];  % go from response key to location index
        reliabilitymanipulationidx = 2;
        
    case 'Contrast'
        % Nset1,Nset2,Cset,Rmat
        reliabilitymanipulationidx = 3;
        contrastVec = unique(designMat(:,reliabilitymanipulationidx))';
        nLowRelVec = [0 2 2 4];            % number of stimuli in low rel condition
        targetCondVec = [fliplr(contrastVec) fliplr(contrastVec)];           % which reliability condition the target was in
        respVec = [3 0 4 0 0 0 2 0 1];  % go from response key to location index
end
setsize = 4;                    % total stimuli presented
deltaVec = 2.5:5:87.5;              % possible change values
nCond = length(nLowRelVec);        % number of conditions

% for each condition, make the whatever
data = cell(1,nCond);
for icond = 1:nCond
    nLowRel = nLowRelVec(icond);
    targetCond = targetCondVec(icond);
   
    data{icond}.Nset1 = nLowRel;             % set size for low reliability condition
    data{icond}.Nset2 = setsize - nLowRel;   % set size for high reliability condition
    data{icond}.Cset = find(contrastVec == targetCond);               % which display change target was presented
    
    % finding indices of trials in the correct condition (reliability and target location)
    switch exptype
        case 'Delay'
            idx = (designMat(:,reliabilitymanipulationidx) == nLowRel);
        case 'Contrast'
            idx = (sum(designMat(:,reliabilitymanipulationidx:reliabilitymanipulationidx+3) == contrastVec(1),2) == nLowRel);
    end
    % getting the condition number of the target
    idx = idx & (diag(stimuliMat(:,stimuliMat(:,targetlocationidx))) == targetCond);
    
    Rmat = nan(length(deltaVec),3);
    % for each change value, delta...
    for idelta = 1:length(deltaVec)
        delta = deltaVec(idelta); % current delta value
        
        idxx = idx & (abs(designMat(:,deltaidx)) == delta); % how many of that delta
        Rmat(idelta,1) = sum(designMat(idxx,correctidx) == 1); % how many correct
        Rmat(idelta,2) = sum((designMat(idxx,correctidx) == 0) & ... % how many incorrect & 
            (diag(stimuliMat(idxx,respVec(designMat(idxx,responseidx)))) == targetCond)); % same rel as correct response
        Rmat(idelta,3) = sum((designMat(idxx,correctidx) == 0) & ... % how many incorrect & 
            (diag(stimuliMat(idxx,respVec(designMat(idxx,responseidx)))) ~= targetCond)); % different rel as correct response
    end
    
    data{icond}.Rmat = Rmat;
end

% save
save(sprintf('%s/processeddata_ChangeLocalization_%s_subj%s.mat', filepath,exptype,subjid),'data');
display('done')