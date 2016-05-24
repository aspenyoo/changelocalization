function processdata_luigi(subjid)

subjid = upper(subjid);
projectpath('CL')
filepath = 'Experiment and Data/output_mat/';
filenamestart = ['Exp_ChangeDetection2_Seq_subj' subjid];
varnames = {'designMat','stimuliMat'};

[concatvars] = concatcode(filepath,filenamestart,varnames); % getting all the variables from the files

v2struct(concatvars); % taking the variables out of the struct

% deleting incompleted trials
stimuliMat(isnan(designMat(:,end-2)),:) = [];
designMat(isnan(designMat(:,end-2)),:) = [];


% Nset1,Nset2,Cset,Rmat
NsetVec = [0 2 2 4];            % number of stimuli in first display
CsetVec = [2 1 2 1];            % which display change target was presented
setsize = 4;                    % total stimuli presented
deltaVec = 5:5:90;              % possible change values
nCond = length(NsetVec);        % number of conditions
respVec = [3 0 4 0 0 0 2 0 1];  % go from response key to location index

% for each condition, make the whatever
data = cell(1,nCond);
for icond = 1:nCond;
    Nset1 = NsetVec(icond);
    Cset = CsetVec(icond);
   
    data{icond}.Nset1 = Nset1;             % set size for first display
    data{icond}.Nset2 = setsize - Nset1;   % set size 2nd display
    data{icond}.Cset = Cset;               % which display change target was presented
    
    % finding which trials have the matching Nset1 & Cset (Cset necessary for the two conditions with Nset = 2)
    idx = (designMat(:,2) == Nset1) & (diag(stimuliMat(:,stimuliMat(:,9))) == Cset);
    
    Rmat = nan(length(deltaVec),3);
    % for each change value, delta...
    for idelta = 1:length(deltaVec); 
        delta = deltaVec(idelta); % current delta value
        
        idxx = idx & (abs(designMat(:,1)) == delta); % how many of that delta
        Rmat(idelta,1) = sum(designMat(idxx,7) == 1); % how many correct
        Rmat(idelta,2) = sum((designMat(idxx,7) == 0) & ... % how many incorrect & 
            (diag(stimuliMat(idxx,respVec(designMat(idxx,6)))) == Cset)); % same display as correct response
        Rmat(idelta,3) = sum((designMat(idxx,7) == 0) & ... % how many incorrect & 
            (diag(stimuliMat(idxx,respVec(designMat(idxx,6)))) ~= Cset)); % different display as correct response
    end
    
    data{icond}.Rmat = Rmat;
end

% save
c = clock;
save(sprintf('data_ChangeLocalization_subj%s_%02d%02d%04d.mat', subjid,c(2),c(3),c(1)),'data');

