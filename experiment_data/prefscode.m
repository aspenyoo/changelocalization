function prefs = prefscode(exptype, subjid, sessionnum, nTrialsPerCond)

if nargin < 2; subjid = []; end
if nargin < 3; sessionnum = 0; end
if nargin < 4; nTrialsPerCond = 60; end
if isempty(nTrialsPerCond); nTrialsPerCond = 60; end

if nargin < 1; exptype = input('Delay/Contrast: ','s');end
if isempty(subjid); subjid = input('enter subject ID: ', 's'); end

subjid = upper(subjid);

%=========================================================================
% THINGS THAT ALL EXPERIMENTS SHOULD HAVE IN COMMON
% ========================================================================

% colors
red                   = [200 0 0];
green                 = [0 200 0];
white                 = 255;
lightgrey             = 200;  % foreground
grey                  = 128;  % background
black                 = 0;

% experiment timing (sec)
initTrialDur = 0;     % larger fixation cross that initializes
fix1Dur = 0.5;          % fixation screen in the beginning
pres1Dur = 0.15;         % presentation 1
pres1ISIDur = 3;    % inter-stimulus interval w/in pres1 (only for seq)
pres2Dur = .100;        % will not exist when respInPres2 == 1
ITIDur = 0.5;           % inter-trial interval

% display settings
screenHeight    = 30;               % in cm (participant screen in curtislab)
bgColor         = grey;             % background color
stimColor       = lightgrey;        % stimulus color
fixLength       = 0.15;             % fixation cross length (dva)
fixColor        = black;            % fixation cross color
jitter          = 0.5;              % amount of x/y-jitter (dva)
stimecc         = 5;                % stimulus eccentricity (dva)
stimArea        = 2.25; %(1.5^2)    % stimulus area (deg).
stimType        = 'gabor';          % stimulus type

% ========================================================================
% EXPERIMENT SPECIFIC
% ========================================================================

% experiment file name
fileName = ['Exp_ChangeLocalization_' exptype '_subj' num2str(subjid)];


if strcmp(exptype(end-4:end),'Delay'); % detection task (yes/no)
    
    % yes(1)/no(0)
    vmprior = 8.742;             % kappa of vm prior. 0 if uniform prior. (used to be sd of gauss prior)
    %             if ~(vmprior); nTrialsPerCond = 10; end % 10 per method of constant stimuli
    screenshot = 0;         % do you want a screenshot?
    feedback = 0;           % do you want to give feedback on performance?
    allStimInPres2 = 1;     % all stimuli to be in pres2 (vs.just target)?
%     directionChange = 0;    % task clockwise/counterclockwise (vs. yes/no)?
    respInPres2 = 1;        % does 2nd presentation stay up until S responds?
    stimPresent = 0;        % are all stimuli presented simultaneous in first presentation?
    permLocInPres1 = 1;     % are stimuli locations in pres1 permuted?
    lineAtPres2 = 0;        % oriented line at presentation 2 (vs. a gabor or ellipse)
    fixationinITI =  0;     % if fixation cross shown during ITI
    %     pres2contrast = 0;      % randomly change contrast of presentation 2 gabors
    
    % experimental design
    deltaNum = [1]; %
    f1 = length(deltaNum);
    % if directionChange = 1, magnitude and direction of change
    % if directionChange = 0, should be 0 or 1 for yes/no change
    
    pres1stimNum =[0 2 2 2 2 4]; % stimulus set size in first presentation
    f2 = length(pres1stimNum);
    
    reliabilityNum = [ones(1,max(pres1stimNum))]; % contrast for gabor
    f3 = size(reliabilityNum,1);
    
    ISIdelayNum = [1];   % ISI delay time (sec)
    f4 = length(ISIdelayNum);
    
    % breaks and feedback
    blocknum = 1;           % number of blocks ( 1 + number of breaks )
    breakDuration = 20;     % duration of each break between blocks (sec)
    feedbacktrial = 50;%nTrialsPerCond*size(fullfact([f1 f2 f3 f4]),1);     % every feedbacktrialth trial, feedback will be given
    
elseif strcmp(exptype(end-7:end),'Contrast')
    
    % yes(1)/no(0)
    vmprior = 0;            % kappa of vm prior. 0 if uniform prior. (used to be sd of gauss prior)
    screenshot = 0;         % do you want a screenshot?
    feedback = 0;           % do you want to give feedback on performance?
    allStimInPres2 = 1;     % all stimuli to be in pres2 (vs.just target)?
%     directionChange = 0;    % task clockwise/counterclockwise (vs. yes/no)?
    respInPres2 = 1;        % does 2nd presentation stay up until S responds?
    stimPresent = 1;        % are all stimuli presented simultaneous in first presentation?
    permLocInPres1 = 1;     % are stimuli locations in pres1 permuted?
    lineAtPres2 = 0;        % oriented line at presentation 2 (vs. a gabor or ellipse)
    fixationinITI =  0;     % if fixation cross shown during ITI
    %     pres2contrast = 0;      % randomly change contrast of presentation 2 gabors
    
    % breaks and feedback
    blocknum = 1;           % number of blocks ( 1 + number of breaks )
    breakDuration = 20;     % duration of each break between blocks (sec)
    feedbacktrial = 24;     % every feedbacktrialth trial, feedback will be given
    
    % figure out low contrast value
    try
        lowcontrastval = calculate_bestcontrast(subjid);
    catch
        YN = input('Could not calculate low contrast value for subject. Continue with 0.3? Y/N: ','s');
        if strcmp(upper(YN),'Y');
            lowcontrastval = 0.1;
        else
            error('need to specify low contrat value')
        end
    end
    
    % experimental design
    deltaNum = [1]; %
    f1 = length(deltaNum);
    % if directionChange = 1, magnitude and direction of change
    % if directionChange = 0, should be 0 or 1 for yes/no change
    
    pres1stimNum =[4]; % stimulus set size in first presentation (or total set size if no delay time manipulation)
    f2 = length(pres1stimNum);
    
    reliabilityNum = [1 1 1 1; ...
        1 1 lowcontrastval lowcontrastval; ...
        1 1 lowcontrastval lowcontrastval; ...
        lowcontrastval*ones(1,4)]; % contrasts of gabors
    %     reliabilityNum = [ones(1,4); 1 1 .2 .2; .2*ones(1,4)]; % contrast for gabor
    f3 = size(reliabilityNum,1);
    
    ISIdelayNum = [1];   % ISI delay time (sec)
    f4 = length(ISIdelayNum);
    
elseif strcmp(exptype,'sausagestim')
    
    % info for current experiment
    expName = ['Exp_ChangeDiscrim_' exptype '_subj' num2str(subjid)];
    
    % yes(1)/no(0)
    gaussprior = 20;      % gaussian or stepwise
    screenshot = 0;      % do you want a screenshot?
    feedback = 1;        % do you want to give feedback on performance?
    allStimInPres2 = 0;  % all stimuli to be in pres2 (vs.just target)?
    directionChange = 1; % task clockwise/counterclockwise (vs. yes/no)?
    respInPres2 = 1;     % does 2nd presentation stay up until S responds?
    simultPres1 = 1;        % are all stimuli presented simultaneous in first presentation?
    permLocInPres1 = 0;  % are stimuli locations in pres1 permuted?
    
    
    % breaks and feedback
    blocknum = 1;     % number of blocks ( 1 + number of breaks )
    breakDuration = 20;   % duration of each break between blocks (sec)
    feedbacktrial = []; % every feedbacktrialth trial, feedback will be given
    
    % experimental design
    directionNum = [1]; % always changing
    f1 = length(directionNum);
    
    setSizeNum = [1]; % stimulus set size
    f2 = length(setSizeNum);
    
    % reliabilityName = ['low', 'medium', 'high'];
    reliabilityNum = [0 pi/6 pi/2]; % ellipse eccentricity [.65 .97];%
    f3 = length(reliabilityNum);
    
    ISIdelayNum = [0.981]; % [.5 1 2];   % ISI delay time (sec)
    f4 = length(ISIdelayNum);
    
    pres1stimNum = 1;
end


if strcmp(exptype(1:5),'Pract'); % if this is practice
    % info for current experiment
    fileName = ['Pract' fileName(4:end)]; % changing the name to have "pract" instead of "exp"
    feedbacktrial = 1000;     % no feedback!
end

if sum(pres1stimNum == 1)==1;
    stimecc = 0;
    jitter = 0;
end

% making struct
varNames = who;
fieldNames1 = ['fieldNames' varNames'];
prefs = v2struct(fieldNames1);


