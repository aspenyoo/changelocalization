function prefs = prefscode(afctype, exptype, subjid, sessionnum, nTrialsPerCond)

if nargin < 3; subjid = []; end
if nargin < 4; sessionnum = 0; end
if nargin < 5; nTrialsPerCond = []; end

if nargin < 1; afctype = input('Discrim/Detect: ','s');end
if nargin < 2; exptype = input('Simult/Seq: ' , 's'); end
if isempty(subjid); subjid = input('enter subject ID: ', 's'); end
if isempty(nTrialsPerCond);
    if strcmp(afctype,'Discrim'); nTrialsPerCond = 30;end
    if strcmp(afctype,'Detect'); nTrialsPerCond = 1; end
end

subjid = upper(subjid);
dateStr = datestr(now, 30);

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
initTrialDur = 0.25;     % larger fixation cross that initializes 
fix1Dur = 0.5;          % fixation screen in the beginning
pres1Dur = 0.1;         % presentation 1
pres1ISIDur = 2;    % inter-stimulus interval w/in pres1 (only for seq)
pres2Dur = .100;        % will not exist when respInPres2 == 1
ITIDur = 0.5;           % inter-trial interval

% display settings
screenHeight    = 30.5;               % in cm (Dell@T115A: ~48cm; Dell@T101C: ~40 cm)
bgColor         = grey;             % background color
stimColor       = lightgrey;        % stimulus color
fixLength       = 7;                % fixation cross length
fixColor        = black;            % fixation cross color
jitter          = 0;                % amount of x/y-jitter (deg)
stimecc         = 5;                % stimulus eccentricity (deg)
stimArea        = 2.25; %(1.5^2)    % stimulus area (deg). 
stimType        = 'gabor';          % stimulus type

% ========================================================================
% EXPERIMENT SPECIFIC
% ========================================================================

if strcmp(afctype(end-5:end),'Detect'); % detection task (yes/no)
    
    % info for current experiment
    expName = ['Exp_ChangeLocalization_' exptype '_subj' num2str(subjid) '_session' num2str(sessionnum)];
    
    switch exptype
        case 'Simult' % simultaenous presentation
            
            error('not set up')
            
        case 'Seq' % sequential presentation
            
            % yes(1)/no(0)
            vmprior = 0;            % kappa of vm prior. 0 if uniform prior. (used to be sd of gauss prior)
%             if ~(vmprior); nTrialsPerCond = 10; end % 10 per method of constant stimuli
            screenshot = 0;         % do you want a screenshot?
            feedback = 0;           % do you want to give feedback on performance?
            allStimInPres2 = 1;     % all stimuli to be in pres2 (vs.just target)?
            directionChange = 0;    % task clockwise/counterclockwise (vs. yes/no)?
            respInPres2 = 1;        % does 2nd presentation stay up until S responds?
            stimPresent = 0;        % are all stimuli presented simultaneous in first presentation?
            permLocInPres1 = 1;     % are stimuli locations in pres1 permuted?
            
            % breaks and feedback
            blocknum = 6;           % number of blocks ( 1 + number of breaks )
            breakDuration = 20;     % duration of each break between blocks (sec)
            feedbacktrial = 24;     % every feedbacktrialth trial, feedback will be given
            
            % experimental design
            deltaNum = [1]; %
            f1 = length(deltaNum);
            % if directionChange = 1, magnitude and direction of change
            % if directionChange = 0, should be 0 or 1 for yes/no change
            
            pres1stimNum =[0 2 2 4]; % stimulus set size in first presentation
            f2 = length(pres1stimNum);
            
            reliabilityNum = [1]; % contrast for gabor
            f3 = length(reliabilityNum);
            
            ISIdelayNum = [1];   % ISI delay time (sec)
            f4 = length(ISIdelayNum);
    end
end


if strcmp(afctype(1:5),'Pract'); % detection task (yes/no)
    
    % info for current experiment
    expName = ['Pract' expName(4:end)]; % changing the name to have "pract" instead of "exp"
    
    feedback = 1;        % feedback
    pres1Dur = pres1Dur*2;      % twice as long ellipse presentation time
    blocknum = 1;     % number of blocks ( 1 + number of breaks )
end

if sum(pres1stimNum == 1)==1;
    stimecc = 0;
    jitter = 0;
end

% making struct
varNames = who;
fieldNames1 = ['fieldNames' varNames'];
prefs = v2struct(fieldNames1);


