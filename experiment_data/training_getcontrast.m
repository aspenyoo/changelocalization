function training_getcontrast(subjid)

% ==================================================
%   INITIAL EXPERIMENTAL STUFF
% ==================================================

changeVec = [11.25 33.75 56.25 78.75];
changeVecNames = {'cond1','cond2','cond3','cond4'};
nCond = length(changeVec);
nTrialsperCond = 150;
nTrials = nCond*nTrialsperCond;
names.designMat = {'amnt change','logconstrast1','logcontrast2','target loc','resp','RT','correct?'};

designMat = nan(nTrials,length(names.designMat));  % all the experimental stuff
temp = [repmat(changeVec,1,nTrialsperCond/2) -repmat(changeVec,1,nTrialsperCond)]; % half are negative
designMat(:,1) = temp(randperm(nTrials));  % pseudorandom shuffling
designMat(:,4) = ceil(rand(nTrials,1)*4);

% initializing stuff for each condition
psy = [];
method = 'ent';     % Minimize the expected posterior entropy
vars = [1 1 1];     % Minimize posterior entropy of mean and sigma
for icond = 1:nCond;
    condstring = changeVecNames{icond};  % condition string used for structs
    
    % get list of trials in that condiiton
    trials.(condstring) = find(abs(designMat(:,1)) == changeVec(icond));
    
    % ========== initialize psy struct for each condition ==========
    % Initialize PSY structure
    psy.(condstring) = [];
    
    % Set chance level (for PCORRECT psychometric functions)
    psy.(condstring).gamma = 0.25;
    
    % Specify user-defined psychometric function (as a string)
    psy.(condstring).psychofun{1} = '@(x,mu,sigma,lambda,gamma) psyfun_pcorrect(x,mu,sigma,lambda,gamma,@psygumbelcdf);';
    psy.(condstring).psychofun{2} = '@(x,mu,sigma,lambda,gamma) psyfun_pcorrect(x,mu,sigma,lambda,gamma,@psynormcdf);';
    psy.(condstring).psychofun{3} = '@(x,mu,sigma,lambda,gamma) psyfun_pcorrect(x,mu,sigma,lambda,gamma,@psylogicdf);';
    
    % Define range for stimulus and for parameters of the psychometric function
    % (lower bound, upper bound, number of points)
    MinContrast = 1;     % Minimum contrast
    psy.(condstring).range.x = [log(MinContrast),log(100),21];      % Stimulus range in log Hz
    psy.(condstring).range.mu = [log(MinContrast),log(50),31];     % Psychometric function mean in log Hz
    psy.(condstring).range.sigma = [0.1,4,19];                   % The range for sigma is automatically converted to log spacing
    psy.(condstring).range.lambda = [0,0.5,21];
    % psy.range.lambda = [0.05-eps,0.05+eps,2];  % This would fix the lapse rate to 0.05
    
    % Define priors over parameters
    psy.(condstring).priors.mu = [log(MinContrast),4];   % mean and std of (truncated) Gaussian prior over MU
    psy.(condstring).priors.logsigma = [0,1];        % mean and std of (truncated) Gaussian prior over log SIGMA (Inf std means flat prior)
    psy.(condstring).priors.lambda = [1.3 4];         % alpha and beta parameters of beta pdf over LAMBDA
    
    % Units -- used just for plotting in axis labels and titles
    psy.(condstring).units.x = 'log';
    psy.(condstring).units.mu = 'log';
    psy.(condstring).units.sigma = 'log';
    psy.(condstring).units.lambda = [];
    psy.(condstring).units.psychofun = {'Gumbel','Normal','Logistic'};
    
    % Refractory time before presenting same stimulus again
    psy.(condstring).reftime = 2;        % Expected number of trials (geometric distribution)
    psy.(condstring).refradius = 0;      % Refractory radius around stimulus (in x units)
    
    
    % initialize random contrast for first trial
    [x, psy.(condstring)] = psybayes(psy.(condstring), method, vars, [], []);
    designMat(trials.(condstring)(1),2) = x;
    trials.(condstring)(1) = [];
end

% ========================================================
%    DO EXPERIMENT!!
% ========================================================

% initial stimulus orientations
names.stimuliMat = {'orientation1','orientation2','orientation3','orientation4','pos1_x','pos2_x','pos3_x','pos4_x','pos1_y','pos2_y','pos3_y','pos4_y'};
setsize = length(names.stimuliMat)/3;
stimuliMat = nan(nTrials,setsize*3);
stimuliMat(:,1:setsize) = round(rand(nTrials,setsize).*180);
stimuliMat(stimuliMat == 0,1:setsize) = 180;

% response buttons
keys = [KbName('9') KbName('7') KbName('1') KbName('3') KbName('esc')];
keysNum = [9 7 1 3 0];

% trial timing (seconds)
length_fixation = 0.5;
length_int1 = 0.1;
length_ISI = 1;
length_int2 = Inf;  % on until response
length_ITI = 0.5;

% display settings
screenHeight    = 30.5;               % in cm (Dell@T115A: ~48cm; Dell@T101C: ~40 cm)
screenDistance  = 40;
bgColor         = 128;              % background color
stimColor       = 200;              % stimulus color
fixLength       = 0.2;                % fixation cross length (DVA)
fixColor        = 0;                % fixation cross color
jitter          = 0.5;                % amount of x/y-jitter (dva)
stimecc         = 5;                % stimulus eccentricity (dva)
stimArea        = 2.25; %(1.5^2)    % stimulus area (deg).

% file saving stuff
if nargin < 1; subjid = input('enter subject ID: ', 's'); end
subjid = upper(subjid);
fidmat = fullfile('output_mat',['Training_ChangeLocalization_' subjid '.mat']);

% screen info (visual)
screenNumber =  max(Screen('Screens'));       % use external screen if exists
[w, h] = Screen('WindowSize', screenNumber);  % screen resolution of display
screenResolution = [w h];                 % screen resolution
screenCenter = screenResolution/2;       % screen center
screenAngle = atand((screenHeight/2) / screenDistance) ; % total visual angle of screen (height)
screen_ppd = h / screenAngle;  % pixels per degree
stimArea = stimArea * screen_ppd^2; % stimulus area (pixels)
stimLength = round(sqrt(stimArea));
srcrect = [0 0 stimLength stimLength];
fixLength = round(fixLength*screen_ppd);

% open screen
% if isempty(sessionnum)
%     windowPtr = Screen('OpenWindow',screenNumber,prefs.grey,[],32,2);
% else
    windowPtr = 10;
% end
Screen(windowPtr,'BlendFunction',GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
HideCursor;

% positions of stimuli
locAngles = pi/4-(1:setsize)*(2*pi)/setsize;
[X_pos, Y_pos] = pol2cart(locAngles, screen_ppd * stimecc);
X_pos = X_pos + screenCenter(1);
Y_pos = Y_pos + screenCenter(2);

dy = 30;
textx = w/2-300;

% ====== DETECTION INFO SCREEN =======
texty = screenCenter(2) - 200;
Screen('TextSize',windowPtr,20);
Screen('TextFont',windowPtr,'Helvetica');
Screen('DrawText',windowPtr,'Indicate which orientation changed.',textx,texty,[255 255 255]);
Screen('Flip', windowPtr);
pause
% waitForKey;


for itrial = 1:nTrials
    
    orientationchange = designMat(itrial,1);
    condstring = changeVecNames{changeVec == abs(orientationchange)}; % current condition
    x = designMat(itrial,2); % current level
    
    % jitter positions
    x_pos = X_pos + round((rand(1,setsize)-.5)*jitter*screen_ppd);
    y_pos = Y_pos + round((rand(1,setsize)-.5)*jitter*screen_ppd);
    stimuliMat(itrial,setsize+1:2*setsize) = x_pos;
    stimuliMat(itrial,2*setsize+1:3*setsize) = y_pos;
    
    % current trial stuff
    pres1orientations = stimuliMat(itrial,1:setsize);
    pres2orientations = pres1orientations;
    locChange = designMat(itrial,4);
    pres2orientations(locChange) = pres1orientations(locChange) + designMat(itrial,1);
    
    % contrast of second presentation
    if rand > 0.5;
        designMat(itrial,3) = x; % same contrast
    else
        designMat(itrial,3) = log(100); % full contrast
    end
    
    % ========================================
    % ================ TRIAL =================
    % ========================================
    
    % ------------ fixation ----------------
    Screen('fillRect',windowPtr,bgColor);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),fixColor,fixLength);
    t0 = GetSecs();
    % inittirladur = toc
    Screen('flip',windowPtr); % tic;
    while (GetSecs()-t0) < length_fixation;
        % do nothing
    end
    
    % ----------------- int 1 --------------------
    Screen('fillRect',windowPtr,bgColor);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),fixColor,fixLength);
    t0 = GetSecs();
    
    % making stimulus
    im = makeGabor(stimLength, exp(x)/100);
    stimpatch = Screen('MakeTexture',windowPtr,im);
    
    % drawing stimuli9
    for istim = 1:setsize;
        destrect = CenterRectOnPoint(srcrect,x_pos(istim),y_pos(istim));
        Screen('DrawTexture', windowPtr, stimpatch, srcrect,destrect, 180-pres1orientations(istim));
    end
    
    Screen('flip',windowPtr); % tic;
    while (GetSecs()-t0) < length_int1;
        % do nothing
    end
    
    % --------------- ISI ----------------------
    Screen('fillRect',windowPtr,bgColor);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),fixColor,fixLength);
    t0 = GetSecs();
    % inittirladur = toc
    Screen('flip',windowPtr); % tic;
    while (GetSecs()-t0) < length_ISI;
        % do nothing
    end
    
    % ----------------- int 2 --------------------
    Screen('fillRect',windowPtr,bgColor);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),fixColor,fixLength);
    t0 = GetSecs();
    
    % making stimulus
    im = makeGabor(stimLength, exp(designMat(itrial,3))/100);
    stimpatch = Screen('MakeTexture',windowPtr,im);
    
    % drawing stimuli
    for istim = 1:setsize;
        destrect = CenterRectOnPoint(srcrect,x_pos(istim),y_pos(istim));
        Screen('DrawTexture', windowPtr, stimpatch, srcrect,destrect, 180-pres2orientations(istim));
    end
    
    Screen('flip',windowPtr); % tic;
    if ~isinf(length_int2)
        while (GetSecs()-t0) < length_int2;
            % do nothing
        end
    end
    
    % ----------------- response -------------------
    [pressedKey, designMat(itrial,6)] = waitForKeys(keys,GetSecs());
    if pressedKey == length(keys);
        sca;
        ShowCursor;
        fclose('all');
        clear all;
    else
        designMat(itrial,5) = keysNum(pressedKey);
    end
    
    r = keysNum(pressedKey) == keysNum(designMat(itrial,4));
    designMat(itrial,7) = r; 
    
    % ------------------ ITI ----------------------
    Screen('fillRect',windowPtr,bgColor);
    t0 = GetSecs();
    % inittirladur = toc
    Screen('flip',windowPtr); % tic;
    while (GetSecs()-t0) < length_ITI;
        % do nothing
    end
    
    % ------------------------------------------------------
    % get next contrast for this condition
    [x, psy.(condstring)] = psybayes(psy.(condstring), method, vars, x, r);
    designMat(trials.(condstring)(1),2) = x;
    try trials.(condstring)(1) = []; end
    
    % save current stuff
    save(fidmat,'designMat','stimuliMat','psy');
end

end

function [pressedKey, RT] = waitForKeys(keys, tstart)

pressedKey=0;
while (1)
    
    [keyIsDown, secs, keyCode] = KbCheck();
    if  any(keyCode(keys))
        RT = GetSecs - tstart;
        pressedKey = find(keyCode(keys));
        break;
    end
    
end
end

function drawfixation(windowPtr,x,y,color,sizee)
Screen('DrawLine',windowPtr,color,x-sizee,y,x+sizee,y,2);
Screen('DrawLine',windowPtr,color,x,y-sizee,x,y+sizee,2);
end

