function get_contrastresponsefunction(subjid)

% open screen
screenNumber = max(Screen('Screens'));       % use external screen if exists
[w, h] = Screen('WindowSize', screenNumber);  % screen resolution of smaller display
screenCenter = [w h]/2;    
windowPtr = Screen('OpenWindow',screenNumber,128*ones(1,3),[],32,2);
Screen(windowPtr,'BlendFunction',GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
Screen('TextSize',windowPtr,20);
Screen('TextFont',windowPtr,'Helvetica');
% HideCursor;

% =========== SETTING EXPERIMENTAL STUFF ===========
prefs = prefscode('Delay',subjid,2,1);
setsize = 4;
filename = ['output_mat/ContrastFn_' upper(subjid) '_'  datestr(now, 30) ...
    '.mat'];

% response keys
prefs.keys = [KbName('c') KbName('v') KbName('esc')];
prefs.keysNum = [1 -1 0];

screenDistance = 40;                      % distance between observer and screen (in cm)
screenAngle = atand((prefs.screenHeight/2) / screenDistance) ; % total visual angle of screen (height)
screen_ppd = h / screenAngle;  % pixels per degree
prefs.stimArea = prefs.stimArea * screen_ppd^2; % stimulus area (pixels)
stimLength = round(sqrt(prefs.stimArea));
StimSizes = [stimLength stimLength];

% stimulus locations (fixed at each quadrant)
locAngles = pi/4-(1:setsize)*(2*pi)/setsize;
[X, Y] = pol2cart(locAngles, screen_ppd * prefs.stimecc);
x_positions = X(1:setsize) + screenCenter(1)...
    + round((rand(1,setsize)-.5)*prefs.jitter*screen_ppd);
y_positions = Y(1:setsize) + screenCenter(2)...
    + round((rand(1,setsize)-.5)*prefs.jitter*screen_ppd);

% setting up design and stimulus matrices
nlapseTrials = 200;
nTrials = 250;
names.designMat = {'contrast','delta','correct'};
names.stimuliMat = {'stim1pres1','stim2pres1','stim3pres1','stim4pres1','target_position'};

lapsedesignMat = [ones(nlapseTrials,1) 30*ones(nlapseTrials,1) nan(nlapseTrials,1)];
lapsedesignMat(rand(nlapseTrials,1) > 0.5,2) = -30;
lapsestimuliMat = [round(rand(nlapseTrials,setsize)*180) randi(setsize,nlapseTrials,1)];

designMat = [nan(nTrials,1) 30*ones(nTrials,1) nan(nTrials,1)];
designMat(rand(nTrials,1) > 0.5,2) = -30;
stimuliMat = [round(rand(nTrials,setsize)*180) randi(setsize,nTrials,1)];

% ================ LINE STIMULI ===================
res = 1; % resolution
numDegrees = ceil(180/res);
% prefs.lineArea = prefs.stimArea; % so line is same area as ellipse
prefs.lineWidth = 3; % pixels
prefs.lineLength = stimLength;
lineCoordinates = nan(2,numDegrees);
for j = 1:numDegrees;
    [lineCoordinates(1,j), lineCoordinates(2,j)] = lineCoord(prefs.lineLength,j);
end

% ===========================================================
% =========== FULL CONTRAST (TO GET LAPSE RATE) ===========
% ===========================================================

% make gabor patch with full contrast
reliability = 1; % contrast of 1
im = makeGabor(stimLength, reliability);
StimPatch= Screen('MakeTexture',windowPtr,im);

for itrial = 1:nlapseTrials;
    
    % ========== FIXATION FOR BEGINNING OF TRIAL ==========
    Screen('fillRect',windowPtr,prefs.bgColor);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    t0 = GetSecs();
    % inittirladur = toc
    Screen('flip',windowPtr); % tic;
    while (GetSecs()-t0)<prefs.fix1Dur;
        % do nothing
    end
    
    % ========== PRESENT STIMULI FOR FIRST TIME ===========
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    pres1orientations = lapsestimuliMat(itrial,1:setsize);
    for istim = 1:setsize;
        srcrect = [0 0 StimSizes];
        destrect = CenterRectOnPoint(srcrect,x_positions(istim),y_positions(istim));
        Screen('DrawTexture', windowPtr, StimPatch, srcrect,destrect, 180-pres1orientations(istim));
    end
    Screen('flip',windowPtr); % tic;
    t0 = GetSecs();
    while (GetSecs()-t0)<prefs.pres1Dur;
        % do nothing
    end
    
    % ========== ISI ===========
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    Screen('flip',windowPtr); % tic;
    t0 = GetSecs();
    while (GetSecs()-t0)<prefs.ISIdelayNum;
        % do nothing
    end
    
    % ========== ONE PROBE SHOWS UP. RESPONSE ===========
    stim = lapsestimuliMat(itrial,end);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    if (prefs.lineAtPres2)
        pres2orientation = mod(90 - pres1orientations(stim) - lapsedesignMat(itrial,2),180);
        if (pres2orientation == 0); pres2orientation = 180; end
        lineStim = lineCoordinates(:,pres2orientation);
        xy = [-lineStim lineStim ];
        Screen('DrawLines',windowPtr, xy, prefs.lineWidth,prefs.stimColor,[x_positions(stim) y_positions(stim)],1);
    else
        pres2orientation = pres1orientations(stim) + lapsedesignMat(itrial,2);
        destrect = CenterRectOnPoint(srcrect,x_positions(stim),y_positions(stim));
        Screen('DrawTexture', windowPtr, StimPatch, srcrect,destrect, 180-pres2orientation);
    end
    

    Screen('flip',windowPtr); % tic;
    
    % check response
    [pressedKey] = waitForKeys(prefs.keys,GetSecs());
    if pressedKey == length(prefs.keys);
        sca;
        ShowCursor;
        fclose('all');
        clear all;
    else
        %         designMat(itrial,5+setsize) = prefs.keysNum(pressedKey);
    end
    
    % check whether it is correct or not
    if (sign(lapsedesignMat(itrial,2)/prefs.keysNum(pressedKey)) == 1) % if correct
        lapsedesignMat(itrial,3) = 1;
    else
        lapsedesignMat(itrial,3) = 0;
    end
    
    % ========== ITI ===========
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    Screen('flip',windowPtr); % tic;
    t0 = GetSecs();
    while (GetSecs()-t0)<prefs.ITIDur;
        % do nothingqs
    end
    
    % save file
    save(filename,'stimuliMat','designMat','lapsedesignMat','lapsestimuliMat','prefs','subjid')
end

lapse = 1-mean(lapsedesignMat(:,3))

% ===========================================
% INTERMEDIATE SCREEN
% ===========================================

textx = screenCenter(1) - 200;
texty = screenCenter(2) - 200;
Screen('TextSize',windowPtr,20);
Screen('TextFont',windowPtr,'Helvetica');
Screen('DrawText',windowPtr,'Now it is going to get a little harder!',textx,texty,[255 255 255]);
Screen('Flip', windowPtr);
pause;

% =====================================
% ADAPTIVE METHOD
% =====================================

% ========== INITIALIZE PSYBAYES ==========
psy = [];       % Initialize PSY structure

% Set chance level (for PCORRECT psychometric functions)
psy.gamma = 0.5;
% psyinit.gamma = [];   % Leave it empty for YES/NO psychometric functions

% Define range for stimulus and for parameters of the psychometric function
% (lower bound, upper bound, number of points)
psy.range.x = [5,100,61];
psy.range.mu = [15,40,51];
psy.range.sigma = [5,50,25];      % The range for sigma is automatically converted to log spacing
psy.range.lambda = [lapse,lapse,1];

% Define priors over parameters
psy.priors.mu = [23,5];                  % mean and std of (truncated) Gaussian prior over MU
psy.priors.logsigma = [log(20),Inf];   % mean and std of (truncated) Gaussian prior over log SIGMA (Inf std means flat prior)
psy.priors.lambda = [1 3];             % alpha and beta parameter of beta pdf over LAMBDA

% Units -- used just for plotting in axis labels and titles
psy.units.x = 'percent contrast';
psy.units.mu = 'percent contrast';
psy.units.sigma = 'percent contrast';
psy.units.lambda = [];

method = 'ent';     % Minimize the expected posterior entropy
% vars = [1 0 0];   % Minimize posterior entropy of the mean only
vars = [1 1 0];     % Minimize joint posterior entropy of mean, sigma and lambda
plotflag = 1;       % Plot visualization

[x,psy] = psybayes(psy, method, vars, [], []); % initial point recommendation

for itrial = 1:nTrials;
    
    % ========== MAKE STIM PATCH ==============
    x/100
    designMat(itrial,1) = x/100;
    im = makeGabor(stimLength, designMat(itrial,1));
    StimPatch= Screen('MakeTexture',windowPtr,im);
    
    % ========== FIXATION FOR BEGINNING OF TRIAL ==========
    Screen('fillRect',windowPtr,prefs.bgColor);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    t0 = GetSecs();
    % inittirladur = toc
    Screen('flip',windowPtr); % tic;
    while (GetSecs()-t0)<prefs.fix1Dur;
        % do nothing
    end
    
    % ========== PRESENT STIMULI FOR FIRST TIME ===========
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    pres1orientations = stimuliMat(itrial,1:setsize);
    for istim = 1:setsize;
        srcrect = [0 0 StimSizes];
        destrect = CenterRectOnPoint(srcrect,x_positions(istim),y_positions(istim));
        Screen('DrawTexture', windowPtr, StimPatch, srcrect,destrect, 180-pres1orientations(istim));
    end
    Screen('flip',windowPtr); % tic;
    t0 = GetSecs();
    while (GetSecs()-t0)<prefs.pres1Dur;
        % do nothing
    end
    
    % ========== ISI ===========
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    Screen('flip',windowPtr); % tic;
    t0 = GetSecs();
    while (GetSecs()-t0)<prefs.ISIdelayNum;
        % do nothing
    end
    
    % ========== ONE PROBE SHOWS UP. RESPONSE ===========
    
    stim = stimuliMat(itrial,end);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    if (prefs.lineAtPres2)
        pres2orientation = mod(90 - pres1orientations(stim) - designMat(itrial,2),180);
        if (pres2orientation == 0); pres2orientation = 180; end
        lineStim = lineCoordinates(:,pres2orientation);
        xy = [-lineStim lineStim ];
        Screen('DrawLines',windowPtr, xy, prefs.lineWidth,prefs.stimColor,[x_positions(stim) y_positions(stim)],1);
    else
        pres2orientation = pres1orientations(stim) + designMat(itrial,2);
        destrect = CenterRectOnPoint(srcrect,x_positions(stim),y_positions(stim));
        Screen('DrawTexture', windowPtr, StimPatch, srcrect,destrect, 180-pres2orientation);
    end
    
    % check response
    [pressedKey] = waitForKeys(prefs.keys,GetSecs());
    if pressedKey == length(prefs.keys);
        sca;
        ShowCursor;
        fclose('all');
        clear all;
    else
        %         designMat(itrial,5+setsize) = prefs.keysNum(pressedKey);
    end
    
    % check whether it is correct or not
    if (sign(designMat(itrial,2)/prefs.keysNum(pressedKey)) == 1) % if correct
        designMat(itrial,3) = 1;
    else
        designMat(itrial,3) = 0;
    end
    
    % ========== ITI ===========
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    Screen('flip',windowPtr); % tic;
    t0 = GetSecs();
    while (GetSecs()-t0)<prefs.ITIDur;
        % do nothingqs
    end
    
    % pick next contrast point and save file
    [x, psy] = psybayes(psy, method, vars, x, designMat(itrial,3));
    save(filename,'stimuliMat','designMat','lapsedesignMat','lapsestimuliMat','prefs','subjid','psy')
    
    if plotflag
        psybayes_plot(psy);
        drawnow;
    end
    
end

% ======================================
% save file and exit
% ======================================

save(filename,'stimuliMat','designMat','lapsedesignMat','lapsestimuliMat','prefs','subjid','psy')

textx = screenCenter(1) - 200;
texty = screenCenter(2) - 200;
Screen('TextSize',windowPtr,20);
Screen('TextFont',windowPtr,'Helvetica');
Screen('DrawText',windowPtr,'This part of the experiment is over',textx,texty,[255 255 255]);
Screen('DrawText',windowPtr,'Please get the experimenter',textx,texty+100,[255 255 255]);
Screen('Flip', windowPtr);
pause;

% % exit
% sca;
% ShowCursor;
% fclose('all');
% clear all;

% =========== HELPER FUNCTIONS =============

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
        %     Screen('DrawLines',windowPtr,[x-sizee,y; x+sizee,y],2,color);
        %     Screen('DrawLines',windowPtr,[x,y-sizee; x,y+sizee],2,color);
        
        %     Screen('DrawLines',windowPtr, xy, prefs.lineWidth,prefs.stimColor,[xpositions(j) ypositions(j)],1)
    end
end