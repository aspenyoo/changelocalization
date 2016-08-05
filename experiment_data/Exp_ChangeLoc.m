function Exp_ChangeLoc(subjid, exptype, sessionnum, nTrialsPerCond, nSessions, ispractice)
% runs 2AFC experiment: orientation change detection task.
%
% STRUCT EXPLANATIONS
% prefs: experimental preferences
% D: related to data (things you may analyze about the experiment in the
% future)

% ========================================================================
% SETTING THINGS UP (changing this section is usually unnecessary)
% ========================================================================
if nargin < 3; sessionnum = []; end
if nargin < 4; nTrialsPerCond = []; end
if nargin < 5; nSessions = 8; end
if nargin < 6; ispractice = 0; end

commandwindow;

% random number generator
rng('shuffle');

% ==================================================================
% PREFERENCES (should be changed to match experiment preferences)
% ==================================================================

screenDistance = 40;                      % distance between observer and screen (in cm)

% importing preferences for simultaneous/sequential experiment
if (ispractice)
    prefs = prefscode(['Pract_' exptype],subjid,sessionnum,nTrialsPerCond);
else
    prefs = prefscode(exptype,subjid,sessionnum,nTrialsPerCond);
end

% response keys
% prefs.keys = [KbName('9') KbName('7') KbName('1') KbName('3') KbName('esc')];
% prefs.keysNum = [9 7 1 3 0];
prefs.keys = [KbName('9') KbName('7') KbName('1') KbName('3') KbName('esc')];
prefs.keysNum = [9 7 1 3 0];

prefs.keys2 = [KbName('a') KbName('s') KbName('d') KbName('f')];
prefs.keysNum2 = [1 2 3 4];

% Data file
prefs.fidmat = fullfile('output_mat',[prefs.fileName '.mat']);

% calculating full experimental design and pseudo-randomizing order
prefs.design = fullfact([prefs.f1 prefs.f2 prefs.f3 prefs.f4]);
prefs.nCond = size(prefs.design,1);

% ========================================================================
% CALCULATIONS BASED ON PREFERENCES (change not necessary)
% ========================================================================
% skipping sync tests
% Screen('Preference', 'SkipSyncTests', 1);

% screen info (visual)
screenNumber =  max(Screen('Screens'));       % use external screen if exists
[w, h] = Screen('WindowSize', screenNumber);  % screen resolution of display
screenResolution = [w h];                 % screen resolution
screenCenter = screenResolution/2;       % screen center
screenAngle = atand((prefs.screenHeight/2) / screenDistance) ; % total visual angle of screen (height)
screen_ppd = h / screenAngle;  % pixels per degree
prefs.stimArea = prefs.stimArea * screen_ppd^2; % stimulus area (pixels)
stimLength = round(sqrt(prefs.stimArea));

% open screen
if isempty(sessionnum)
    windowPtr = Screen('OpenWindow',screenNumber,prefs.grey,[],32,2);
else
    windowPtr = 10;
end
Screen(windowPtr,'BlendFunction',GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
HideCursor;

% ========================================================================
% EXPERIMENT TRIAL AND STIMULI INFORMATION/CALCULATIONS
% ========================================================================
setsize = max(prefs.pres1stimNum);


% CREATE/LOAD TRIAL VARIABLES
% -------------------------------------------------------------------------

% great designMat for all sessions
if (sessionnum == 1)
    
    % DESIGN MAT: conditions and responses of experiment
    designMat = [];
    designMat(:,1) = prefs.deltaNum(prefs.design(:,1));
    designMat(:,2) = prefs.pres1stimNum(prefs.design(:,2));
    designMat(:,3:2+setsize) = prefs.reliabilityNum(prefs.design(:,3),:);
    designMat(:,3+setsize) = prefs.ISIdelayNum(prefs.design(:,4));
    designMat(:,4+setsize) = 1:prefs.nCond;
    if (prefs.vmprior)
        designMat = repmat(designMat, [prefs.nTrialsPerCond 1]);
        designMat(:,1) = designMat(:,1).*circ_vmrnd(0,prefs.vmprior, [size(designMat,1),1])*180/pi; %sort(repmat([-45:10:45]',[prefs.nCond 1]));
    else  % method of constant stimuli
        if (ispractice)
            constantStim = ceil(rand(3,1)*180-90); % 3 levels
            constantStim(constantStim == 0) = 10;
        else
            constantStim = [-87.5:5:87.5]';
            % constantStim = [-90:5:-5 5:5:85]'; % 35 levels
        end
        designMat = sort(repmat(designMat,[length(constantStim),1]));
        designMat(:,1) = designMat(:,1).*repmat(constantStim,[prefs.nCond 1]);
        designMat = repmat(designMat, [prefs.nTrialsPerCond 1]);
    end
    designMat = [designMat nan(size(designMat,1),3)];
    prefs.nTrials = size(designMat,1);
    prefs.ncurrTrials = round(prefs.nTrials/nSessions);
    designMat = designMat(randperm(prefs.nTrials),:);
    designMat = [designMat nan(prefs.nTrials,2)];
    nPC = floor(prefs.nTrials/prefs.feedbacktrial);
    
    % names of variables in designMat
    names.designMat = cell(1,8);
    names.designMat{1} = 'delta'; names.designMat{2} = 'set size';
    names.designMat{3+setsize} = 'delay time';
    names.designMat{4+setsize} = 'condition number'; names.designMat{5+setsize} = 'response';
    names.designMat{6+setsize} = 'Correct?'; names.designMat{7+setsize} = 'RT';
    
    switch prefs.stimType
        case 'ellipse'
            names.designMat{3} = 'ellipse reliability';
        case 'gabor'
            names.designMat{3} = 'gabor contrast';
    end
    
    % calculating where breaks will occur
    if prefs.blocknum >1;
        prefs.breakpoints = round((1:(prefs.blocknum-1)).*(prefs.nTrials/prefs.blocknum));
    else
        prefs.breakpoints = prefs.nTrials+1;
    end
    
    % STIMULI MAT: screen positions, orientations, target locations
    stimuliMat = nan(prefs.nTrials, 2*setsize + 1);
    % names of columns
    names.stimuliMat = cell(1,(2*setsize+1));
    for j = 1:setsize;
        switch exptype
            case 'Delay'
                names.stimuliMat{j} = ['stimulus ' num2str(j) ' presentation interval'];
            case 'Contrast'
                names.stimuliMat{j} = ['stimulus' num2str(j) 'contrast'];
        end
        names.stimuliMat{j+setsize} = ['stimulus' num2str(j) ' orientation (pres 1)'];
    end
    names.stimuliMat{2*setsize+1} = 'target location number';
    
    stimuliMat(:,end) = 1-(designMat(:,1)==0);
    stimuliMat(:,end) = ceil(rand(prefs.nTrials,1).*stimuliMat(:,end).*setsize);
    
    % initial stimulus orientations
    stimuliMat(:,setsize+1:2*setsize) = round(180*rand(prefs.nTrials,setsize));
    
    % stimulus location
    randloc = 0;
    if (randloc)
        %     for itrial = 1:prefs.nTrials
        %         % Pick angle of first one rand and add the rest in counterclockwise
        %         % fashion with angle spacing = 2pi/max(E.setSizeNum)
        %         locAngles = rand*2*pi+(1:setSize)*(2*pi)/setSize;
        %
        %         [X, Y] = pol2cart(locAngles, screen_ppd * prefs.stimecc);
        %         % positions with jitter
        %         % ASPEN: OUTDATEDD
        % %         stimuliMat(itrial,1:setSize) = X(1:setSize) + screenCenter(1)...
        % %             + round((rand(1,setSize)-.5)*prefs.jitter*screen_ppd);
        % %         stimuliMat(itrial,setSize+1:setSize+setSize) = Y(1:setSize) + screenCenter(2)...
        % %             + round((rand(1,setSize)-.5)*prefs.jitter*screen_ppd);
        %     end
    else
        locAngles = pi/4-(1:setsize)*(2*pi)/setsize;
        [X, Y] = pol2cart(locAngles, screen_ppd * prefs.stimecc);
        x_positions = X(1:setsize) + screenCenter(1)...
            + round((rand(1,setsize)-.5)*prefs.jitter*screen_ppd);
        y_positions = Y(1:setsize) + screenCenter(2)...
            + round((rand(1,setsize)-.5)*prefs.jitter*screen_ppd);
    end
    trial = 1; % starting on the first trial
else
    load(prefs.fidmat); % load previous session
    trial = find(isnan(designMat(:,end)),1,'first'); % figure out what trial you are on
end

% matrix of all possible stimuli
% -----------------------------------------------------------------------
% (res degrees between each possible stimulus orientations)

clear StimPatches;

% Set the number of possible orientations based on resolution

% ellipse stuff?
res = 1; % resolution
numDegrees = ceil(180/res);

uniqueReliabilities = unique(prefs.reliabilityNum);
switch prefs.stimType
    case 'ellipse'
        StimPatches = zeros(length(uniqueReliabilities),numDegrees); % holds the stimulus image patches
        StimSizes = zeros(length(uniqueReliabilities),numDegrees+1,2); % holds stimulus image sizes
        
        % Fill StimPatches and StimSizes
        for irel = 1:length(prefs.reliabilityNum)
            % Eccentricity = reliability for now
            reliability = uniqueReliabilities(irel);
            % Draw a patch for each orientation
            for ideg = 1:numDegrees
                im = drawEllipse(prefs.stimArea,reliability,ideg*res,prefs.stimColor,prefs.bgColor);
                StimPatches(irel,ideg) = Screen('MakeTexture',windowPtr,im);
                StimSizes(irel,ideg,:) = size(im);
                StimSizes(irel,ideg,:) = StimSizes(irel,ideg,[2 1]);
            end
        end
    case 'gabor'
        StimPatches = zeros(1,length(uniqueReliabilities)); % make a stimpatch for each unique reliability value
        StimSizes = [stimLength stimLength];
        for irel = 1:length(uniqueReliabilities)
            reliability = uniqueReliabilities(irel);
            im = makeGabor(stimLength, reliability);
            StimPatches(irel) = Screen('MakeTexture',windowPtr,im);
            %             StimPatches(i) = CreateProceduralGabor(windowPtr, stimLength, stimLength ,[],[0.5 0.5 0.5 0.0], [], reliability);
        end
end

if (prefs.lineAtPres2)
    % ================ LINE SIMULI ===================
    res = 1; % resolution
    numDegrees = ceil(180/res);
    % prefs.lineArea = prefs.stimArea; % so line is same area as ellipse
    prefs.lineWidth = 3; % pixels
    prefs.lineLength = stimLength;
    lineCoordinates = nan(2,numDegrees);
    for j = 1:numDegrees;
        [lineCoordinates(1,j), lineCoordinates(2,j)] = lineCoord(prefs.lineLength,j);
    end
end
% ========================================================================
% RUN THE EXPERIMENT
% ========================================================================
dy = 30;
textx = w/2-300;
expTime = GetSecs();
maxexpTime = 2*60; % 30 minutes (in seconds)

if (sessionnum)
    % ====== DETECTION INFO SCREEN =======
    texty = screenCenter(2) - 200;
    Screen('TextSize',windowPtr,20);
    Screen('TextFont',windowPtr,'Helvetica');
    Screen('DrawText',windowPtr,'Indicate which orientation changed.',textx,texty,[255 255 255]);
    Screen('Flip', windowPtr);
    waitForKey;
end


% run a trial
% -------------------------------------------------------------------------
fbtrial = 0; progPC = [];
for itrial = trial:prefs.ncurrTrials*nSessions;
    
    % setting values for current trial
    %     condition = designMat(itrial,5); % current condition
    pres1orientations = stimuliMat(itrial,setsize+1:2*setsize);
    pres2orientations = pres1orientations;
    locChange = stimuliMat(itrial,end);
    pres2orientations(locChange) = pres1orientations(locChange) + designMat(itrial,1);
    ISI2delay = designMat(itrial,3+setsize); % delay time between intervals that transition between first and second presentations...(poorly explained)
    
    % adjusting number to be between 1-180 for stimulus presentations
    pres1orientations = mod(round(pres1orientations),180);
    pres1orientations(pres1orientations == 0) = 180;
    pres2orientations = mod(round(pres2orientations),180);
    pres2orientations(pres2orientations == 0) = 180;
    if (prefs.lineAtPres2)
        pres2orientations = mod(90 - pres2orientations,180);
        if (pres2orientations == 0); pres2orientations = 180; end
    end
    
    %     lineStim = lineCoordinates(:,pres2orientations);
    
    % trial initiation screen
    Screen('fillRect',windowPtr,prefs.bgColor);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength.*1.5); % make slightly longer fixation cross
    t0 = GetSecs();
    Screen('flip',windowPtr); % tic;
    while (GetSecs()-t0)<prefs.initTrialDur;
        % do nothing
    end
    
    % blank fixation screen
    Screen('fillRect',windowPtr,prefs.bgColor);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    t0 = GetSecs();
    % inittirladur = toc
    Screen('flip',windowPtr); % tic;
    while (GetSecs()-t0)<prefs.fix1Dur;
        % do nothing
    end
    
    % first stimulus presentation
    % =========================================================
    
    stimPerInterval = [designMat(itrial,2) setsize];
    
    % if nothing in first interval...
    if ~(prefs.stimPresent) && (stimPerInterval(1) == 0) % sequential presentation & nothing in first interval
        
        % blank screen for stimulus duration + ISI
        Screen('fillRect',windowPtr,prefs.bgColor);
        drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
        % fix1 = toc
        Screen('flip',windowPtr); % tic;
        t0 = GetSecs();
        while (GetSecs()-t0) < (prefs.pres1Dur + prefs.pres1ISIDur);
            % do nothing
        end
    end
    
    % the stimulus presentations
    Screen('fillRect',windowPtr,prefs.bgColor);
    if (prefs.permLocInPres1) % if you want to permute the locations in presentation 1
        k = randperm(setsize);
    else
        k = 1:setsize;
    end
    stimuliMat(itrial,k) = designMat(itrial,3:2+setsize);  % position index. from top right counterclockwise. only for contrast exptype
    
    for istim = 1:setsize;
        switch prefs.stimType
            case 'ellipse'
                srcrect = [0 0 squeeze(StimSizes(uniqueReliabilities == designMat(itrial,2+istim),pres1orientations(k(istim)),:))'];
                destrect = CenterRectOnPoint(srcrect,x_positions(k(istim)),y_positions(k(istim)));
                Screen('drawtexture',windowPtr,StimPatches(uniqueReliabilities == designMat(itrial,istim),pres1orientations(k(istim))),srcrect,destrect,0);
            case 'gabor'
                srcrect = [0 0 StimSizes];
                destrect = CenterRectOnPoint(srcrect,x_positions(k(istim)),y_positions(k(istim)));
                Screen('DrawTexture', windowPtr, StimPatches(uniqueReliabilities == designMat(itrial,2+istim)), srcrect,destrect, 180-pres1orientations(k(istim)));
                if strcmp(exptype,'Delay')
                    stimuliMat(itrial,k(istim)) = 2-(istim <= stimPerInterval(1));% index of which interval
                end
        end
        if ~prefs.stimPresent % sequential presentation
            if sum(istim == stimPerInterval); % flip window when the proper number of stimuli are on the screen.
                if (prefs.stimecc)
                    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
                end
                % blank = toc
                Screen('flip',windowPtr);
                % tic;
                t0 = GetSecs();
                while (GetSecs()-t0) < prefs.pres1Dur;
                    % do nothing
                end
                
                %blank screen in between
                Screen('fillRect',windowPtr,prefs.bgColor);
                drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
                % stim1 = toc
                Screen('flip',windowPtr);
                % tic;
                t0 = GetSecs();
                
                if (istim == setsize); % if it's last interval for first presentation
                    while (GetSecs()-t0) < ISI2delay; % ISI between first and second presentation intervals
                        % do nothing
                    end
                else % if there is going to be another interval for first presentation
                    while (GetSecs()-t0) < prefs.pres1ISIDur; % ISI between first presentation intervals
                        % do nothing
                    end
                end
            end
        end
    end
    
    if prefs.stimPresent % if simultaneous presentation
        if (prefs.stimecc)
            drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
        end
        Screen('flip',windowPtr);
        % tic
        t0 = GetSecs();
        while (GetSecs()-t0) < prefs.pres1Dur;
            % do nothing
        end
        
        
        % ISI
        Screen('fillRect',windowPtr,prefs.bgColor);
        drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
        % stim1 = toc
        Screen('flip',windowPtr);
        % tic;
        t0 = GetSecs();
        while (GetSecs()-t0) < ISI2delay; % ISI between first and second presentation intervals
            % do nothing
        end
    end
    
    % if need a blank screen for a pres1 period of time at the end
    if ~(prefs.stimPresent) && (stimPerInterval(1) == setsize) % simultaneous presentation
        % blank screen
        Screen('fillRect',windowPtr,prefs.bgColor);
        drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
        % blank = toc
        Screen('flip',windowPtr); % tic;
        t0 = GetSecs();
        while (GetSecs()-t0) < (prefs.pres1Dur + prefs.pres1ISIDur);
            % do nothing
        end
    end
    
    % second stimulus presentation
    Screen('fillRect',windowPtr,prefs.bgColor);
    if (prefs.stimecc)
        drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    end
    srcrect = [0 0 StimSizes];
    if (prefs.allStimInPres2) % if presenting all stimuli at once
        for istim= 1:setsize % for each stimulus, draw appropriate stimulus
            if (prefs.lineAtPres2)
                lineStim = lineCoordinates(:,pres2orientations(k(istim)));
                xy = [-lineStim lineStim ];
                Screen('DrawLines',windowPtr, xy, prefs.lineWidth,prefs.stimColor,[x_positions(k(istim)) y_positions(k(istim))],1);
            else
                destrect = CenterRectOnPoint(srcrect,x_positions(k(istim)),y_positions(k(istim)));
                Screen('DrawTexture', windowPtr, StimPatches(uniqueReliabilities == max(uniqueReliabilities)), srcrect,destrect, 180-pres2orientations(k(istim)));
            end
        end
    else
        Screen('DrawTexture', windowPtr, StimPatches(uniqueReliabilities == max(uniqueReliabilities)), srcrect,destrect, 180-pres2orientations(k(istim)));
    end
    % ISI = toc
    Screen('flip',windowPtr);
    % fprintf('ISI: %f \n', toc)
    t0 = GetSecs();
    
    if (prefs.respInPres2) % if response is in 2nd stim presentation
        % leave pres2 stimuli on screen until response
    else % if response is not during 2nd presentation
        % leave pres2 stimuli on screen for fixed duration
        while (GetSecs()-t0)<prefs.pres2Dur;
            % do nothing
        end
        % blank screen (waiting for response)
        Screen('fillRect',windowPtr,prefs.bgColor);
        drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
        Screen('flip',windowPtr);
    end
    
    % check response
    [pressedKey, designMat(itrial,7+setsize)] = waitForKeys(prefs.keys,GetSecs());
    if pressedKey == length(prefs.keys);
        sca;
        ShowCursor;
        fclose('all');
        clear all;
    else
        designMat(itrial,5+setsize) = prefs.keysNum(pressedKey);
    end
    
    % calculating correct resp
    designMat(itrial,6+setsize) = pressedKey == stimuliMat(itrial,end);
    
    % confidence judgment
    textx = screenCenter(1) - 3*screen_ppd;
    texty = screenCenter(2) - screen_ppd;
    Screen('TextSize',windowPtr,24);
    Screen('TextFont',windowPtr,'Helvetica');
    Screen('DrawText',windowPtr,'Confidence?',textx,texty,[255 255 255]);
    Screen('Flip', windowPtr);
    [designMat(itrial,8+setsize), designMat(itrial,9+setsize)] = waitForKeys(prefs.keys2,GetSecs());
    
    % ITI
    if (prefs.feedback) % if sound feedback
        if (designMat(itrial,6+setsize))
            beep = MakeBeep(1200,.2);
        else
            beep = MakeBeep(500,.2);
        end
        Snd('Open');
        Snd('Play',beep);
    end
    
    Screen('fillRect',windowPtr,prefs.bgColor);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
    Screen('flip',windowPtr); tic
    t0 = GetSecs();
    while (GetSecs()-t0)<prefs.ITIDur;
        % do nothing
    end
    % ITI = toc
    
    save(prefs.fidmat)
    
    % ======= BREAK/BLOCK SECTION =======
    % if big block break
    if ~isempty(intersect(prefs.breakpoints,itrial))
        Screen('fillRect',windowPtr,prefs.bgColor);
        %tic
        
        while toc<prefs.breakDuration
            Screen('TextSize',windowPtr,20);
            Screen('DrawText',windowPtr,['You have finished ' num2str(find(prefs.breakpoints == itrial)) '/' num2str(length(prefs.breakpoints)+1) ' of the trials.'],100,screenCenter(2)-80,[255 255 255]);
            Screen('DrawText',windowPtr,['Please take a short break now. You can continue in ' num2str(round(prefs.breakDuration-toc)) ' seconds.'],100,screenCenter(2)-120,[255 255 255]);
            Screen('Flip', windowPtr);
        end
    end
    
    % if little feedback block break
    if ~mod(itrial,prefs.feedbacktrial) % every FEEDBACKTRIALth trial
        
        progPC = [progPC round(100*mean(designMat(fbtrial+1:itrial,6+setsize)))];
        fbtrial = itrial;
        %             progPC = mean(designMat(1:i,7));
        Screen('TextSize',windowPtr,20);
        Screen('fillRect',windowPtr,prefs.bgColor);
        Screen('DrawText',windowPtr,['You have earned ' num2str(progPC(end)) ' points. Press any key to continue'],250,screenCenter(2)-80,[255 255 255]);
        
        % GRAPH OF PROGRESS
        graphwidth = round(w/3); % half of graph width
        graphheight = round(h/5); % half of graph height
        origin = [w/2-graphwidth/2 h/2+graphheight];
        Screen('DrawLine',windowPtr,255*ones(1,3),origin(1), origin(2),origin(1)+graphwidth,origin(2)); % x-axis
        Screen('DrawLine',windowPtr,255*ones(1,3),origin(1), origin(2),origin(1),origin(2)-graphheight); % y-axis
        Screen('TextSize',windowPtr,12);
        Screen('DrawText',windowPtr,'points per block',w/2-50,origin(2)+15,[255 255 255]); % xlabel
        
        % draw lines to make graph
        pc = (progPC./75 - 1/3)*graphheight; % rescaling PC to size of graph
        dx = graphwidth./nPC; % equally spaced out to fill size of graph at end of exp.
        og = origin;
        dc = nan(2,length(pc));
        dc(:,1) = [og(1),og(2)-pc(1)];
        for ipc = 2:length(pc); % draw the lines
            dc(:,ipc) = [og(1)+dx og(2)-pc(ipc)];
            Screen('DrawLine',windowPtr,255*ones(1,3),dc(1,ipc-1),dc(2,ipc-1),dc(1,ipc),dc(2,ipc)); % y-axis
            og(1) = og(1) + dx;
        end
        Screen('DrawDots',windowPtr,dc,5,255*ones(1,3),[],1)
        
        Screen('Flip', windowPtr);
        waitForKey;
        
        drawfixation(windowPtr,screenCenter(1),screenCenter(2),prefs.fixColor,prefs.fixLength);
        Screen('Flip', windowPtr);
        WaitSecs(prefs.fix1Dur);
    end
    
    if (GetSecs() - expTime) > maxexpTime;
        save(prefs.fidmat)
        break;
    end
    
end

% final screen
textx = screenCenter(1) - 3*screen_ppd;
texty = screenCenter(2) - 200;
Screen('TextSize',windowPtr,20);
Screen('TextFont',windowPtr,'Helvetica');
Screen('DrawText',windowPtr,['End of this task. You earned ' num2str(round(100*nanmean(designMat(:,6+setsize)))) ' points!'],textx,texty,[255 255 255]); texty = texty + 5*dy;
Screen('DrawText',windowPtr,'Press any key to end.',textx,texty,[255 255 255]); texty = texty + 5*dy;
Screen('Flip', windowPtr);
waitForKey;

Screen('Flip', windowPtr);
save(prefs.fidmat)

% ========================================================================
% ---------------------------HELPER FUNCTIONS-----------------------------
% ========================================================================

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

%----------------------------------------------------
    function keyCode = waitForKey
        keyCode = ones(1,256);
        while sum(keyCode(1:254))>0
            [keyIsDown,secs,keyCode] = KbCheck;
        end
        while sum(keyCode(1:254))==0
            [keyIsDown,secs,keyCode] = KbCheck;
        end
        keyCode = min(find(keyCode==1));
    end

%----------------------------------------------------
    function drawfixation(windowPtr,x,y,color,sizee)
        Screen('DrawLine',windowPtr,color,x-sizee,y,x+sizee,y,2);
        Screen('DrawLine',windowPtr,color,x,y-sizee,x,y+sizee,2);
        %     Screen('DrawLines',windowPtr,[x-sizee,y; x+sizee,y],2,color);
        %     Screen('DrawLines',windowPtr,[x,y-sizee; x,y+sizee],2,color);
        
        %     Screen('DrawLines',windowPtr, xy, prefs.lineWidth,prefs.stimColor,[xpositions(j) ypositions(j)],1)
    end

end
% ---------------------------------------------------
% function writeHeaderToFile(D, fid)
%
% h = fieldnames(D);
%
% for i=1:length(h)
%     fprintf(fid, '%s\t', h {i });
% end
% fprintf(fid, '\n');
%
%
% % ------------------------------------------------------
% function writeTrialToFile(D, trial, fid)
%
% h = fieldnames(D);
% for i=1:length(h)
%     data = D.(h {i })(trial);
%     if isnumeric(data)
%         fprintf(fid, '%s\t', num2str(data));
%     elseif iscell(data)
%         fprintf(fid, '%s\t', char(data));
%     else
%         error('wrong format!')
%     end
% end
% fprintf(fid, '\n');
