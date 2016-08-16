function MOC_getcontrast(subjid)

mixedcontrast = 0;  % if you want the 2 full contrast or all same contrast

% open screen
screenNumber = max(Screen('Screens'));       % use external screen if exists
[w, h] = Screen('WindowSize', screenNumber);  % screen resolution of smaller display
screenCenter = [w h]/2;    
windowPtr = Screen('OpenWindow',screenNumber,128*ones(1,3),[],32,2);
Screen(windowPtr,'BlendFunction',GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
Screen('TextSize',windowPtr,20);
Screen('TextFont',windowPtr,'Helvetica');
HideCursor;
commandwindow;

% =========== SETTING EXPERIMENTAL STUFF ===========
prefs = prefscode('Delay',subjid,2,1);
setsize = 4;
filename = ['output_mat/ContrastFn_MOC_' upper(subjid) '_'  datestr(now, 30) ...
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
nTrialsperContrast = 20;
if (mixedcontrast); nTrialsperContrast = 40; end
contrastVec = exp(linspace(log(.05),log(1),10)); % 10 contrasts linearly spaced in log space
nTrials = 10*nTrialsperContrast;

designMat = [repmat(contrastVec',nTrialsperContrast,1) 30*ones(nTrials,1) nan(nTrials,1)];
designMat(rand(nTrials,1) > 0.5,2) = -30;
designMat(randperm(nTrials),:) = designMat;  % permuting trials
stimuliMat = [round(rand(nTrials,setsize)*180) randi(setsize,nTrials,1)];

names.designMat = {'contrast','delta','correct'};
names.stimuliMat = {'stim1pres1','stim2pres1','stim3pres1','stim4pres1','target_position'};

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

% make gabor patch with full contrast
reliability = 1; % contrast of 1
im = makeGabor(stimLength, reliability);
StimPatch1= Screen('MakeTexture',windowPtr,im);

% ===========================================
% STARTING SCREEN
% ===========================================

textx = screenCenter(1) - 200;
texty = screenCenter(2) - 200;
Screen('TextSize',windowPtr,20);
Screen('TextFont',windowPtr,'Helvetica');
Screen('DrawText',windowPtr,'Push any key to begin',textx,texty,[255 255 255]);
Screen('Flip', windowPtr);
pause;

% =====================================
% ADAPTIVE METHOD
% =====================================

for itrial = 1:nTrials;
    
    % ========== MAKE STIM PATCH ==============
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
    k = randperm(setsize);
    for istim = 1:setsize/2;
        srcrect = [0 0 StimSizes];
        destrect = CenterRectOnPoint(srcrect,x_positions(k(istim)),y_positions(k(istim)));
        Screen('DrawTexture', windowPtr, StimPatch, srcrect,destrect, 180-pres1orientations(k(istim)));
    end
    for istim = (setsize/2)+1:setsize; % making half of them full contrast
        srcrect = [0 0 StimSizes];
        destrect = CenterRectOnPoint(srcrect,x_positions(k(istim)),y_positions(k(istim)));
        if (mixedcontrast)
            Screen('DrawTexture', windowPtr, StimPatch1, srcrect,destrect, 180-pres1orientations(k(istim)));
        else
            Screen('DrawTexture', windowPtr, StimPatch, srcrect,destrect, 180-pres1orientations(k(istim)));
        end
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
        pres2orientation = mod(90 - pres1orientations(k(stim)) - designMat(itrial,2),180);
        if (pres2orientation == 0); pres2orientation = 180; end
        lineStim = lineCoordinates(:,pres2orientation);
        xy = [-lineStim lineStim ];
        Screen('DrawLines',windowPtr, xy, prefs.lineWidth,prefs.stimColor,[x_positions(k(stim)) y_positions(k(stim))],1);
    else
        pres2orientation = pres1orientations(k(stim)) + designMat(itrial,2);
        destrect = CenterRectOnPoint(srcrect,x_positions(k(stim)),y_positions(k(stim)));
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
    
    % save file
    save(filename,'stimuliMat','designMat','prefs','subjid')

    
end

% ======================================
% save file and exit
% ======================================

save(filename,'stimuliMat','designMat','prefs','subjid')

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