function train_loc_delay(subjid, sessionnum, nTrialsPerCond)
if nargin < 2; nTrialsPerCond = 5; end

% training for orientation discrimination task
% prefs = prefscode('Detect','Seq',subjid);
prefs = prefscode('Delay',subjid,sessionnum);

% Screen('Preference', 'SkipSyncTests', 1);

% screen info (visual)
screenNumber = max(Screen('Screens'));       % use external screen if exists
[w, h] = Screen('WindowSize', screenNumber);  % screen resolution of smaller display
screenResolution = [w h];       % screen resolution
screenCenter = screenResolution/2;       % screen center
screenDistance = 40;                      % distance between observer and screen (in cm)
screenAngle = atand((prefs.screenHeight/2) / screenDistance) ; % total visual angle of screen
screen_ppd = h / screenAngle;  % pixels per degree
prefs.stimArea = prefs.stimArea * screen_ppd^2; % ellipse area (pixels)
stimLength = round(sqrt(prefs.stimArea));

% prefs.lineArea = prefs.stimArea; % so line is same area as ellipse
prefs.lineWidth = 3;%round(prefs.lineArea*.01); % pixels
prefs.lineLength = stimLength;

% open screen
windowPtr = 10;
% windowPtr = Screen('OpenWindow',screenNumber,prefs.grey,[],32,2);
Screen(windowPtr,'BlendFunction',GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

Screen('TextSize',windowPtr,20);
Screen('TextFont',windowPtr,'Helvetica');
textx = 4*screen_ppd;
texty = screenCenter(2) - 2*screen_ppd;
dy = 30;

if sessionnum == 1;
    
    % ========== EXPLANATION OF TASK ==========
    texty = screenCenter(2) - 150;
    Screen('DrawText',windowPtr,'Four oriented blobs will flash up in some order, and all together after',textx,texty,[255 255 255]);
    Screen('Flip', windowPtr);
    pause;
    
    
    %  ========== EXAMPLE TRIAL ============
    % stuff that are reused.
    reliability = prefs.reliabilityNum(round(length(prefs.reliabilityNum)/2));
    im = makeGabor(stimLength, reliability);
    StimPatch = Screen('MakeTexture',windowPtr,im);
    StimSize = size(im);
    yshift = 50;
    
    % FIXATION
    Screen('DrawText',windowPtr,'Four oriented blobs will flash up in some order, and all together after',textx,texty,[255 255 255]);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2)+yshift,prefs.fixColor,prefs.fixLength);
    Screen('Flip', windowPtr);
    pause(0.5);
    
    % STIMULUS PRESENTATION 1
    setSize = 4; % max(prefs.pres1stimNum);
    locAngles = pi/4+2*pi+(1:setSize)*(2*pi)/setSize;
    stimecc = 2;
    [X, Y] = pol2cart(locAngles, screen_ppd * stimecc);
    X = X+screenCenter(1);
    Y = Y+screenCenter(2)+yshift;
    
    % 1
    Screen('DrawText',windowPtr,'Four oriented blobs will flash up in some order, and all together after',textx,texty,[255 255 255]);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2)+yshift,prefs.fixColor,prefs.fixLength);
    deg1 = round(rand(1,setSize).*180);
    for istim = 1:ceil(setSize/2);
        cuesrcrect = [0 0 StimSize];
        destrect = CenterRectOnPoint(cuesrcrect,X(istim),Y(istim));
        Screen('DrawTexture',windowPtr,StimPatch,cuesrcrect,destrect,180-deg1(istim));
    end
    Screen('Flip', windowPtr);
    pause(0.3);
    
    % stim1 ISI
    Screen('DrawText',windowPtr,'Four oriented blobs will flash up in some order, and all together after',textx,texty,[255 255 255]);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2)+yshift,prefs.fixColor,prefs.fixLength);
    Screen('Flip', windowPtr);
    pause(1);
    
    % 2
    Screen('DrawText',windowPtr,'Four oriented blobs will flash up in some order, and all together after',textx,texty,[255 255 255]);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2)+yshift,prefs.fixColor,prefs.fixLength);
    for istim = ceil(setSize/2)+1:setSize;
        cuesrcrect = [0 0 StimSize];
        destrect = CenterRectOnPoint(cuesrcrect,X(istim),Y(istim));
        Screen('DrawTexture',windowPtr,StimPatch,cuesrcrect,destrect,180-deg1(istim));
    end
    Screen('Flip', windowPtr);
    pause(0.3);
    
    
    % ISI
    Screen('DrawText',windowPtr,'Four oriented blobs will flash up in some order, and all together after',textx,texty,[255 255 255]);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2)+yshift,prefs.fixColor,prefs.fixLength);
    Screen('Flip', windowPtr);
    pause(1);
    
    % STIMULUS 2
    Screen('DrawText',windowPtr,'Four oriented blobs will flash up in some order, and all together after',textx,texty,[255 255 255]);
    drawfixation(windowPtr,screenCenter(1),screenCenter(2)+yshift,prefs.fixColor,prefs.fixLength);
    for istim = 1:setSize;
        cuesrcrect = [0 0 StimSize];
        destrect = CenterRectOnPoint(cuesrcrect,X(istim),Y(istim));
        Screen('DrawTexture',windowPtr,StimPatch,cuesrcrect,destrect,180-deg1(istim));
    end
    Screen('Flip', windowPtr);
    pause;
    
    
    % ========== UNIFORM PRIOR OVER CHANGE ==========
    texty = screenCenter(2) - 150;
    Screen('DrawText',windowPtr,'The orientation of one of the blobs will change in EVERY TRIAL.',textx,texty,[255 255 255]); texty = texty + 3*dy;
    Screen('DrawText',windowPtr,'The change will be random.',textx,texty,[255 255 255]); texty = texty + 5*dy;
    
    % EXAMPLES
    Screen('DrawText',windowPtr,'Examples:',textx,texty,[255 255 255]);
    Screen('Flip', windowPtr);
    pause;
    
    nExamples = 3; nSamples = 10;
    newx = textx;
    for iexamp = 1:nExamples;
        texty = screenCenter(2) - 150;
        Screen('DrawText',windowPtr,'The orientation of one of the blobs will change in EVERY TRIAL.',textx,texty,[255 255 255]); texty = texty + 3*dy;
        Screen('DrawText',windowPtr,'The change will be random.',textx,texty,[255 255 255]); texty = texty + 5*dy;
        
        % EXAMPLES
        Screen('DrawText',windowPtr,'Examples:',textx,texty,[255 255 255]);
        
        % draw one gabor
        deg1 = round(rand*180);
        reliability = prefs.reliabilityNum(randi(length(prefs.reliabilityNum)));
        im = makeGabor(stimLength, reliability);
        StimPatch = Screen('MakeTexture',windowPtr,im);
        StimSize = size(im);
        cuesrcrect = [0 0 fliplr(StimSize)];
        newx = newx + 200; newy = texty + 150;
        destrect = CenterRectOnPoint(cuesrcrect,newx,newy);
        Screen('DrawTexture',windowPtr,StimPatch,cuesrcrect,destrect,180-deg1);
        Screen('Flip', windowPtr);
        pause;
        
        for isamp = 1:nSamples;
            
            % redrawing the same stuff that was already there
            texty = screenCenter(2) - 150;
            Screen('DrawText',windowPtr,'The orientation of one of the blobs will change in EVERY TRIAL.',textx,texty,[255 255 255]); texty = texty + 3*dy;
            Screen('DrawText',windowPtr,'The change will be random.',textx,texty,[255 255 255]); texty = texty + 5*dy;
            Screen('DrawText',windowPtr,'Examples:',textx,texty,[255 255 255]);
            Screen('DrawTexture',windowPtr,StimPatch,cuesrcrect,destrect,180-deg1);
            
            % draw a bunch of lines from that one ellipse
            deg2 = deg1+round(rand*180-90);
            [x,y] = lineCoord(prefs.lineLength, deg2);
            Screen('DrawLine',windowPtr, prefs.stimColor-20, newx+x, newy+y, newx-x, newy-y, prefs.lineWidth);
            
            Screen('Flip', windowPtr);
            pause(0.2);
        end
        pause;
    end
    
    
    % ========== RESPONSE BUTTONS ==========
    texty = screenCenter(2) - 150;
    Screen('DrawText',windowPtr,'Your task is to indicate which blob changed orientation',textx,texty,[255 255 255]); texty = texty + 5*dy;
    Screen('DrawText',windowPtr,'using the number keypad',textx,texty,[255 255 255]); texty = texty + 3*dy;
    
    % show the gabors again but with numbers instead of gabors
    deg1 = round(rand(1,setSize).*180);
    shift = 50;
    X = X + shift;
    Y = Y + shift;
    numbers = {'1','7','9','3'};
    drawfixation(windowPtr,screenCenter(1)+shift,screenCenter(2)+yshift+shift,prefs.fixColor,prefs.fixLength);
    for istim = 1:setSize;
        cuesrcrect = [0 0 StimSize];
        destrect = CenterRectOnPoint(cuesrcrect,X(istim),Y(istim));
        Screen('DrawTexture',windowPtr,StimPatch,cuesrcrect,destrect,180-deg1(istim));
        Screen('DrawText',windowPtr,numbers{istim},X(istim)-10,Y(istim)-10,[255 255 255]);
    end
    
    Screen('Flip', windowPtr);
    pause;
    
    % ========== CONFIDENCE RESPONSE BUTTONS ===========
    texty = screenCenter(2) - 150;
    Screen('DrawText',windowPtr,'You will also report your confidence in your choice',textx,texty,[255 255 255]); texty = texty + 5*dy;
    
    keys = {'A','S','D','F'};
    conf = {'low',' ',' ', 'high'};
    keysize = 75;
    shift = 100;
    x1 = screenCenter(1)-200;
    y1 = screenCenter(2)+shift;
    for ikeys = 1:4;
        Screen('DrawText',windowPtr,keys{ikeys},x1+keysize/2.7,y1+keysize/2.7,255*ones(1,3));
        Screen('DrawText',windowPtr,conf{ikeys},x1+keysize/5,y1+100,255*ones(1,3));
        Screen('FrameRect',windowPtr, 255*ones(1,3),[x1 y1 x1+keysize y1+keysize]);
        x1 = x1 + shift;
    end
    Screen('Flip', windowPtr);
    pause;
    
    % ========== ANY QUESTIONS? ==========
    textx = screenCenter(1) - 75;
    texty = screenCenter(2) - 70;
    Screen('DrawText',windowPtr,'Any questions?',textx,texty,[255 255 255]);
    Screen('Flip', windowPtr);
    pause;
end

% ========== RUNNING TASK FOR A LITTLE ==========
Screen('Flip', windowPtr);
Exp_ChangeLoc(subjid, 'Delay', sessionnum, nTrialsPerCond,1);

end

% =====================================================================
% HELPER FUNCTIONS
% =====================================================================

function drawfixation(windowPtr,x,y,color,size)
Screen('DrawLine',windowPtr,color,x-size,y,x+size,y,2);
Screen('DrawLine',windowPtr,color,x,y-size,x,y+size,2);
end