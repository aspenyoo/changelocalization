function runsession_CL(subjid,sessionnum,redo)
% run a session of the second change detection task (CD2)
% redo: if someone made an error and you want to pick up where you left
% off. best for first session (other ones you don't need to redo to pick up
% where you left off)

clc;

if nargin < 1; subjid = input('enter subject ID: ', 's'); end
if nargin < 2; sessionnum = input('enter session number: '); end
if nargin < 3; redo = 0; end

screenNumber = max(Screen('Screens'));       % use external screen if exists
[w, h] = Screen('WindowSize', screenNumber);  % screen resolution of smaller display
windowPtr = Screen('OpenWindow',screenNumber,128*ones(1,3),[],32,2);
Screen(windowPtr,'BlendFunction',GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
Screen('TextSize',windowPtr,20);
Screen('TextFont',windowPtr,'Helvetica');
HideCursor;

% pick out whether contrast or delay is first
if mod(sessionnum,2)
    exptype1 = 'Contrast';
    exptype2 = 'Delay';
else
    exptype1 = 'Delay';
    exptype2 = 'Contrast';
end

% ========== TRAINING SESSION ON FIRST DAY ==========
textx = w/2-300;
texty = h/2 - 150;
    
if sessionnum == 1;
    training_getcontrast(subjid)

else
    expsession = sessionnum - 1;
    if (expsession == 1); nTrialsPerCond = 10; else nTrialsPerCond = []; end
    nSessions = 6;
    
    % ========== FIRST TASK (CONTRAST AND DELAY COUNTERBALANCED) ============
    % practice
    if (expsession == 1);
        Exp_ChangeLoc(subjid, ['Pract' exptype1], expsession, 1, 6)
        Screen('Flip', windowPtr);
    end
    
    Exp_ChangeLoc(subjid, exptype1, expsession, nTrialsPerCond, nSessions,1-redo)
    Screen('Flip', windowPtr);
    
    % ========== BREAK SCREEN TO MAKE THE PARTICIPANT MOVE =============
    
    Screen('DrawText',windowPtr,'You have completed the first task.',textx,texty,[255 255 255]);
    Screen('DrawText',windowPtr,'Please find your experimenter.',textx,texty+50,[255 255 255]);
    Screen('Flip', windowPtr);
    pause;
    
    % =========== SECOND TASK (DELAY AND CONTRAST COUNTERBALANCED) ==========
    % practice
    if (expsession == 1);
        Exp_ChangeLoc(subjid, ['Pract' exptype2], expsession, 1, 6)
        Screen('Flip', windowPtr);
    end
    
    Exp_ChangeLoc(subjid, exptype2, expsession, nTrialsPerCond, nSessions,1-redo)
    Screen('Flip', windowPtr);
end

% =========== THANK YOU SCREEN ==========
Screen('DrawText',windowPtr,'Thanks for participating!',textx,texty,[255 255 255]);

Screen('Flip', windowPtr);
pause;


sca;