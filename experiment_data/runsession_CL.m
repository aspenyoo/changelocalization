function runsession_CL(subjid,sessionnum) 
% run a session of the second change detection task (CD2)

clc;

if nargin < 1; subjid = input('enter subject ID: ', 's'); end
if nargin < 2; sessionnum = input('enter session number: '); end

beepVec = [1200 1000 800 1000 1200];

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

switch sessionnum
    case 1; nTrialsPerCond = 3;
    case 2; nTrialsPerCond = 3;
    case 3; nTrialsPerCond = 3;
    case 4; nTrialsPerCond = 3;
    case 98; nTrialsPerCond = 1;
end

% ========== TRAINING SESSION ON FIRST DAY ==========
% if sessionnum == 1; 
%     switch exptype1
%         case 'Delay'
%             train_loc_delay(subjid,sessionnum,1); % 30 trials of training
%         case 'Contrast'
%             train_loc_contrast(subjid,sessionnum,1);
%     end
%     
%     Snd('Open');
%     for i = 1:length(beepVec);
%         beep = MakeBeep(beepVec(i),0.1);
%         Snd('Play',beep);
%     end
% 
% % ========== TASK TO GET CONTRAST VALUE =========
% 
% get_contrastresponsefunction(subjid)
% end
Screen('Flip', windowPtr);

% ========== FIRST TASK (CONTRAST AND DELAY COUNTERBALANCED) ============

Exp_ChangeLoc(subjid, exptype1, sessionnum, nTrialsPerCond)
Screen('Flip', windowPtr);

% ========== BREAK SCREEN TO MAKE THE PARTICIPANT MOVE =============
textx = 400;
texty = h/2 - 150;
Screen('DrawText',windowPtr,'Please find your experimenter for caffeine.',textx,texty,[255 255 255]);
Screen('Flip', windowPtr);
pause;

% =========== SECOND TASK (DELAY AND CONTRAST COUNTERBALANCED) ==========

Exp_ChangeLoc(subjid, exptype2, sessionnum, nTrialsPerCond)
Screen('Flip', windowPtr);

% =========== THANK YOU SCREEN ==========
Screen('DrawText',windowPtr,'Thanks for participating!',textx,texty,[255 255 255]);

Screen('Flip', windowPtr);
pause;

sca;