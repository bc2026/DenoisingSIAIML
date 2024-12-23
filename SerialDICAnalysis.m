%% Run a serial analysis on DIC images
%YW 06/16/21
%The entire analysis only uses auto-tracking mode
%% INPUT (the absolute path of an excel sheet with X/Y coordinates and stream path)
%Excel Columns
%Stream:the absolute path of the stream
%Particle ID: depend on the file naming system (better to use a systematic
%naming)
%X: X coordinate of the particle
%Y: Y coordinate of the particle
%SFrame: Starting frame of tracking (default value is 1 if left empty)
%EFrame: Ending frame of tracking (default value is inf if left empty)

ExcelPath = '/Users/Yuhao/Desktop/TeriOdom/Experiments/DICAnalysis/RunFile20221119_LD_N87.xlsx';

%% Load Excel sheet
Table = readtable(ExcelPath,'Format','auto');
Table = Table(~ismissing(Table(:,1)),:);
trials = height(Table); 
RunTime = zeros(1,trials); %used to estimate the run time
Message = sprintf('Total %d runs.',trials);
disp(Message);
JumpFrames = zeros(trials,1);
%% Loop analysis
for i = 1:trials
    tic
    Path = Table.Stream{i};
    ID = Table.ParticleID{i};
    X = Table.X(i);
    Y = Table.Y(i);
    SFrame = Table.SFrame(i);
    if isnan(SFrame)
        SFrame = 1;
    end
    EFrame = Table.EFrame(i);
    if isnan(EFrame)
        EFrame = inf;
    end
    JumpFrames(i) = DICAnalysis(Path,ID,X,Y,SFrame,EFrame);
    RunTime(i) = toc;
    AvgRunTime = mean(RunTime(1:i));
    EstimatedFinishTime = round(AvgRunTime*(trials-i)/60);
    Message = sprintf('%d runs finished, %d trials left. Estimated time needed %.0f mins.',i,trials-i,EstimatedFinishTime);
    disp(Message);
end
%convert data into a table
T_JumpFrames = table(JumpFrames);
%concatenate the tables
NewT = [Table,T_JumpFrames];
%write the new table into the existing excel
writetable(NewT,ExcelPath,'WriteMode','overwritesheet')
   



