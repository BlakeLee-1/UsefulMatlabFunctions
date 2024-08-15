% Sarah
% 3/29/2022
% KVO Data Analysis log

%% Plot a single set of raw data measurements
Sarah_PlotSingleKappa()

%% Checking MS(850) data delta T vs time, delta T vs current

figure('name','MS850 time check');
   
            subplot(1,2,1); 
            plot(MS(850).Time,MS(850).DeltaT,'or');
            Label_Plot('Time','s','Delta T','K');
            
            subplot(1,2,2); 
            plot(MS(850).HeaterPower,MS(850).DeltaT,'ob');
            Label_Plot('Power','W','Delta T','K');
            
            drawnow
            
%% Checking MS(850) data 
figure('name','MS850');
   
            subplot(1,3,1); 
            plot(MS(850).HeaterPower,MS(850).DeltaT,'or');
            Label_Plot('Power','W','Delta T','K');
            
            subplot(1,3,2); 
            plot(MS(850).HeaterPowerStepsOnly,MS(850).DeltaTStepsOnly,'ob');
            Label_Plot('Power','W','Delta T','K');
            
            subplot(1,3,3); 
            plot(MS(850).HeaterPowerStepsCut,MS(850).DeltaTStepsCut,'og');
            Label_Plot('Power','W','Delta T','K');
 
            drawnow
            
 
%%  Comparing different MS(i) data
figure('name','Sorted Comparing Currents');
   
            plot(MS(1).HeaterPowerStepsCut,MS(1).DeltaTStepsCut,'or'); hold on;
            plot(MS(150).HeaterPowerStepsCut,MS(150).DeltaTStepsCut,'om'); hold on;
            plot(MS(300).HeaterPowerStepsCut,MS(300).DeltaTStepsCut,'oy'); hold on;
            plot(MS(450).HeaterPowerStepsCut,MS(450).DeltaTStepsCut,'og'); hold on;
            plot(MS(600).HeaterPowerStepsCut,MS(600).DeltaTStepsCut,'oc'); hold on;
            plot(MS(850).HeaterPowerStepsCut,MS(850).DeltaTStepsCut,'ob'); hold on;
            Label_Plot('Power','W','Delta T','K');
 
            drawnow
            
%% Getting sample avg temp and sorting it
for i = length(MS)
    MS(i).TotalAvgerageTemp = mean(MS(i).SampleAverageTemp);     % Create new field for single Total Average Temp per .mat data file
end
TotalAverageTemp = vertcat(MS.TotalAvgerageTemp);                % Create Total Average Temp vector from structure field
[Temp,Index] = sort(TotalAverageTemp);                           % Sort Total Average Temp vector and index values 
MS = MS(Index)

%% Use this to save figures, Saahhh
SavePngFig(figfr,'KVO_300K_Ex100_Sens100_400SPS')

%% Recalibrating Cernox's
% Load the SD cernox calibrations. 
% BathCell = {SD1, SD2}; 

% make new calibrations using x = bath temp, y = hot/cold cernox temps
BathCell = load('C:\Users\LeeLabLaptop\Documents\Minhyea Lee Research\Maglab_March2019\CXfits_March2019\SDCXcell.mat');
BathFit = BathCell.CXcell{1};

deg = 7;
CX_fr = 'C:\Users\LeeLabLaptop\Documents\Sarah\CXData\MN_3_2022\'; 
cd(CX_fr);
TLim = 70; 
Dathot = load('CX_M_CrI3_Hot.csv');  % Hot Cx Fit
x = Dathot(:,1); y = Dathot(:,2);
KeepInds = y>=TLim; 
x = x(KeepInds); y = y(KeepInds);
HotFit = RtoT_cal_inputRT(x,y,'fitRvT',deg,'Yes','Polynomial',[]);
Datcold = load('CX_N_CrI3_Trans.csv');   % Cold Cx Fit
x = Datcold(:,1); y = Datcold(:,2);
KeepInds = y>=TLim; 
x = x(KeepInds); y = y(KeepInds);
ColdFit = RtoT_cal_inputRT(x,y,'plotRvT',deg,'Yes','Polynomial',[]);