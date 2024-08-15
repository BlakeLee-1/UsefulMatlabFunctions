% Sarah and Ian Measurement Log 3/21/2022 K2V3O8

%% DAQ Initialization
DAQ.Bswpfileroot = 'C:\Users\LeeLabLaptop\Documents\Sarah\KVO_Nitrogen\Bswp'; % This is where field sweep data is saved.
DAQ.Tswpfileroot = 'C:\Users\LeeLabLaptop\Documents\Sarah\KVO_Nitrogen\Tswp'; % This is where T-sweep data is saved.
DAQ.Testsfileroot = 'C:\Users\LeeLabLaptop\Documents\Sarah\KVO_Nitrogen\tests\';
DAQ.SampleInfoStr = 'Sarah_KVO'; % Besides datetime, temperature, and field, what do you name the file? 

DAQ.avgTol = .001; % Used for tolerance in PID control to new temperature. 
DAQ.varTol = .001; % % Used for tolerance in PID control to new temperature. 
DAQ.SampleHeaterRange = 3; % What heater power range do you want for PID heater? 
DAQ.TempControlLoop = 1; % Which PID control loop

DAQ.pauseTime = 0; % Add any pause to measurement?
DAQ.SamplesPerStep = 200; % How many data points per heater current. 
DAQ.HeaterCurrentSteps = 3; % How many heater current steps. 
DAQ.HeaterCurrentStopFunc = @(T) 2e-3; % What is the max current used during the last step? Can be fit to a function of temperature to give variable gradients with varying thermal conductivity. 
DAQ.HeaterCurrentStart = 0; % What heater current do you start at?

DAQ.ls340_gpib = 13; % Measurement thermometers, assign channels. 
DAQ.HotChannel = 'A'; %Hot CX M
DAQ.ColdChannel = 'B'; % Cold CX N
DAQ.BathChannel = 'C'; % SD1

DAQ.CurrentSrc_gpib = 11;
DAQ.K2182a_Htr_gpib = 18;

% Load the SD cernox calibrations. 
% BathCell = {SD1, SD2}; 
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
HotFit = RtoT_cal_inputRT(x,y,'plotRvT',deg,'Yes','Polynomial',[]);
Datcold = load('CX_N_CrI3_Trans.csv');   % Cold Cx Fit
x = Datcold(:,1); y = Datcold(:,2);
KeepInds = y>=TLim; 
x = x(KeepInds); y = y(KeepInds);
ColdFit = RtoT_cal_inputRT(x,y,'plotRvT',deg,'Yes','Polynomial',[]);

% DAQ.CXcell = {M Cal, N Cal, SD1}
DAQ.CXcell = {HotFit, ColdFit, BathFit};

optimizeCX = updateCXcruves(DAQ.CXcell,{39.162,38.252,39.425},BathFit(39.425));
DAQ.CXcell = optimizeCX;

%% Start the Measurements here

Sarah_KappaOnePoint_v1(DAQ)


% Upped cernox excitations to 100muA and 100ohm range, recalibrate to bath:
optimizeCX = updateCXcruves(DAQ.CXcell,{39.093,38.34,39.407},BathFit(39.407));
DAQ.CXcell = optimizeCX;

% Check for improved T stability b/c of higher R sensitivity. 
Sarah_KappaOnePoint_v1(DAQ)

% Change samples per step from 200 to 400 to get more points at /Delta T
% equilibration
DAQ.SamplesPerStep = 400;
Sarah_KappaOnePoint_v1(DAQ)

% Test the vacuum gauge. 
% Ambient: 
% VacCell_GaugeTest_KeithleyYoko(Fileroot,InfoString,HeaterCurrent); 
VacCell_GaugeTest_KeithlyYoko(DAQ.Testsfileroot ,'KVO_Ambient',3e-3);
% Input file to extract time constants. \tau_{ambient} = 0.94s; 
VacGauge_Process_Nov2018(DAQ.Testsfileroot,'KVO_Ambient_21_03_22_14_06_15','Fits','New'); 

% Pumped: 
VacCell_GaugeTest_KeithlyYoko(DAQ.Testsfileroot ,'KVO_Pumped',3e-3);
% \tau_{pumped} = 1.65s
VacGauge_Process_Nov2018(DAQ.Testsfileroot,'KVO_Pumped_21_03_22_16_05_01','Fits','New'); 

DAQ.SamplesPerStep = 300;
Temperatures = [300,200,77]; 
Currents     = [1.5e-3,1e-3,.7e-3];
HeaterPoly = polyfit(Temperatures,Currents,2); 
DAQ.HeaterCurrentStopFunc = @(T) polyval(HeaterPoly,T);

DAQ.SamplesPerStep = 300;
MsgHandle = msgbox('Stop Measuring?');
while ishandle(MsgHandle)
    Sarah_KappaOnePoint_v1(DAQ);
end
