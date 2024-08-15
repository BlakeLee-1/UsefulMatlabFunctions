%% 12022/09/14

% Cernox E, Heater 1 
% Plus CsYbSe2 sample, mass 0.5mg. 
% 
% DAQ Initialization
DAQ.Bswpfileroot = 'C:\Users\LeeLab Laptop 2\Documents\Blake\HeatCapacity\HomeMeasurements_WithSample\Bswp'; % This is where field sweep data is saved.
DAQ.Tswpfileroot = 'C:\Users\LeeLab Laptop 2\Documents\Blake\HeatCapacity\HomeMeasurements_WithSample\Tswp'; % This is where T-sweep data is saved.
DAQ.Testsfileroot = 'C:\Users\LeeLab Laptop 2\Documents\Blake\HeatCapacity\HomeMeasurements_WithSample\tests\'; % Test data files go here. 
DAQ.Fswpfileroot = 'C:\Users\LeeLab Laptop 2\Documents\Blake\HeatCapacity\HomeMeasurements_WithSample\FrequencySweeps'; % This is where Frequency sweep data is saved.
DAQ.Aswpfileroot = 'C:\Users\LeeLab Laptop 2\Documents\Blake\HeatCapacity\HomeMeasurements_WithSample\AmplitudeSweeps'; % This is where Amplitude sweep data is saved.
DAQ.SampleInfoStr = 'HeatCapacity_CxE_Htr1_CsYbSe2_HT_CD1'; % Besides datetime, temperature, and field, what do you name the file? 

DAQ.SRS830_Thermometer_gpib = 8;
DAQ.SRS830_Heater_gpib = 9;
DAQ.K2182_Thermometer_gpib = 7; 
DAQ.K6220_Thermometer_gpib = 14;
DAQ.K6221_Heater_gpib      = 12;
DAQ.LS340_gpib = 13;

DAQ.Current_Thermometer = 1e-4; 
DAQ.Current_Heater = 1.5e-3; 
DAQ.LI_Gain = 1; 

DAQ.FrequencySweepVector = [0.08,0.03];
DAQ.PointsPerFrequency = 500;

DAQ.SetHeaterCurrent = 1.5e-3;
DAQ.SetHeaterFrequency = 1; 

DAQ.BathChannel = 'C'; 

% Load Cernox Calibrations
BathCell = load('C:\Users\LeeLab Laptop 2\Documents\Minhyea Lee Research\Maglab_March2019\CXfits_March2019\SDCXcell.mat');
BathFit = BathCell.CXcell{1};

CX_E_Fileroot = 'C:\Users\LeeLab Laptop 2\Documents\Blake\HeatCapacity\CX_E_CalibrationData\CX_E_ComboCal.csv'; 
CX_E_Data = readmatrix(CX_E_Fileroot);  % Hot Cx Fit
Res = CX_E_Data(:,1); Temp = CX_E_Data(:,2);
TLimit = 70; 
KeepInds = Temp>=TLimit; 
Res = Res(KeepInds); Temp = Temp(KeepInds); 
PolynomialFitDegree = 7; 
CX_E_Fit = RtoT_cal_inputRT(Res',Temp','fitRvT',PolynomialFitDegree,...
    'Yes','Polynomial',[]);
%Was 0.005
InterpRes = 44:.05:650;
InterpTemp = CX_E_Fit(InterpRes)'; 
Sensitivity = diff(InterpRes)./diff(InterpTemp); 
SensitivityFit = RtoT_cal_inputRT(InterpRes(1:end-1),abs(Sensitivity),'fitRvT',PolynomialFitDegree,...
    'Yes','Polynomial',[]);
SensitivityFit = @(R) -1.*SensitivityFit(R); 

% {Hot Cernox Calibration Curve, Bath Cernox Calibration Curve}
DAQ.CXCell = {CX_E_Fit,BathFit};
DAQ.SensitivityCurve = SensitivityFit; 


%% Square Wave DC Tests
SRS830_Thermometer = OpenGPIBObject_BufferSize(DAQ.SRS830_Thermometer_gpib);
SRS830_Heater      = OpenGPIBObject_BufferSize(DAQ.SRS830_Heater_gpib);
K2182_Thermometer  = OpenGPIBObject(DAQ.K2182_Thermometer_gpib);
K6220_Thermometer  = OpenGPIBObject(DAQ.K6220_Thermometer_gpib);
K6221_Heater       = OpenGPIBObject(DAQ.K6221_Heater_gpib);
%LS340              = OpenGPIBObject(DAQ.LS340_gpib);
K6221_WaveModeSetup(K6221_Heater);

%Wrote some quick programs stored in the InstrumentCommands folder

K6221_Output(K6221_Heater, 'ON')
myhandle = msgbox('Stop the measurement?'); 
ii = 1;
MyClock = tic;
datacell.Current_Thermometer = 1e-4;
% K6221_DC_SetAmplitude(K6220_Thermometer, datacell.Current_Thermometer)
% K6221_DC_SetAmplitude(K6220_Thermometer, 'ON');
realstarttime = toc(MyClock);
figure; xlabel('Time (s)'); ylabel('Voltage (\muV)'); title('Square Wave Response'); box on; grid on;
datacell.setVoltage = 1.8e-3;
datacell.setTime = 10;
datacell.Time = [];
datacell.Thermometer_Voltage = [];
while ishandle(myhandle)

    K6221_DC_SetAmplitude(K6221_Heater, datacell.setVoltage)
    timestart = toc(MyClock);
   
    while (toc(MyClock) - timestart) < datacell.setTime/2
        datacell.Time(ii) = toc(MyClock);
        datacell.Thermometer_Voltage(ii) = readonevoltage(K2182_Thermometer);
        ii = ii + 1;
    end
    hold off;
    plot(datacell.Time, datacell.Thermometer_Voltage, 'ob', 'DisplayName', 'Thermometer Voltage');
    drawnow;
    while (toc(MyClock) - timestart) < datacell.setTime
        datacell.Time(ii) = toc(MyClock);
        datacell.Thermometer_Voltage(ii) = readonevoltage(K2182_Thermometer);
        ii = ii + 1;
    end
    hold off;
    plot(datacell.Time, datacell.Thermometer_Voltage, 'ob', 'DisplayName', 'Thermometer Voltage');
    drawnow; 
    K6221_DC_SetAmplitude(K6221_Heater, 0)
    timestart = toc(MyClock);
    while (toc(MyClock) - timestart) < datacell.setTime/2
        datacell.Time(ii) = toc(MyClock);
        datacell.Thermometer_Voltage(ii) = readonevoltage(K2182_Thermometer);
        ii = ii + 1;
    end
    hold off;
    plot(datacell.Time, datacell.Thermometer_Voltage, 'ob', 'DisplayName', 'Thermometer Voltage');
    drawnow;
    while (toc(MyClock) - timestart) < datacell.setTime
        datacell.Time(ii) = toc(MyClock);
        datacell.Thermometer_Voltage(ii) = readonevoltage(K2182_Thermometer);
        ii = ii + 1;
    end
    hold off;
    plot(datacell.Time, datacell.Thermometer_Voltage, 'ob', 'DisplayName', 'Thermometer Voltage');
    drawnow; 

end

%Though written here, I just saved this script as DC_TimeConstantTest(DAQ, heatercurrent)

%i wrote the processing code
DC_TimeConstantTest(DAQ,2e-3)
%Time const 1.8809, std 0.63
DC_TimeConstantTest(DAQ,1e-3)
%Time const 1.05, std 0.09


%%12022/09/16
%Time for me to do a frequency sweep. We have the time constant of around
%1.5s, corresponding to a frequency of 0.66Hz. 
HeatCapacity_FrequencySweepSep22(DAQ)

%12022/09/16-19
%Over the weekend I ran a frequency sweep at 2 hours per frequency. A few
%things went wrong and I had to restart a couple times, but I wound up with
%some data. Now I need to see if there are any straggler values I missed so
%that I can perform the calculations Ian did in his version of the code.
%I will use his version in the future so there will be nothing missed, but
%this was mostly practice. Code I ran is below.

HeatCapacity_PreFrequencySweep(DAQ);

%In Ian's Frequency Sweep code, he uses the following calculations
%Note when I say don't have it, I had it, but it wasn't saved and is easily
%recoverable.

datacell.Heater_LI_Sensitivity = Heater_LI_Sensitivity; %Don't have it
datacell.Heater_LI_Timeconstant = LI_TimeConstant; %Don't have it
datacell.Thermometer_LI_Sensitivity = Thermometer_LI_Sensitivity; %Don't have it
datacell.Thermometer_LI_Timeconstant = LI_TimeConstant; %Don't have it
datacell.Frequency = SetFrequency; %Have it
datacell.HeaterAmplitude = SetAmplitude; %Have it, same for all at 1.5mA
datacell.R_Heater  = sqrt(datacell.X_Heater.^2+datacell.Y_Heater.^2); %Have it
datacell.R_Thermometer  = sqrt(datacell.X_Thermometer.^2+datacell.Y_Thermometer.^2); %Have it
%datacell.BathTemperature = CXCell{2}(datacell.BathResistance); %Did not record temperature, will assume it was constant. 340 was borrowed
datacell.Thermometer_Resistance = datacell.Thermometer_Voltage_DC./Current_Thermometer; %Have it
datacell.SampleTemperature = CXCell{1}(datacell.Thermometer_Resistance); %Have it
datacell.TOsc_Rough = datacell.R_Thermometer./(RoughSensitivity.*Current_Thermometer.*LI_Gain); %Don't have rough sensitivity
datacell.VOsc_Times_Frequency = datacell.R_Thermometer.*datacell.Frequency./LI_Gain; %Have it


%All in all, I need to grab 5 quantities, which should only take 10-15
%minutes
oopsie(DAQ) %Write this real quick
        


