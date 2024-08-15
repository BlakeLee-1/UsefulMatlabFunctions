% Blake Lee
% 12022/09/16
% Heat Capacity Pre Frequency Sweep Analysis

function HeatCapacity_PreFrequencySweep(DAQ)

savefolder = 'C:\Users\LeeLab Laptop 2\Documents\Blake\HeatCapacity\HomeMeasurements\FrequencySweeps';

% Open GPIB Objects and extract DAQ measurement parameters.
SRS830_Thermometer = OpenGPIBObject_BufferSize(DAQ.SRS830_Thermometer_gpib);
SRS830_Heater      = OpenGPIBObject_BufferSize(DAQ.SRS830_Heater_gpib);
K2182_Thermometer  = OpenGPIBObject(DAQ.K2182_Thermometer_gpib);
K6220_Thermometer  = OpenGPIBObject(DAQ.K6220_Thermometer_gpib);
K6221_Heater       = OpenGPIBObject(DAQ.K6221_Heater_gpib);
LS340              = OpenGPIBObject(DAQ.LS340_gpib);
K6221_WaveModeSetup(K6221_Heater);

FrequencySweepVector = [0.08,0.03];
PointsPerFrequency = DAQ.PointsPerFrequency;
for i = 1:length(FrequencySweepVector)


freq = FrequencySweepVector(i);
datacell.freq = freq;
filename = ['PreFreqSweep_', num2str(freq),'Hz.mat'];

CXCell = DAQ.CXCell;
SensitivityCurve = DAQ.SensitivityCurve;
BathChannel = DAQ.BathChannel;
R_unit = 'SRDG?';

datacell.Current_Thermometer = DAQ.Current_Thermometer;

Current_Thermometer = datacell.Current_Thermometer;
datacell.Current_Heater      = DAQ.Current_Heater;
Current_Heater = datacell.Current_Heater;
datacell.LI_Gain = DAQ.LI_Gain;
LI_Gain = datacell.LI_Gain;
datacell.PointsPerFrequency = PointsPerFrequency; 
K6221_WaveMode_SetAmplitude(K6221_Heater,Current_Heater);
K6221_WaveMode_SetFrequency(K6221_Heater,freq);
K6221_WaveMode_OnOff(K6221_Heater,'On');



LI_TimeConstant = SRS830_FindTimeConstant(freq);
datacell.LI_TimeConstant = SRS830_FindTimeConstant(freq);

SRS830_Set_Sens_Tau(SRS830_Heater,'',LI_TimeConstant,'Time Constant Only')
SRS830_Set_Sens_Tau(SRS830_Thermometer,'',LI_TimeConstant,'Time Constant Only');
pause(round(1/freq*3*5)); 
%



% Measure preliminary temperature for filename;
% BathResistance = read340_obj(LS340,BathChannel,R_unit);
% CurrentTemperature = CXCell{2}(BathResistance);
%Filename = GenerateFilename(DAQ,CurrentTemperature);

% Calculate a rough sensitivity for to estimate Tosc during measurement. 
Thermometer_Resistance = readonevoltage(K2182_Thermometer)./Current_Thermometer;
RoughSensitivity = SensitivityCurve(Thermometer_Resistance); % In ohm/K. To get gradient, divide resistance reating by sensitivity.
%43 <-This was just here

% Sweep frequency.
DataFigure = figure;
K6221_DC_SetAmplitude(K6220_Thermometer, datacell.Current_Thermometer)
K6221_Output(K6220_Thermometer, 'ON');

%myhandle = msgbox('Stop the measurement?'); 
    
ii = 1;

titlee = ['Trial Data at ',num2str(freq), 'Hz'];
MyClock = tic;
    while toc(MyClock) < 1/freq*100
    [datacell.Xheater(ii),datacell.Yheater(ii)] = ReadSRS830_XY(SRS830_Heater);
    datacell.CX_Voltage(ii) = readonevoltage(K2182_Thermometer);
    datacell.Timee(ii) = toc(MyClock);
    [datacell.Xcernox(ii),datacell.Ycernox(ii)] = ReadSRS830_XY(SRS830_Thermometer);
    ii = ii + 1;
        if mod(ii,100) == 0
            plot(datacell.Timee, datacell.CX_Voltage, 'ob', 'DisplayName', 'Thermometer Resistance')
            xlabel('Time (s)'); ylabel('Resistance (\Omega)'); hold off; grid on; box on;
            title(titlee)
            drawnow;
        end
    end

filenamefig = [filename(1:end-3),'fig'];
saveas(gcf, filenamefig);

save(filename,'-STRUCT','datacell');

clear('datacell')
end
end
