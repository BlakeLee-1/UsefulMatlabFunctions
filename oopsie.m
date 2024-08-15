function oopsie(DAQ)

    Current_heater = 1.5e-3;
    Current_Thermometer = 1e-4;
    SRS830_Thermometer = OpenGPIBObject_BufferSize(DAQ.SRS830_Thermometer_gpib);
SRS830_Heater      = OpenGPIBObject_BufferSize(DAQ.SRS830_Heater_gpib);
K2182_Thermometer  = OpenGPIBObject(DAQ.K2182_Thermometer_gpib);
K6220_Thermometer  = OpenGPIBObject(DAQ.K6220_Thermometer_gpib);
K6221_Heater       = OpenGPIBObject(DAQ.K6221_Heater_gpib);
% LS340              = OpenGPIBObject(DAQ.LS340_gpib);
K6221_WaveModeSetup(K6221_Heater);

    FrequencySweepVector = [0.1,0.2,0.3,0.5,1,2,3,5];
    for i = 1:length(FrequencySweepVector)
        cd('C:\Users\LeeLab Laptop 2\Documents\Blake\HeatCapacity\HomeMeasurements\FrequencySweeps')
        load(['PreFreqSweep_',num2str(FrequencySweepVector(i)),'Hz.mat'])
        datacell(i).X_Heater = Xheater;
        datacell(i).Y_Heater = Yheater;
        datacell(i).X_Thermometer = Xcernox;
        datacell(i).Y_Thermometer = Ycernox;
        datacell(i).Time = Timee;
        datacell(i).SetAmplitude = Current_heater;
        datacell(i).Current_Thermometer = Current_Thermometer;
        datacell(i).Frequency = FrequencySweepVector(i);
        datacell(i).BathTemperature = 295;
        datacell(i).Thermometer_Voltage_DC = CX_Voltage;
        datacell(i).CXCell = DAQ.CXCell;
        %Now that that's all loaded, I need some of the other pieces

        PointsPerFrequency = DAQ.PointsPerFrequency;


        freq = FrequencySweepVector(i);
        datacell(i).freq = freq;
        filename = ['PreFreqSweep_', num2str(freq),'Hz.mat'];
        
        CXCell = DAQ.CXCell;
        SensitivityCurve = DAQ.SensitivityCurve;
        BathChannel = DAQ.BathChannel;
        R_unit = 'SRDG?';
        
        datacell(i).Current_Thermometer = DAQ.Current_Thermometer;
        
        Current_Thermometer = datacell(i).Current_Thermometer;
        datacell(i).Current_Heater      = DAQ.Current_Heater;
        Current_Heater = datacell(i).Current_Heater;
        datacell(i).LI_Gain = DAQ.LI_Gain;
        datacell(i).PointsPerFrequency = PointsPerFrequency; 
        K6221_WaveMode_SetAmplitude(K6221_Heater,Current_Heater);
        K6221_WaveMode_SetFrequency(K6221_Heater,freq);
        K6221_WaveMode_OnOff(K6221_Heater,'On');
        
        
        pause(7*5); 

        datacell(i).LI_TimeConstant = SRS830_FindTimeConstant(freq);
        SRS830_Set_Sens_Tau(SRS830_Heater,'',datacell(i).LI_TimeConstant,'Time Constant Only')
        SRS830_Set_Sens_Tau(SRS830_Thermometer,'',datacell(i).LI_TimeConstant,'Time Constant Only');
        
        datacell(i).Heater_LI_Sensitivity = 1; %Don't have it
        datacell(i).Heater_LI_Timeconstant = datacell(i).LI_TimeConstant; %Don't have it
        datacell(i).Thermometer_LI_Sensitivity = 20*10^-6; %Don't have it
        datacell(i).Thermometer_LI_Timeconstant = datacell(i).LI_TimeConstant; %Don't have it
        datacell(i).HeaterAmplitude = datacell(i).SetAmplitude; %Have it, same for all at 1.5mA
        datacell(i).R_Heater  = sqrt(datacell(i).X_Heater.^2+datacell(i).Y_Heater.^2); %Have it
        datacell(i).R_Thermometer  = sqrt(datacell(i).X_Thermometer.^2+datacell(i).Y_Thermometer.^2); %Have it
        datacell(i).Thermometer_Resistance = datacell(i).Thermometer_Voltage_DC./datacell(i).Current_Thermometer; %Have it
        RoughSensitivity = SensitivityCurve(datacell(i).Thermometer_Resistance(1)/datacell(i).Current_Thermometer);

        datacell(i).SampleTemperature = datacell(i).CXCell{1}(datacell(i).Thermometer_Resistance); %Have it
        datacell(i).TOsc_Rough = datacell(i).R_Thermometer./(RoughSensitivity.*datacell(i).Current_Thermometer.*datacell(i).LI_Gain); %Don't have rough sensitivity
        datacell(i).VOsc_Times_Frequency = datacell(i).R_Thermometer.*datacell(i).Frequency./datacell(i).LI_Gain; %Have it
    

    end
    for i = 1:8
        switch datacell(i).LI_TimeConstant
            case 10
                datacell(i).TimeConstant = 1
    
            case 11
                datacell(i).TimeConstant = 3
    
            case 12
                datacell(i).TimeConstant = 10
    
            case 13
                datacell(i).TimeConstant = 30
        end
    end
end
