% Sarah, Ian
% 3/21/2022
% Measure Thermal Conductivity at a fixed temperature.
function Sarah_KappaOnePoint_v1(DAQ)

%% Premeasurement Declarations
% Open GPIB Objects for measurement.
K2182a_Heater = OpenGPIBObject(DAQ.K2182a_Htr_gpib);
ls340 = OpenGPIBObject(DAQ.ls340_gpib);
CurrentSource = OpenGPIBObject(DAQ.CurrentSrc_gpib);

% Migrate multiple structure calls to locals:
CXcell = DAQ.CXcell;
HotChannel = DAQ.HotChannel;
ColdChannel = DAQ.ColdChannel;
BathChannel = DAQ.BathChannel;
SamplesPerStep = DAQ.SamplesPerStep;

% Parameters
R_unit = 'SRDG?';

% Make a filename.
BathResistance = read340_obj(ls340,BathChannel,R_unit);
CurrentTemperature = CXcell{3}(BathResistance);
filename = GenerateFilename(DAQ,CurrentTemperature);

% Make vector of Current
CurrentVector = GenerateCurrentVector(DAQ,CurrentTemperature);

%% Start Measurement Sequence

% Measurement parameters to save in datacell.
datacell.SampPerStep = SamplesPerStep;

% Make a figure to monitor data acquisition.
DataPlot = figure('Position',[-845.4 -175.8 840 928]);

% Print out when the current measurement started.
display(['Measuring at ',datestr(now,'mm_dd_yyyy_HH_MM_SS')]);

% Start current
yokoampset_obj(0,CurrentSource);
yoko_on_off_obj(CurrentSource,'ON');

% Initialize loop
ii =0; tic;
for ActiveCurrent = CurrentVector
    yokoampset_obj(ActiveCurrent,CurrentSource);
    for k = 1:SamplesPerStep
        ii=ii+1;
        
        datacell.HeaterCurrent(ii) = ActiveCurrent;
        datacell.HeaterVoltage(ii) = readonevoltage(K2182a_Heater);
        datacell.Time(ii) = toc;
        datacell.HotRes(ii) = read340_obj(ls340,HotChannel,R_unit);
        datacell.ColdRes(ii) = read340_obj(ls340,ColdChannel,R_unit);
        datacell.BathRes(ii) = read340_obj(ls340,BathChannel,R_unit);
        datacell.HotTemp(ii) = CXcell{1}(datacell.HotRes(ii));
        datacell.ColdTemp(ii) = CXcell{2}(datacell.ColdRes(ii));
        datacell.BathTemp(ii) = CXcell{3}(datacell.BathRes(ii));
        datacell.DeltaT(ii) = datacell.HotTemp(ii) - datacell.ColdTemp(ii); 
        datacell.HeaterPower(ii) = datacell.HeaterCurrent(ii).^2.*1000;
        datacell.SampleAverageTemp(ii) = (datacell.HotTemp(ii)+...
            datacell.ColdTemp(ii))./2;
        
        if mod(ii,100)==0
            figure(DataPlot);
            % System Temperatures
            subplot(3,2,1); 
            plot(datacell.Time,datacell.BathTemp,'-b','DisplayName','Bath'); 
            hold on; plot(datacell.Time,datacell.SampleAverageTemp,'-r',...
                'DisplayName','Sample Average'); hold off; 
            Label_Plot('Time','s','T','K'); 
            
            subplot(3,2,2); 
            yyaxis left; hold off; plot(datacell.Time,datacell.HeaterCurrent,'-k',...
                'DisplayName','Heater Current'); 
            Label_Plot('Time','s','I_{htr}','A'); 
            yyaxis right; hold off; plot(datacell.Time,datacell.HeaterVoltage,'-r',...
                'DisplayName','Heater Voltage'); 
            Label_Plot('Time','s','V_{htr}','V'); 
         
            subplot(3,2,3); 
            plot(datacell.Time,datacell.HotRes,'-r','DisplayName','R_{Hot}'); 
            hold on; plot(datacell.Time,datacell.ColdRes,'-b',...
                'DisplayName','R_{cold}'); hold off; 
            Label_Plot('Time','s','R_{cernox}','\Omega'); 

            subplot(3,2,4); 
            plot(datacell.Time,datacell.HotTemp,'-r','DisplayName','T_{Hot}'); 
            hold on; plot(datacell.Time,datacell.ColdTemp,'-b',...
                'DisplayName','T_{cold}'); hold off; 
            Label_Plot('Time','s','T_{cernox}','K'); 
            
            subplot(3,2,5); 
            plot(datacell.Time,datacell.DeltaT,'-b','DisplayName','Delta T'); 
            hold off; 
            Label_Plot('Time','s','\Delta T','K'); 
            
            subplot(3,2,6); 
            plot(datacell.HeaterPower,datacell.DeltaT,'-ob','DisplayName','Delta T'); 
            hold off; 
            Label_Plot('Power','W','\Delta T','K'); 
            drawnow
        end
    end
    
end

% % Try to close the figure if it is still open. 
try
    close(DataPlot)
end

% Save our data
save(filename,'-STRUCT','datacell');
yokoampset_obj(0,CurrentSource)
yoko_on_off_obj(CurrentSource,'OFF');



end

function filename = GenerateFilename(DAQ,CurrentTemperature)
% DAQ is the unaltered DAQ input.
% CurrentTemperature -- String.
CurrentTemperature = num2str(round(CurrentTemperature,2));
filename = [DAQ.Tswpfileroot,'\',regexprep([DAQ.SampleInfoStr,'_',...
    CurrentTemperature,'K_',datestr(now,'mm_dd_yyyy_HH_MM_SS')],...
    '\.','p'),'.mat'];
end

function CurrentVector = GenerateCurrentVector(DAQ,CurrentTemperature)

% Defining local variables from DAQ structure
HeaterCurrentStop = DAQ.HeaterCurrentStopFunc(CurrentTemperature);  % Sets ending current value
HeaterCurrentStart = DAQ.HeaterCurrentStart;                        % Sets starting current value
HeaterCurrentSteps = DAQ.HeaterCurrentSteps;                        % Sets current step size
CurrentVector = sqrt([linspace(HeaterCurrentStart.^2,...            % Create current vector that scales linearly with I^2 and ends with zero
    HeaterCurrentStop.^2,HeaterCurrentSteps+1),0]);

% Too high current protection.
if HeaterCurrentStop>4.5e-3
    error('Too big current');
end

end
