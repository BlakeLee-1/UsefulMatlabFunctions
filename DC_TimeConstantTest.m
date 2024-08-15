function DC_TimeConstantTest(DAQ,setCurrent)
% Square Wave DC Tests
SRS830_Thermometer = OpenGPIBObject_BufferSize(DAQ.SRS830_Thermometer_gpib);
SRS830_Heater      = OpenGPIBObject_BufferSize(DAQ.SRS830_Heater_gpib);
K2182_Thermometer  = OpenGPIBObject(DAQ.K2182_Thermometer_gpib);
K6220_Thermometer  = OpenGPIBObject(DAQ.K6220_Thermometer_gpib);
K6221_Heater       = OpenGPIBObject(DAQ.K6221_Heater_gpib);
LS340              = OpenGPIBObject(DAQ.LS340_gpib);
K6221_WaveModeSetup(K6221_Heater);

%Wrote some quick programs stored in the InstrumentCommands folder

K6221_Output(K6221_Heater, 'ON')
myhandle = msgbox('Stop the measurement?'); 
ii = 1;
MyClock = tic;
datacell.Current_Thermometer = 1.0e-4;
K6221_DC_SetAmplitude(K6220_Thermometer, datacell.Current_Thermometer)
K6221_Output(K6220_Thermometer, 'ON');
realstarttime = toc(MyClock);
figure; xlabel('Time (s)'); ylabel('Voltage (\muV)'); title('Square Wave Response'); box on; grid on;
datacell.setCurrent = setCurrent;
datacell.setTime = 10;
datacell.Time = [];
datacell.Thermometer_Voltage = [];
datacell.sectionvec = [];
ll = 1;
while ishandle(myhandle)
    
    K6221_DC_SetAmplitude(K6221_Heater, datacell.setCurrent)
    timestart = toc(MyClock);
    datacell.sectionvec(ll) = ii;
    ll = ll + 1;
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
    datacell.sectionvec(ll) = ii;
    ll = ll + 1;
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

%Start to fit the data
intervals = length(datacell.sectionvec);


for i = 1:intervals
    if i ~= length(datacell.sectionvec)
        sections(i).Time = datacell.Time(datacell.sectionvec(i):datacell.sectionvec(i+1));
        sections(i).Voltage = datacell.Thermometer_Voltage(datacell.sectionvec(i):datacell.sectionvec(i+1));
        minnT = min(sections(i).Time);
        minnV = min(sections(i).Voltage);
        for j = 1:length(sections(i).Time)
            sections(i).Time(j) = sections(i).Time(j) - minnT;
            sections(i).Voltage(j) = sections(i).Voltage(j) - minnV;
        end
    elseif i == length(datacell.sectionvec)
        sections(i).Time = datacell.Time(datacell.sectionvec(i):end);
        sections(i).Voltage = datacell.Thermometer_Voltage(datacell.sectionvec(i):end);
        minnT = min(sections(i).Time);
        minnV = min(sections(i).Voltage);
        for j = 1:length(sections(i).Time)
            sections(i).Time(j) = sections(i).Time(j) - minnT;
            sections(i).Voltage(j) = sections(i).Voltage(j) - minnV;
        end
    end
end
colorz = lines(intervals);
figure; xlabel('Time (s)'); ylabel('Voltage (\muV)'); title('Square Wave Response'); box on; grid on;

for i = 1:intervals
    
    hold on;
    plot(sections(i).Time,sections(i).Voltage, 'o', 'Color', colorz(i,:), 'MarkerFaceColor', colorz(i,:), 'DisplayName', ['Section ', num2str(i)])
    sections(i).Time
    
end

 x0 = [1 1 1]; 
timeconst = [];
for i = 1:length(sections)

    if mod(i,2) == 1
        x = sections(i).Time;
        y = sections(i).Voltage;
        fitfun = fittype( @(a,b,d,x) a*exp(x./b)+d);
        [fitted_curve,gof] = fit(x',y',fitfun,'StartPoint',x0);
        coeffvals = coeffvalues(fitted_curve);
        timeconst = [timeconst,coeffvals(2)];

    elseif mod(i,2) == 0
        x = sections(i).Time;
        y = sections(i).Voltage;
        fitfun = fittype( @(a,b,d,x) a*(1-exp(x./b))+d);
        [fitted_curve,gof] = fit(x',y',fitfun,'StartPoint',x0);
        coeffvals = coeffvalues(fitted_curve);
        timeconst = [timeconst,coeffvals(2)];
    end

end
timeconst
mean(timeconst)
std(timeconst)