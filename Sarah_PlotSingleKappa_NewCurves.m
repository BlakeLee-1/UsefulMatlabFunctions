function Sarah_PlotSingleKappa_NewCurves(InputCurves)

datacell = load(uigetfile('.mat','C:\Users\LeeLabLaptop\Documents\Sarah\KVO_Nitrogen\'));
% DataPlot = figure('Position',[-845.4 -175.8 840 928]);
DataPlot = figure; 
SampleAverageTemp = mean(datacell.SampleAverageTemp);
filename = strcat('KVO_',strrep(num2str(SampleAverageTemp),'.',','),'K_Measurement_Plots');  %%% this keeps cutting off the "_Measurement_Plots" part
figure('name',filename);
% System Temperatures
subplot(3,2,1);
SampleAverageTemp = (InputCurves{1}(datacell.HotRes)+InputCurves{2}(datacell.ColdRes))./2;
plot(datacell.Time,datacell.BathTemp,'-b','DisplayName','Bath');
hold on; plot(datacell.Time,SampleAverageTemp,'-r',...
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
plot(datacell.Time,InputCurves{1}(datacell.HotRes),'-r','DisplayName','T_{Hot}');
hold on; plot(datacell.Time,InputCurves{2}(datacell.ColdRes),'-b',...
    'DisplayName','T_{cold}'); hold off;
Label_Plot('Time','s','T_{cernox}','K');

DeltaT = InputCurves{1}(datacell.HotRes)-InputCurves{2}(datacell.ColdRes);
subplot(3,2,5);
plot(datacell.Time,DeltaT,'-b','DisplayName','Delta T');
hold off;
Label_Plot('Time','s','\Delta T','K');

subplot(3,2,6);
plot(datacell.HeaterPower,DeltaT,'-ob','DisplayName','Delta T');
hold off;
Label_Plot('Power','W','\Delta T','K');

end