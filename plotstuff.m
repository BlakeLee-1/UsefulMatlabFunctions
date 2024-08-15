load('C:\Users\LeeLab Laptop 2\Documents\Blake\HeatCapacity\HomeMeasurements\FrequencySweeps\datacell.mat');
clz = lines(8);
figure
for i = 1:length(datacell)
hold on; grid on; box on;
plot(datacell(i).Time(1:10000),10^6*(datacell(i).Thermometer_Voltage_DC(1:10000)...
    - mean(datacell(i).Thermometer_Voltage_DC(1:10000))), 'o', 'Color',...
    clz(i,:), 'MarkerFaceColor', clz(i,:), 'DisplayName', ...
    [num2str(datacell(i).Frequency), 'Hz DC'])
plot(datacell(i).Time(1:10000),10^6*datacell(i).R_Thermometer(1:10000)*sqrt(2),...
    'o', 'Color',...
    clz(i,:), 'MarkerFaceColor', 'k', 'DisplayName', ...
    [num2str(datacell(i).Frequency), 'Hz Amp'])
xlabel('Time(s)'); ylabel('Voltage(\muV)');title('Thermometer Voltage Response');


end
figure
for i = 1:1
hold on; grid on; box on;
plot(datacell(i).Time(1:90000),10^6*(datacell(i).Thermometer_Voltage_DC(1:90000)...
    - mean(datacell(i).Thermometer_Voltage_DC(1:10000))), 'o', 'Color',...
    clz(i,:), 'MarkerFaceColor', clz(i,:), 'DisplayName', ...
    [num2str(datacell(i).Frequency), 'Hz DC'])
plot(datacell(i).Time(1:90000),10^6*datacell(i).R_Thermometer(1:90000)*sqrt(2),...
    'o', 'Color',...
    clz(i,:), 'MarkerFaceColor', 'k', 'DisplayName', ...
    [num2str(datacell(i).Frequency), 'Hz Amp'])
xlabel('Time(s)'); ylabel('Voltage(\muV)');title('Thermometer Voltage Response');


end


figure
delV = [];
freqs = [];
for i = 1:length(datacell)
delV = [delV, mean(datacell(i).R_Thermometer)];
freqs = [freqs, datacell(i).Frequency];
end
plot(freqs, delV*10^6*sqrt(2), '-ob', 'MarkerFaceColor', 'b', 'DisplayName', '\DeltaV v Freqs')
xlabel('Frequency (Hz)'); ylabel('Delta Voltage(\muV)');title('Thermometer Voltage Response v Frequency');


for i = 1:length(datacell)
curveFitter(datacell(i).Time, datacell(i).Thermometer_Voltage_DC)
end

