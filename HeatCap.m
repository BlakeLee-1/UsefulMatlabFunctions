function HeatCap()
    M = LoadData();
    hold off;
    fieldd = [5,10,18,14]
    for i = 1:4
        Mp = calcC(M(i));
        Md = plotData(Mp, [int2str(fieldd(i)),'T']);
    end
end

function [M] = LoadData()
    grid on;
    fileroot = '/home/blake/Documents/MATLAB/SCM2/Angleswp/Just1p8/'; %Parent file location for our data
    cd(fileroot); %We need to get there
    filedata = dir('*.mat'); 
    filedata
    for i = 1:length(filedata)  %for all of these elements
       M(i) = load(filedata(i).name);  % Load all data first before adding fields so they are all the same length
    end
    for i = 1:length(filedata) %straight up addressing all of those named elements from last loop and giving them these attributes
       M(i).name = filedata(i).name;
       M(i).Fileroot = fileroot;
       M(i).DateStr = filedata(i).date;
       M(i).Datenum = filedata(i).datenum;
    end
   
   
end

function [M] = calcC(M)
load("/home/blake/Documents/MATLAB/SCM2/SensitiveityFitCX_E.mat");
    for i = 1:length(M.Time)
        M.HeatCap(i) = M.Current_Heater*M.R_Heater(i)/(4*pi*M.Frequency*(M.R_Thermometer(i)/M.Current_Thermometer)*abs(SensitivityFit(M.Thermometer_Resistance(i))));

    end
end

function [M] = plotData(M,nam)
%Deriv = smooth(diff(M.HeatCap)./diff(M.SampleTemperature),100);
% for i = 2:length(M.HeatCap)
% deriv = [deriv, (M.HeatCap(i)-M.HeatCap(i-1))/(M.SampleTemperature(i)-M.SampleTemperature(i-1))];
% end
% ld = length(deriv)
% lt = length(M.SampleTemperature)
%angle = M.RotatorNumber/102
%plot(M.RotatorNumber./108, M.HeatCap - (max(M.HeatCap)+min(M.HeatCap))/2, 'o', 'DisplayName', nam)
curveFitter(M.RotatorNumber./103,M.HeatCap - (max(M.HeatCap)+min(M.HeatCap))/2)
hold on;
end
