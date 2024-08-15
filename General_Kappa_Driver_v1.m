function MS = General_Kappa_Driver_v1(LoadOption,PlotOption,NewCalibrationCurves)
[MS,PP] = LoadKappaData(LoadOption);
MS = ProcessKappaData(MS,PP,NewCalibrationCurves);
switch PlotOption
    case 'Kappa vs T'
        xarg = 'MeanSampleTemperature';
        yarg = 'Kappa';
        figure; hold on; GeneralPlot1(MS,xarg,yarg);
        xlabel('T_{SD2} [K]'); ylabel('\kappa [W/m\cdotK]');
    case 'Kappa vs T RSquared'
        xarg = 'MeanBathTemperature';
        yarg = 'Kappa';
        RSquaredTolerance = 0.999; 
        figure; hold on; GeneralPlot2(MS,xarg,yarg,RSquaredTolerance);
        xlabel('T_{SD2} [K]'); ylabel('\kappa [W/m\cdotK]');
    case 'Kappa vs Probe T'
        xarg = 'MeanProbeTemperature';
        yarg = 'Kappa';
        figure; hold on; GeneralPlot1(MS,xarg,yarg);
        xlabel('T_{Probe} [K]'); ylabel('\kappa [W/m\cdotK]');

    otherwise
        disp('Invalid PlotOption in General_Kappa_Driver_v1!');
end

end


function [outstruct,ProcessParameters] = LoadKappaData(LoadOption)
switch LoadOption
    case 'MgCrO 110 HParallel'
        ProcessParameters.Length = 2e-3;
        ProcessParameters.Width = 3e-3;
        ProcessParameters.Thickness = 0.4e-3;
        ProcessParameters.HeaterCurve = load('E:\IanComputer\Documents\Physics\Minhyea Lee Research\Maglab_May2022\SCM2\MgCr2O4_110\MgCrO_110_HeaterFit_LogSpacePoly_O7.mat');
        ProcessParameters.HeaterCurve = ProcessParameters.HeaterCurve.HtrFit;
        Fileroot = 'E:\IanComputer\Documents\Physics\Minhyea Lee Research\Maglab_May2022\SCM2\MgCr2O4_110\Tswp\';
        ProcessParameters.Kappa_vs_T_Fields = [0,18];
    case 'MgCrO 110 HParallel NewLT'
        ProcessParameters.Length = 2e-3;
        ProcessParameters.Width = 3e-3;
        ProcessParameters.Thickness = 0.4e-3;
        ProcessParameters.HeaterCurve = load('E:\IanComputer\Documents\Physics\Minhyea Lee Research\Maglab_May2022\SCM2\MgCr2O4_110\MgCrO_110_HeaterFit_LogSpacePoly_O7.mat');
        ProcessParameters.HeaterCurve = ProcessParameters.HeaterCurve.HtrFit;
        Fileroot = 'E:\IanComputer\Documents\Physics\Minhyea Lee Research\Maglab_May2022\SCM2\MgCr2O4_110\Tswp_LTNewPoints\';
        ProcessParameters.Kappa_vs_T_Fields = [0,18];
    otherwise
        disp('Invalid LoadOption in LoadKappaData')
end
ProcessParameters.A_Over_L = ProcessParameters.Width.*ProcessParameters.Thickness./ProcessParameters.Length;
cd(Fileroot);

% Load all Data at different fixed fields into structre.
for i=1:length(ProcessParameters.Kappa_vs_T_Fields)
    CurrentFieldString = [num2str(ProcessParameters.Kappa_vs_T_Fields(i)),'T'];
    DStore = dir(['*',CurrentFieldString,'*.mat']);
    for j=1:length(DStore)
        outstruct(i,j) = load(DStore(j).name);
    end
end

% Append the file names. Do this after to preserve the likeness of the load
% structures.
for i=1:length(ProcessParameters.Kappa_vs_T_Fields)
    CurrentFieldString = [num2str(ProcessParameters.Kappa_vs_T_Fields(i)),'T'];
    DStore = dir(['*',CurrentFieldString,'*.mat']);
    for j=1:length(DStore)
        outstruct(i,j).FileName = DStore(j).name;
    end
end

end


function outstruct = ProcessKappaData(instruct,ProcessParameters,NewCalibrationCurves)
outstruct = instruct;
OutstructDimensions  = size(outstruct);
for i=1:OutstructDimensions(1)
    disp(['Now Processing Field Value ' ,num2str(ProcessParameters.Kappa_vs_T_Fields(i)),'T']);
    Updatebar = textprogressbar(OutstructDimensions(2));
    for j=1:OutstructDimensions(2)
        if length(NewCalibrationCurves)==2
            outstruct(i,j).HotTemp = NewCalibrationCurves{1}(outstruct(i,j).HotRes)';
            outstruct(i,j).ColdTemp= NewCalibrationCurves{2}(outstruct(i,j).ColdRes)';
        end
        outstruct(i,j).DeltaT  = outstruct(i,j).HotTemp - outstruct(i,j).ColdTemp;
        outstruct(i,j).MeanBathTemperature = mean(outstruct(i,j).BathTemp);
        outstruct(i,j).MeanSampleTemperature = mean(outstruct(i,j).HotTemp+outstruct(i,j).ColdTemp)./2;
        outstruct(i,j).MeanProbeTemperature = mean(outstruct(i,j).ProbeTemp); 
        outstruct(i,j).HeaterPower = (outstruct(i,j).HtrCurrent.^2).*...
            ProcessParameters.HeaterCurve(outstruct(i,j).MeanBathTemperature);
        
        HeaterCurrentSteps = length(outstruct(i,j).Time)./outstruct(i,j).SampPerStep;
        Tolerance = .4;
        SingleStepLogical = logical([zeros(1,outstruct(i,j).SampPerStep*(1-Tolerance)),...
            ones(1,outstruct(i,j).SampPerStep*(Tolerance))]);
%         Tolerance = .99;
%         SingleStepLogical = logical([zeros(1,floor(outstruct(i,j).SampPerStep*(1-Tolerance))),...
%             ones(1,floor(outstruct(i,j).SampPerStep*(Tolerance)))]);

        CutSteps = logical(repmat(SingleStepLogical,1,...
            HeaterCurrentSteps));
        
        ZeroCurrentInds = 1:outstruct(i,j).SampPerStep;
        outstruct(i,j).ZeroCurrentHotRes = mean(outstruct(i,j).HotRes(ZeroCurrentInds));
        outstruct(i,j).ZeroCurrentColdRes = mean(outstruct(i,j).ColdRes(ZeroCurrentInds));
        outstruct(i,j).ZeroCurrentBathTemp = mean(outstruct(i,j).BathTemp(ZeroCurrentInds));
        
        CutFields = {'Time','BathTemp','DeltaT','HeaterPower'};
        for q = 1:length(CutFields)
            outstruct(i,j).(['C_',CutFields{q}]) = outstruct(i,j).(CutFields{q})(CutSteps);
        end
        
        [FitResult,GOF] = polyfit(outstruct(i,j).C_HeaterPower,outstruct(i,j).C_DeltaT,1);
        outstruct(i,j).Kappa_RSquared = 1 - (GOF.normr/norm(outstruct(i,j).C_DeltaT - mean(outstruct(i,j).C_DeltaT)))^2;
        outstruct(i,j).Kappa = FitResult(1).*ProcessParameters.A_Over_L;
        outstruct(i,j).Kappa_Intercept = FitResult(2);
    end
    
end
end

function GeneralPlot1(instruct,xarg,yarg)
InstructDimensions = size(instruct);
if InstructDimensions(1)==1
    colorz = [1 0 0];
else
    colorz = pmkmp(InstructDimensions(1),'Swtth');
end
for i=1:InstructDimensions(1)
    plot([instruct(i,:).(xarg)],[instruct(i,:).(yarg)],'o','Color',...
        colorz(i,:),'MarkerFaceColor',brc(colorz(i,:),.2),'DisplayName',...
        [num2str(mean([instruct(i,:).Bfield])),'T']);
end

end


function GeneralPlot2(instruct,xarg,yarg,tolerance)
InstructDimensions = size(instruct);
if InstructDimensions(1)==1
    colorz = [1 0 0];
else
    colorz = pmkmp(InstructDimensions(1),'Swtth');
end
for i=1:InstructDimensions(1)
    XData = [instruct(i,:).(xarg)];
    YData = [instruct(i,:).(yarg)];
    RSquaredData = [instruct(i,:).Kappa_RSquared];
    KeepIndices = RSquaredData>=tolerance; 
    plot(XData(KeepIndices),YData(KeepIndices),'o','Color',...
        colorz(i,:),'MarkerFaceColor',brc(colorz(i,:),.2),'DisplayName',...
        [num2str(mean([instruct(i,:).Bfield])),'T']);
end

end