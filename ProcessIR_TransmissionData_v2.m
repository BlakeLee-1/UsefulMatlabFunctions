% ProcessIR_TransmissionData_v2
% Created 6/30/22, Edited 7/1/22
% Ian Leahy
%
% Input Arguments:
% LoadOption (String)
%   'Load 1' -- Data from first IR load, H||plane
%   'Load 2' -- Data from second IR load, H||c
%
% SampleOption (String) -- Different Loads have different samples
% available. The options are split into the different loads below. To see
% how the samples are aligned, please check the powerpoint 'Delafossite
% Sample Orientations'.
%   Load 1 Available:
%   'CsYbSe2'  'CsErSe2'  'CsErSe2 Detail'  'CsPrSe2'  'CsYbSe2 Sample 2'
%
%   Load 2 Available:
%   'CsYbSe2'  'CsErSe2'  'CsPrSe2'  'CsYbSe2 Sample 2'  'KYbSe2'
%
% PlotOption (String) -- Different options cause different plots. Use in
% combo with NormalizationOption  to plot different things.
%   'Transmission' -- Plots 2D plot of data.
%   'Transmission Vertical Offset' -- Plots 2D plot with vertical offset to
%   look like a waterfall.
%   'Transmission Colorplot' -- Plots a surface plot.
%   'Transmission Colorplot Slider' -- Plots colorplot with a slider
%   mapped to the minimum colorbar limit. Helpful to resolve small
%   features! Bounds of the colorbar (greater than one on the high end...)
%   can be set in this switch case option.
%   'Add Lines' -- plots the overlay lines in 3 dimensions on active plot.
%
% NormalizationOption (String) -- Select which normalization to be plotted.
%   'Raw' -- The as measured transmission data.
%   'MikeNorm' -- Mike Ozerov's normalization, given in the 'NORM'
%   datafiles.
%   'SingleMaxNorm' -- All data normalized by the single largest raw data
%   value at any field and wavenumber.
%   '0TNorm' -- All data normalized (division) to the 0T data.
%   'MaxFieldNorm' -- Normalization data set created from the maximum raw
%   data value at each wave number selected from data at all fields.
%   'BGSub' -- This normalization is only valid for Load2 samples. If an
%   the BG is not valid, the function will error. The normalization is
%   currently set to a straight division of the raw data by the pinhole
%   data, but that can be changed in the LoadIRData internal function.
%
% Example Uses:
%   ProcessIR_TransmissionData_v2('Load 1','CsYbSe2','Transmission','MikeNorm')
%   ProcessIR_TransmissionData_v2('Load 1','CsYbSe2','Transmission Vetical
%   Offset','MikeNorm');
%   ProcessIR_TransmissionData_v2('Load 1','CsYbSe2','Transmission
%   Colorplot','MikeNorm')
%   To overlay lines if available:
%   ProcessIR_TransmissionData_v2('Load 1','CsPrSe2','Transmission
%   Colorplot Slider', 'MikeNorm');
%   ProcessIR_TransmissionData_v2('Load 1','CsPrSe2','Add Lines','Mike
%   Norm');
function ProcessIR_TransmissionData_v2(LoadOption,SampleOption,PlotOption,NormalizationOption)

[WhereIsTheData] = SetFilePathHere(); 
[MS] = LoadIRData(LoadOption,SampleOption,WhereIsTheData);

switch PlotOption
    case 'Transmission'
        figure; hold on;
        xarg = 'Wavenumber';
        yarg = ['Transmission_',NormalizationOption];
        GenericIRPlot(MS,xarg,yarg);
        xlabel('Wavenumber [cm^{-1}]');
        ylabel(['Transmission ',NormalizationOption]);
        xlim([0 720]);
    case 'Transmission Vertical Offset'
        figure; hold on;
        xarg = 'Wavenumber';
        yarg = ['Transmission_',NormalizationOption];
        GenericIRPlot_VOffset(MS,xarg,yarg,0.1);
        xlabel('Wavenumber [cm^{-1}]');
        ylabel(['Offset Transmission ',NormalizationOption]);
        title(SampleOption);
        xlim([0 720]);
    case 'Transmission Colorplot'
        figure; hold on;
        zarg = ['Transmission_',NormalizationOption];
        GenericIRSurfacePlot(MS,zarg);
        title([SampleOption,' ',NormalizationOption]);
    case 'Transmission Colorplot Slider'
        figure; hold on;
        zarg = ['Transmission_',NormalizationOption];
        [ColorbarHandle] = GenericIRSurfacePlot(MS,zarg);
        title([SampleOption,' ',NormalizationOption]);
        CBPosition=get(ColorbarHandle,'position');
        
        % S = ['set(gca,''CLim'' get(gca,''CLim'')+['  num2str(0.02) ',0])' ];
        S = ['set(gca,''CLim'',[get(gcbo,''Value''),1])' ];
        uic=uicontrol('style','slider','units','normalized',...
            'position',[CBPosition(1)+0.075 CBPosition(2) 0.04 CBPosition(4)],...
            'min',0.5,'max',.999,'value',0.7,'callback',S);
    case 'Add Lines'
        LoadLines = LoadIRLines(LoadOption,SampleOption,WhereIsTheData);
        PlotIRLines(LoadLines);
    case 'INCOMPLETE PEAK FINDER'
        Tinds = imregionalmin(ZInterp,4);
        plot3(XInterp(Tinds),YInterp(Tinds),ZInterp(Tinds),'ok','MarkerFaceColor',brc([0 0 0],.5))
        Tinds2 = and(Tinds,ZInterp<=.975);
        plot3(XInterp(Tinds2),YInterp(Tinds2),ZInterp(Tinds2),'or','MarkerFaceColor',brc([0 0 0],.5))
        Tinds2 = and(Tinds,ZInterp<=.90);
        
        T1 = clusterdata(X,3);
        
        SizeXI = size(XInterp);
        xv = reshape(XInterp,SizeXI(1).*SizeXI(2),1);
        
        SizeYI = size(YInterp);
        yv = reshape(YInterp,SizeYI(1).*SizeYI(2),1);
        
        SizeZI = size(ZInterp);
        zv = reshape(ZInterp,SizeZI(1).*SizeZI(2),1);
        
        tv = reshape(Tinds2,SizeZI(1).*SizeZI(2),1);
        
        T1 = clusterdata([xv(tv),yv(tv),zv(tv)],16);
    otherwise
        error('Invalid PlotOption! Try one that exists or create a new one!');
end
end


% Internal functions below here.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Where is the data... change this!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% This should be the only thing you need to change.%%%%%%%%%%%%%%
%%%%%%%%% This folder should be where all of the data folders are.%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [WhereIsTheData] = SetFilePathHere()
WhereIsTheData = '/home/blake/Documents/MATLAB/Processing Code/';
end

% Loads IR Data from the June 2022 Maglab.
function [MasterStructure] = LoadIRData(LoadOption,SampleOption,WhereIsTheData)

try
    cd(WhereIsTheData);
catch
    error('Unable to change path -- did you change where the data is?');
end


switch LoadOption
    case 'Load 1'
        LoadFileroot = 'Load1_TrimData';
    case 'Load 2'
        LoadFileroot = 'Load2_TrimData';
    otherwise
        error('Invalid LoadOption in LoadIRData! Make sure LoadOption is either Load 1 or Load 2!');
end
cd(LoadFileroot);
switch SampleOption
    case 'CsYbSe2'
        SearchString = 'CsYb_100_';
    case 'CsErSe2'
        SearchString = 'CsEr_100_';
    case 'CsErSe2 Detail'
        SearchString = 'CsEr_100-FIR';
    case 'CsPrSe2'
        SearchString = 'CsPr_100_';
    case 'CsYbSe2 Sample 2'
        SearchString = 'CsYb_110_';
    case 'KYbSe2'
        SearchString = 'KYb_100_';
    otherwise
        error('Invalid SampleOption in LoadIRData! Did you pick a valid sample?');
end


DStore_Raw = dir(['*',[SearchString,'RAW'],'*.dat']);
DStore_Norm = dir(['*',[SearchString,'Norm'],'*.dat']);

if strcmp(LoadOption,'Load 2')
    % Pinhole File. Point to it here.
    try
        BGMat = readmatrix([WhereIsTheData,'Pin_hole_FrTr_6mu_20Khz_3min_17-20_.0.txt']);
        BGVec = BGMat(:,2);
    catch
        disp('Background not loaded successfully. Background set to unity. Change the filepath for loading LoadIRData in the internal function LoadIRData!');
        BGMat = 1;
        BGVec = 1;
    end
else
    BGMat = 1;
    BGVec = 1;
end
% Load MasterStructure of all data

try
    RawData = readmatrix(DStore_Raw(1).name);
    NormData = readmatrix(DStore_Norm(1).name);
catch
    error(sprintf(['Unable to read datafile in LoadIRData. Either the Fileroot\n',...
        'location of the data must be updated or the selected LoadOption\n',...
        'and SampleOption are incompatible. Make sure the selected sample\n',...
        'was in the selected load (check the comments if you forget)!']));
end
FieldValues_Raw = RawData(1,2:end-1);
FieldValues_Norm = NormData(1,2:end-1);
JustTransmissionData_Raw = RawData(2:end,2:end-1);
JustTransmissionData_Norm = NormData(2:end,2:end-1);
ActualMax = max(JustTransmissionData_Raw,[],2);
NormValue = max(max(JustTransmissionData_Raw));
Norm0T = JustTransmissionData_Raw(:,1);
for i=1:length(FieldValues_Raw)
    MasterStructure(i).Wavenumber = RawData(2:end,1);
    MasterStructure(i).Field = FieldValues_Raw(i);
    MasterStructure(i).Transmission_Raw = RawData(2:end,i+1);
    % Normalization schemes.
    % 'Transmission_MikeNorm' is data taken from the normalized data file
    % sent by Mike.
    % 'Transmission_SingleMaxNorm' is the raw transmission data divided by
    % the single largest transmission value at all measured fields on a
    % given sample.
    % 'Transmission_0TNorm' is the raw transmission data normalized to the
    % 0T data.
    % 'Transmission_MaxFieldNorm' is similar to Mike Ozerov's scheme. This
    % is taking a row of the JustTransmissionData_Raw matrix and finding
    % the max across a single row, corresponding to the max transmission at
    % a given field. This becomes the normalization at that wavenumber.
    % Mike's normalization scheme realistically invovles more detail than
    % this, but this works fairly well for something quick.
    % 'Transmission_BGSub' is supposed to remove the BG pinhole
    % measurement, however this is only valid for Load 2. Otherwise the
    % data is just the transmission data.
    MasterStructure(i).Transmission_MikeNorm = NormData(2:end,i+1);
    MasterStructure(i).Transmission_SingleMaxNorm = MasterStructure(i).Transmission_Raw./NormValue;
    MasterStructure(i).Transmission_0TNorm = MasterStructure(i).Transmission_Raw./Norm0T;
    MasterStructure(i).ML_Average = RawData(2:end,end);
    MasterStructure(i).ML_Average_Correct = ActualMax;
    MasterStructure(i).Transmission_MaxFieldNorm = MasterStructure(i).Transmission_Raw./ActualMax;
    try
        MasterStructure(i).Transmission_BGSub = MasterStructure(i).Transmission_Raw./BGVec;
    catch
        BGInterp = Interp1NonUnique(BGMat(:,1),BGMat(:,2),MasterStructure(i).Wavenumber);
        MasterStructure(i).Transmission_BGSub = MasterStructure(i).Transmission_Raw./BGInterp;
    end
end

end

% GenericIRPlot Function, plots xarg vs yarg in 2D.
function GenericIRPlot(MS,xarg,yarg)
for i=1:length(MS)
    plot(MS(i).(xarg),MS(i).(yarg),'-',...
        'Color',ColormapInterpolate_InBounds([MS(i).Field],[0,18]),...
        'DisplayName',[num2str(MS(i).Field),' T']);
end
end

% GenericIRPlot Function with a vertical offset, plots xarg vs yarg in 2D
% with a vertical offset.
function GenericIRPlot_VOffset(MS,xarg,yarg,VerticalOffset)
for i=1:length(MS)
    plot(MS(i).(xarg),MS(i).(yarg)+VerticalOffset.*(i-1),'-',...
        'Color',ColormapInterpolate_InBounds([MS(i).Field],[0,18]),...
        'DisplayName',[num2str(MS(i).Field),' T']);
end
end

% GenericIRPlot Surface plot function. plots vs. wavenumber (x) and field
% (y).
function [ColorbarHandle] = GenericIRSurfacePlot(MS,zarg)
surf(MS(1).Wavenumber,[MS(:).Field],[MS(:).(zarg)]',...
    'EdgeColor','none','FaceColor','interp');
temp = gca;
temp.CLim = [.5 1];
colormap(flipud(jet(256)));
ylim([0 17.5]);
xlim([0 720]);
ColorbarHandle = colorbar;
xlabel('Wavenumber [cm^{-1}]');
ylabel('\mu_0H [T]');

end

% Function to load IR Lines if available.
function [LoadLines] = LoadIRLines(LoadOption,SampleOption,WhereIsTheData)

switch LoadOption
    case 'Load 1'
        LoadStr = 'Load1_Lines';
    case 'Load 2'
        LoadStr = 'Load2_Lines';
    otherwise
        error('Invalid LoadOption in LoadIRLines! Make sure you input either Load 1 or Load 2!');
end
cd(WhereIsTheData);
cd(LoadStr); 

switch SampleOption
    case 'CsYbSe2'
        SearchString = 'CsYbSe2_S1';
    case 'CsErSe2'
        SearchString = 'CsErSe2';
    case 'CsErSe2 Detail'
        SearchString = 'CsErSe2';
    case 'CsPrSe2'
        SearchString = 'CsPrSe2';
    case 'CsYbSe2 Sample 2'
        SearchString = 'CsYbSe2_S2';
    case 'KYbSe2'
        SearchString = 'KYbSe2';
    otherwise
        disp('Invalid SampleOption in LoadIRData!');
end

DStore = dir(['*',SearchString,'*.mat']);
if isempty(DStore)
    error('File not found. Maybe Ian has not fit the lines yet? Maybe the file is not in the right place or load. Give it a check!');
end
temp = load(DStore(1).name);
FNames = fieldnames(temp);
LoadLines = temp.(FNames{1});
end

function  PlotIRLines(LoadLines)
hold on;
for i=1:length(LoadLines)
    CM = LoadLines{i};
    plot3(CM(:,1),CM(:,2),CM(:,3)+.01,'-r','Linewidth',1.5);
end
end
