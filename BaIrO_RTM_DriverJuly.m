% Ian Leahy, Peter Siegfried, Chris Pocs
% July 19, 2019
% File for processing RTM data

function [SortedMC]=BaIrO_RTM_DriverJuly(LoadOpt,PlotOpt,BPower,FitAngleOpt,orderedplot,derivopt)

inputlist = [];
MC = LoadRTM(LoadOpt,inputlist);

figure; hold on;
colorz = distinguishable_colors(length(MC));
smopt = 20;
for i=1:length(MC)
    switch PlotOpt
        case 'B All'
            
            hold on; plot(MC{i}.Field.^BPower,MC{i}.Frequency,'.','Color',...
                brc(colorz(i,:),.5),'DisplayName',MC{i}.TLeg,'MarkerSize',8)
            hold on; plot(MC{i}.Field.^BPower,smooth(MC{i}.Frequency,smopt),'-','Color',...
                colorz(i,:),'DisplayName',[MC{i}.TLeg,' Smooth'])
            
        case 'B Up'
            lim = -.0001;
            keepinds = diff(MC{i}.Field)>lim;
            hold on; plot(MC{i}.Field(keepinds).^BPower,MC{i}.Frequency(keepinds),'.','Color',...
                brc(colorz(i,:),.5),'DisplayName',MC{i}.TLeg,'MarkerSize',8)
            hold on; plot(MC{i}.Field(keepinds).^BPower,smooth(MC{i}.Frequency(keepinds),smopt),'.','Color',...
                colorz(i,:),'DisplayName',[MC{i}.TLeg,' Smooth'])
           
        case 'B Up Temp'
            lim = -.0001;
            keepinds = diff(MC{i}.Field)>lim;
            hold on; plot(MC{i}.Field(keepinds).^BPower,MC{i}.Frequency(keepinds)-46480,'.','Color',...
                brc(colorz(i,:),.5),'DisplayName',MC{i}.TLeg,'MarkerSize',8)
            hold on; plot(MC{i}.Field(keepinds).^BPower,smooth(MC{i}.Frequency(keepinds)-46480,smopt),'.','Color',...
                colorz(i,:),'DisplayName',[MC{i}.TLeg,' Smooth'])
            
        case 'B Down'
            lim = -.0001;
            keepinds = diff(MC{i}.Field)<lim;
            hold on; plot(MC{i}.Field(keepinds).^BPower,MC{i}.Frequency(keepinds),'.','Color',...
                brc(colorz(i,:),.5),'DisplayName',MC{i}.TLeg,'MarkerSize',8)
            hold on; plot(MC{i}.Field(keepinds).^BPower,smooth(MC{i}.Frequency(keepinds),smopt),'-','Color',...
                colorz(i,:),'DisplayName',[MC{i}.TLeg,' Smooth'])
            
        case 'B Index'
            
            hold on; plot(MC{i}.Frequency,'.','Color',...
                brc(colorz(i,:),.5),'DisplayName',MC{i}.TLeg,'MarkerSize',8)
            hold on; plot(smooth(MC{i}.Frequency,smopt),'-','Color',...
                colorz(i,:),'DisplayName',[MC{i}.TLeg,' Smooth'])
            
    end
end
% ylim([4.544e4,4.546e4]);
legend show;
Label_Plot(['B^',num2str(BPower)],['T^',num2str(BPower)],'f','Hz');

switch FitAngleOpt
    case 'Fit'
        FitAngleDep(MC);
    case 'Fit 20K'
        FitAngleDep_20K(MC);
        
    case 'Fit S3 4K'
        FitAngleDep_4K(MC);
        
    case 'Fit S3 4K Directional'
        FitAngleDep_4K_Directional(MC,'Down')
        
    case 'Fit S3 FFS'
        FitAngleDep_4K_Directional_FFS(MC,'Down','Sample 3 All');
    case 'Fit S4 FFS'
        FitAngleDep_4K_Directional_FFS(MC,'Down','Sample 4 All');
    case 'Fit S4 Angle FFS'
        FitAngleDep_4K_Directional_FFS(MC,'Down','Sample 4 Angles');
    case 'No'
        
end
angles = LU_RTM(LoadOpt);
SortedMC = OrderData(MC,angles);
if strcmp('Ordered',orderedplot)
    figure;
    shift = 0;
    colorz = varycolor(length(SortedMC));
    for i=1:length(SortedMC)
        lim = -.0001;
        switch PlotOpt
            case 'B Up'
                keepinds = diff(SortedMC{i}.Field)>lim;
            case 'B Down'
                keepinds = diff(SortedMC{i}.Field)<lim;
            case 'B All'
                keepinds = ones(1,length(SortedMC{i}.Field));
            otherwise
                disp('Invalid field exclusion')
        end
        hold on; plot(SortedMC{i}.Field(keepinds).^BPower,SortedMC{i}.Frequency(keepinds)+shift,'.','Color',...
            colorz(i,:),'DisplayName',SortedMC{i}.TLeg,'MarkerSize',8)
        temphold = SortedMC{i}.Frequency(keepinds)+shift;
        maxvals(i) = temphold(1);
        shift = shift - 75;
        labellist{i} = [num2str(angles(i)),char(176)];
    end
    Label_Plot(['B^',num2str(BPower)],['T^',num2str(BPower)],'f','Hz');
    %     labelpoints(repmat(65,1,length(SortedMC)),maxvals',labellist','Color',colorz,'FontSize',20);
end

switch derivopt
    case 'Spline'
        SortedMC = FindDeriv(SortedMC,'Spline');
        angles = LU_RTM('Sample 4 All');
        figure;
        shift = 0;
        colorz = varycolor(length(SortedMC));
        for i=1:length(SortedMC)
            hold on; plot(SortedMC{i}.FieldDeriv,SortedMC{i}.FreqDeriv,'.','Color',...
                colorz(i,:),'DisplayName',SortedMC{i}.TLeg,'MarkerSize',8)
            shift = shift - 10;
        end
        Label_Plot(['B'],['T'],'df/dH','Hz/T');
        
        
    case 'None'
        
end

% PlotSurface(SortedMC);
PlotFieldCuts(SortedMC,PlotOpt,LoadOpt);
% PlotFieldCuts_Repeat(SortedMC,PlotOpt,LoadOpt);
end


function outcell = LoadRTM(LoadOpt,inputlist)
strfields = {'Field','Frequency','Amplitude'};

switch LoadOpt
    case 'Phi1 All'
        fr = 'E:\IanComputer\Documents\Physics\Minhyea Lee Research\LosAlamos_July2019\BaIrO_Phi1_FFH_4K\';
        cd(fr);
        DStore = dir('*.dat');
        [Angles,Temperatures] = LU_RTM('Phi1 All');
    case 'Phi1 2K'
        fr = 'E:\IanComputer\Documents\Physics\Minhyea Lee Research\LosAlamos_July2019\BaIrO_Phi1_FFH_1p5K\';
        cd(fr);
        DStore = dir('*.dat');
        [Angles,Temperatures] = LU_RTM('Phi1 2K');

    case 'List'
        fr = 'E:\IanComputer\Documents\Physics\Minhyea Lee Research\LosAlamos_March2019\031819\031819\';
        cd(fr);
        for i=1:length(inputlist)
            DStore(i,1) = dir(inputlist{i});
        end
    case 'Select'
        fr = 'E:\IanComputer\Documents\Physics\Minhyea Lee Research\LosAlamos_March2019\Sample1_March18_FFHits\';
        filename = uigetfile('.dat');
        DStore = dir(filename);
    otherwise
        disp('Invalid LoadOpt! Select a valid Load Option!');
end
for i=1:length(DStore)
    temp = tdfread([fr,DStore(i).name]);
    fnames = fieldnames(temp);
    for j=1:length(fnames)
        outcell{i}.(strfields{j}) = temp.(fnames{j});
    end
    outcell{i}.date = DStore(i).date;
    outcell{i}.name = DStore(i).name;
    try
        outcell{i}.IDNum = str2num(DStore(i).name(2:4));
        outcell{i}.Temp = Temperatures(FileInd==outcell{i}.IDNum);
    catch
        outcell{i}.IDNum = 9999;
        outcell{i}.Temp = 0;
    end
    try
        outcell{i}.IDNum = str2num(DStore(i).name(2:4));
        %         outcell{i}.Temp = Angles(FileInd==outcell{i}.IDNum);
        outcell{i}.Temp = Angles(i);
        
    catch
        outcell{i}.IDNum = 9999;
        outcell{i}.Temp = 0;
        
    end
    
    %         outcell{i}.TLeg = [num2str(outcell{i}.Temp), ' K'];
    outcell{i}.TLeg = [num2str(outcell{i}.Temp), ' deg. df = ',outcell{i}.name];
    
end

end

function FitAngleDep(MC)

for i=1:length(MC)
    temp_pfit = polyfit(MC{i}.Field.^2,MC{i}.Frequency,1);
    slope(i) = temp_pfit(1);
    yint(i) = temp_pfit(2);
end
% angles = [170,180,0,20,30,40,50,60,74,85,100,110,120,130,140,150,160];
angles = LU_RTM('Sample 2');
figure; hold on;
plot(angles,slope,'or','MarkerSize',10,'MarkerFaceColor',brc([1,0,0],.5))
cosfunc = @(a,shift,ang) a.*cosd(2*(ang)+shift);
tempangs = -90:1:200;
hold on; plot(tempangs,cosfunc(.53,75,tempangs),'--k');
% hold on; plot(tempangs,cosfunc(.53,55,tempangs),'--k');
Label_Plot('Angle','deg','Slope','');
tangs = [0,10,20,30,40,50,60,74,85,100,110,120,130,140,150,160,170,180,190];
fangs = [33.24,33.68,36.26,39.07,44.51,100,50.7,47,42.61,32,27.71,26,26.55,26.9,44.2,39.22,36.96,35.06,35.5];
fangs2 =[33.24,33.68,36.26,39.07,48.82,100,50.8,46.9,44.75,100,90,80,70,58.55,44.2,39.22,36.96,35.06,35.5];
yyaxis right; plot(tangs,fangs2,'ob','MarkerSize',10,'MarkerFaceColor',brc([0,0,1],.5));
yyaxis right; plot(tangs,fangs2,'--k');
ylim([25,65]);

end

function FitAngleDep_20K(MC)

for i=1:length(MC)
    clipinds = MC{i}.Field>=7;
    xfitval = MC{i}.Field(clipinds).^2;
    yfitval = MC{i}.Frequency(clipinds);
    temp_pfit = polyfit(xfitval,yfitval,1);
    slope(i) = temp_pfit(1);
    yint(i) = temp_pfit(2);
end
% angles = [170,180,0,20,30,40,50,60,74,85,100,110,120,130,140,150,160];
angles = LU_RTM('Sample 2 20K');
figure; hold on;
plot(angles,slope,'og','MarkerSize',10,'MarkerFaceColor',brc([0,1,0],.5))
cosfunc = @(a,shift,ang) a.*cosd(2*(ang)+shift);
tempangs = -90:1:200;
hold on; plot(tempangs,cosfunc(.53,75,tempangs),'--k');
% hold on; plot(tempangs,cosfunc(.53,55,tempangs),'--k');
Label_Plot('Angle','deg','Slope','');
% tangs = [0,10,20,30,40,50,60,74,85,100,110,120,130,140,150,160,170,180,190];
% fangs = [33.24,33.68,36.26,39.07,44.51,100,50.7,47,42.61,32,27.71,26,26.55,26.9,44.2,39.22,36.96,35.06,35.5];
% fangs2 =[33.24,33.68,36.26,39.07,48.82,100,50.8,46.9,44.75,100,90,80,70,58.55,44.2,39.22,36.96,35.06,35.5];
% yyaxis right; plot(tangs,fangs2,'ob','MarkerSize',10,'MarkerFaceColor',brc([0,0,1],.5));
% yyaxis right; plot(tangs,fangs2,'--k');
% ylim([25,65]);

end

function FitAngleDep_4K(MC)

for i=1:length(MC)
    clipinds = MC{i}.Field>=7;
    xfitval = MC{i}.Field(clipinds).^2;
    yfitval = MC{i}.Frequency(clipinds);
    temp_pfit = polyfit(xfitval,yfitval,1);
    slope(i) = temp_pfit(1);
    yint(i) = temp_pfit(2);
end
% angles = [170,180,0,20,30,40,50,60,74,85,100,110,120,130,140,150,160];
angles = LU_RTM('Sample 3 4K');
figure; hold on;
incolor = [0,0,0];
plot(angles,5.*slope,'o','Color',incolor,'MarkerSize',10,'MarkerFaceColor',brc(incolor,.5))
cosfunc = @(a,shift,ang) a.*cosd(2*(ang)+shift);
tempangs = -90:1:200;
% hold on; plot(tempangs,cosfunc(.53,75,tempangs),'--k');
% hold on; plot(tempangs,cosfunc(.53,55,tempangs),'--k');
Label_Plot('Angle','deg','Slope','');
% tangs = [0,10,20,30,40,50,60,74,85,100,110,120,130,140,150,160,170,180,190];
% fangs = [33.24,33.68,36.26,39.07,44.51,100,50.7,47,42.61,32,27.71,26,26.55,26.9,44.2,39.22,36.96,35.06,35.5];
% fangs2 =[33.24,33.68,36.26,39.07,48.82,100,50.8,46.9,44.75,100,90,80,70,58.55,44.2,39.22,36.96,35.06,35.5];
% yyaxis right; plot(tangs,fangs2,'ob','MarkerSize',10,'MarkerFaceColor',brc([0,0,1],.5));
% yyaxis right; plot(tangs,fangs2,'--k');
% ylim([25,65]);

end

function FitAngleDep_4K_Directional(MC,diropt)

lim = -.0001;
for i=1:length(MC)
    switch diropt
        case 'Up'
            keepinds = diff(MC{i}.Field)>lim;
        case 'Down'
            keepinds = diff(MC{i}.Field)<lim;
        otherwise
            disp('Invalid field exclusion')
    end
    xin = MC{i}.Field(keepinds);
    yin = MC{i}.Frequency(keepinds);
    clipinds = xin>=5;
    xfitval = xin(clipinds).^2;
    yfitval = yin(clipinds);
    temp_pfit = polyfit(xfitval,yfitval,1);
    slope(i) = temp_pfit(1);
    yint(i) = temp_pfit(2);
end
% angles = [170,180,0,20,30,40,50,60,74,85,100,110,120,130,140,150,160];
angles = LU_RTM('Sample 3 4K');
figure; hold on;
incolor = [0,0,0];
[angles,sortind]=sort(angles);
plot(angles,5.*slope(sortind),'o','Color',incolor,'MarkerSize',10,'MarkerFaceColor',brc(incolor,.5))
cosfunc = @(a,shift,ang) a.*cosd(2*(ang)+shift);
tempangs = -90:1:200;
% hold on; plot(tempangs,cosfunc(.53,75,tempangs),'--k');
% hold on; plot(tempangs,cosfunc(.53,55,tempangs),'--k');
Label_Plot('Angle','deg','Slope','');
% tangs = [0,10,20,30,40,50,60,74,85,100,110,120,130,140,150,160,170,180,190];
% fangs = [33.24,33.68,36.26,39.07,44.51,100,50.7,47,42.61,32,27.71,26,26.55,26.9,44.2,39.22,36.96,35.06,35.5];
% fangs2 =[33.24,33.68,36.26,39.07,48.82,100,50.8,46.9,44.75,100,90,80,70,58.55,44.2,39.22,36.96,35.06,35.5];
% yyaxis right; plot(tangs,fangs2,'ob','MarkerSize',10,'MarkerFaceColor',brc([0,0,1],.5));
% yyaxis right; plot(tangs,fangs2,'--k');
% ylim([25,65]);

end

function FitAngleDep_4K_Directional_FFS(MC,diropt,sampleopt)

lim = -.0001;
for i=1:length(MC)
    switch diropt
        case 'Up'
            keepinds = diff(MC{i}.Field)>lim;
        case 'Down'
            keepinds = diff(MC{i}.Field)<lim;
        otherwise
            disp('Invalid field exclusion')
    end
    xin = MC{i}.Field(keepinds);
    yin = MC{i}.Frequency(keepinds);
    clipinds = and(xin>=15,xin<=30);
    %     clipinds = and(xin>=15,xin<=20);
    %     clipinds = and(xin>=57,xin<=60);
    xfitval = xin(clipinds).^2;
    yfitval = yin(clipinds);
    temp_pfit = polyfit(xfitval,yfitval,1);
    slope(i) = temp_pfit(1);
    yint(i) = temp_pfit(2);
end
switch sampleopt
    case 'Sample 3 All'
        angles = LU_RTM('Sample 3 All');
    case 'Sample 4 All'
        angles = LU_RTM('Sample 4 All');
    case 'Sample 4 Angles'
        angles = LU_RTM('Sample 4 Angles');
end
figure; hold on;
incolor = [0,0,0];
[angles,sortind] = sort(angles);
plot(angles,5.*slope(sortind),'o','Color',incolor,...
    'MarkerSize',10,'MarkerFaceColor',brc(incolor,.5))
cosfunc = @(a,shift,ang) a.*cosd(2*(ang)+shift);
tempangs = -90:1:200;
hold on; plot(tempangs,cosfunc(.44,60,tempangs),'--k');
% hold on; plot(tempangs,cosfunc(.53,55,tempangs),'--k');
Label_Plot('Angle','deg','Slope','');
% tangs = [0,10,20,30,40,50,60,74,85,100,110,120,130,140,150,160,170,180,190];
% fangs = [33.24,33.68,36.26,39.07,44.51,100,50.7,47,42.61,32,27.71,26,26.55,26.9,44.2,39.22,36.96,35.06,35.5];
% fangs2 =[33.24,33.68,36.26,39.07,48.82,100,50.8,46.9,44.75,100,90,80,70,58.55,44.2,39.22,36.96,35.06,35.5];
% yyaxis right; plot(tangs,fangs2,'ob','MarkerSize',10,'MarkerFaceColor',brc([0,0,1],.5));
% yyaxis right; plot(tangs,fangs2,'--k');
% ylim([25,65]);

end

function FitAngleDep_Lines(MC,diropt,sampleopt)

lim = -.0001;
for i=1:length(MC)
    switch diropt
        case 'Up'
            keepinds = diff(MC{i}.Field)>lim;
        case 'Down'
            keepinds = diff(MC{i}.Field)<lim;
        otherwise
            disp('Invalid field exclusion')
    end
    xin = MC{i}.Field(keepinds);
    yin = MC{i}.Frequency(keepinds);
    clipinds = and(xin>=15,xin<=30);
    %     clipinds = and(xin>=15,xin<=20);
    %     clipinds = and(xin>=57,xin<=60);
    xfitval = xin(clipinds).^2;
    yfitval = yin(clipinds);
    temp_pfit = polyfit(xfitval,yfitval,1);
    slope(i) = temp_pfit(1);
    yint(i) = temp_pfit(2);
end
switch sampleopt
    case 'Sample 3 All'
        angles = LU_RTM('Sample 3 All');
    case 'Sample 4 All'
        angles = LU_RTM('Sample 4 All');
    case 'Sample 4 Angles'
        angles = LU_RTM('Sample 4 Angles');
end
figure; hold on;
% incolor = [0,0,0];
[angles,sortind] = sort(angles);
plot(angles,5.*slope(sortind),'o','Color',incolor,...
    'MarkerSize',10,'MarkerFaceColor',brc(incolor,.5))
cosfunc = @(a,shift,ang) a.*cosd(2*(ang)+shift);
tempangs = -90:1:200;
hold on; plot(tempangs,cosfunc(.44,60,tempangs),'--k');
% hold on; plot(tempangs,cosfunc(.53,55,tempangs),'--k');
Label_Plot('Angle','deg','Slope','');
% tangs = [0,10,20,30,40,50,60,74,85,100,110,120,130,140,150,160,170,180,190];
% fangs = [33.24,33.68,36.26,39.07,44.51,100,50.7,47,42.61,32,27.71,26,26.55,26.9,44.2,39.22,36.96,35.06,35.5];
% fangs2 =[33.24,33.68,36.26,39.07,48.82,100,50.8,46.9,44.75,100,90,80,70,58.55,44.2,39.22,36.96,35.06,35.5];
% yyaxis right; plot(tangs,fangs2,'ob','MarkerSize',10,'MarkerFaceColor',brc([0,0,1],.5));
% yyaxis right; plot(tangs,fangs2,'--k');
% ylim([25,65]);

end

function PlotFieldCuts(MC,diropt,sampleopt)

lim = -.0001;
fieldvals = [0:5:60];
meanwidth = 2;
for j=1:length(fieldvals)
    for i=1:length(MC)
        switch diropt
            case 'B Up'
                keepinds = diff(MC{i}.Field)>lim;
            case 'B Down'
                keepinds = diff(MC{i}.Field)<lim;
            case 'B All'
                keepinds = ones(1,length(MC{i}.Field));
            otherwise
                disp('Invalid field exclusion')
        end
        xin = MC{i}.Field(keepinds);
        yin = MC{i}.Frequency(keepinds);
        clipinds = and(xin>=(fieldvals(j)-meanwidth),xin<=(fieldvals(j)+meanwidth));
        ang.(['freq_',num2str(fieldvals(j))])(i) = trimmean(yin(clipinds),20);
        ang.(['Field_',num2str(fieldvals(j))])(i) = mean(xin(clipinds));
        
    end
end
switch sampleopt
    case 'Phi1 All'
        [angles,~] = LU_RTM('Phi1 All');
    case 'Phi1 2K'
        [angles,~] = LU_RTM('Phi1 2K');
end
figure; hold on;
[angles,sortind] = sort(angles);
colorz = hsv(length(fieldvals));
for i=1:length(fieldvals)
    hold on;
    %     plot(angles,ang.(['freq_',num2str(fieldvals(i))]),'-o','DisplayName',[num2str(fieldvals(i)),' T'],...
    %         'Color',colorz(i,:),'MarkerFaceColor',brc(colorz(i,:),.8))
    plot(angles,ang.(['freq_',num2str(fieldvals(i))])-ang.freq_0,'-o','DisplayName',[num2str(fieldvals(i)),' T'],...
        'Color',colorz(i,:),'MarkerFaceColor',brc(colorz(i,:),.8))
    
end
Label_Plot('Angle','deg','f','Hz');
title(diropt);
end


function PlotFieldCuts_Repeat(MC,diropt,sampleopt)

lim = -.0001;
fieldvals = [0:5:60];
meanwidth = 2;
for j=1:length(fieldvals)
    for i=1:length(MC)
        switch diropt
            case 'B Up'
                keepinds = diff(MC{i}.Field)>lim;
            case 'B Down'
                keepinds = diff(MC{i}.Field)<lim;
            case 'B All'
                keepinds = ones(1,length(MC{i}.Field));
            otherwise
                disp('Invalid field exclusion')
        end
        xin = MC{i}.Field(keepinds);
        yin = MC{i}.Frequency(keepinds);
        clipinds = and(xin>=(fieldvals(j)-meanwidth),xin<=(fieldvals(j)+meanwidth));
        ang.(['freq_',num2str(fieldvals(j))])(i) = trimmean(yin(clipinds),20);
        ang.(['Field_',num2str(fieldvals(j))])(i) = mean(xin(clipinds));
        
    end
end
switch sampleopt
    case 'Phi4 All'
        [angles,~] = LU_RTM('Phi4 All');
    case 'Phi4 40K'
        [angles,~] = LU_RTM('Phi4 40K');
end
figure; hold on;
[angles,sortind] = sort(angles);
colorz = hsv(length(fieldvals));
for i=1:length(fieldvals)
    hold on;
    plot(angles-180,ang.(['freq_',num2str(fieldvals(i))]),'-o','DisplayName',[num2str(fieldvals(i)),' T'],...
        'Color',colorz(i,:),'MarkerFaceColor',brc(colorz(i,:),.8));
    plot(angles,ang.(['freq_',num2str(fieldvals(i))]),'-o','DisplayName',[num2str(fieldvals(i)),' T'],...
        'Color',colorz(i,:),'MarkerFaceColor',brc(colorz(i,:),.8));
    plot(angles+180,ang.(['freq_',num2str(fieldvals(i))]),'-o','DisplayName',[num2str(fieldvals(i)),' T'],...
        'Color',colorz(i,:),'MarkerFaceColor',brc(colorz(i,:),.8));
    
end
Label_Plot('Angle','deg','f','Hz');
title(diropt);
end


function [angleout,Temperatures] = LU_RTM(option)

angleout = [];
switch option
    case 'Phi1 All'
        fr = 'E:\IanComputer\Documents\Physics\Minhyea Lee Research\LosAlamos_July2019\BaIrO_Phi1_FFH_4K\';
        cd(fr);
        DStore = dir('*.dat');
        for i=1:length(DStore)
            rundate = DStore(i,1).name(8:9);
            numdate = str2num(rundate);
            runID = str2num(DStore(i,1).name(2:4));
            if numdate==25
                Angles =  [180,150];
                FileInd = [94,96];
                Temps = repmat(3.9,1,length(FileInd));
                Temperatures(i) = Temps(FileInd==runID);
                angleout(i) = Angles(FileInd==runID);
            elseif numdate==26
                Angles =  [120,90,60,30];
                FileInd = [8,11,13,15];
                Temps = repmat(3.9,1,length(FileInd));
                Temperatures(i) = Temps(FileInd==runID);
                angleout(i) = Angles(FileInd==runID);
            end
            
            
        end
    case 'Phi1 2K'
        fr = 'E:\IanComputer\Documents\Physics\Minhyea Lee Research\LosAlamos_July2019\BaIrO_Phi1_FFH_1p5K\';
        cd(fr);
        DStore = dir('*.dat');
        for i=1:length(DStore)
            rundate = DStore(i,1).name(8:9);
            numdate = str2num(rundate);
            runID = str2num(DStore(i,1).name(2:4));
            if numdate==26
                Angles =  [180,150,120,90,75,60,45,30,15,135];
                FileInd = [18,23,25,27,29,31,33,35,37,39];
                Temps = repmat(1.5,1,length(FileInd));
                Temperatures(i) = Temps(FileInd==runID);
                angleout(i) = Angles(FileInd==runID);
            end
        end
    otherwise
        disp('Invalid option.')
end

end

function outdata=OrderData(MC,angles)

[sortedangles,sortinds] = sort(angles);

for i=1:length(MC)
    j = sortinds(i);
    outdata{i} = MC{j};
    outdata{i}.Angle = sortedangles(i);
end


end

function MC = FindDeriv(MC,opt)

for i=1:length(MC)
    switch opt
        case 'Spline'
            lim = -.0001;
            keepinds = diff(MC{i}.Field)<lim;
            xtemp = MC{i}.Field(keepinds);
            ytemp = MC{i}.Frequency(keepinds);
            [xData, yData] = prepareCurveData( xtemp, ytemp );
            % Set up fittype and options.
            ft = fittype( 'smoothingspline' );
            opts = fitoptions( 'Method', 'SmoothingSpline' );
            %         opts.SmoothingParam = 0.25894802118215;
            opts.SmoothingParam = 0.026727040805093052;
            % Fit model to data.
            [MC{i}.fitresult, MC{i}.gof] = fit( xData, yData, ft, opts );
            xinterp = 0:.01:65;
            MC{i}.FieldDeriv = xinterp(1:end-1);
            MC{i}.FreqDeriv = diff(MC{i}.fitresult(xinterp))./.01;
            
        otherwise
            disp('Invalid Option Selected');
    end
end
end

function PlotSurface(MC)
Xvals=[]; %Field
Yvals=[]; %Angle
Zvals=[]; %DF

RaworDeriv = 'Raw';
if strcmp(RaworDeriv,'Raw')
    for i=1:length(MC)
        lim = -.0001;
        keepinds = diff(MC{i}.Field)<lim;
        FieldTemp = MC{i}.Field(keepinds);
        FreqHold = MC{i}.Frequency(keepinds);
        FreqTemp = FreqHold - mean(FreqHold(FieldTemp<2));
        
        InterpField = linspace(0,max(FieldTemp),601)';
        InterpFreq = Interp1NonUnique(FieldTemp,FreqTemp,InterpField);
        %     Angles = repmat(MC{i}.Temp,length(InterpField),1);
        Angles = MC{i}.Temp;
        
        Xvals=[Xvals;InterpField];
        Yvals=[Yvals;Angles];
        %     tempz=Axiscell{i}.diffwithang(Axiscell{i}.Field>fieldlim)-Axiscell{subtractval}.diffwithang(Axiscell{i}.Field>fieldlim);
        %     if i<=7
        %     tempz=Axiscell{i}.diffwithang(Axiscell{i}.Field>fieldlim).*(.003076/.005665);
        %     else
        %     tempz=Axiscell{i}.diffwithang(Axiscell{i}.Field>fieldlim);
        %     end
        
        Zvals=[Zvals InterpFreq];
    end
end
if strcmp(RaworDeriv,'Deriv')
    for i=1:length(MC)
        InterpField = linspace(0,max(MC{i}.FieldDeriv),601)';
        InterpDerivFreq = Interp1NonUnique(MC{i}.FieldDeriv,...
            MC{i}.FreqDeriv,InterpField);
        
        %     Xvals=[Xvals Axiscell{i}.Temp];
        Yvals=[Yvals;MC{i}.Temp];
        Zvals=[Zvals InterpDerivFreq];
    end
end
figure;
surf(InterpField,Yvals,Zvals'); shading interp; view(2);
% colormap(hsv);
colormap(flipud(hot));
Label_Plot('\mu_0H','T','Angle','deg'); zlabel('\Delta f');
end

