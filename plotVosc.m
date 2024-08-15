function plotVosc()
    temp = [4,5:1:10, 12:2:20, 22:4:50];
    hold off;
    for i=1:length(temp)
        i;
        M(i) = LoadData(temp(i));
        Plot(M(i),temp(i)) ;
        hold on;
    end

end




   function [M] = LoadData(temp) %Function to create a structure M
    grid on;
    fileroot = '/home/blake/Documents/MATLAB/SCM2/Fswp/'; %Parent file location for our data
    cd(fileroot); %We need to get there
    filedata = dir(['*_',int2str(temp),'K*.mat']); %define structure filedata as all files containing string
     %for all of these elements
       %filedata.name
       M = load(filedata.name);  % Load all data first before adding fields so they are all the same length
    
     %straight up addressing all of those named elements from last loop and giving them these attributes
       M.name = filedata.name;
       M.Fileroot = fileroot;
       M.DateStr = filedata.date;
       M.Datenum = filedata.datenum;
   
   end

function [M] = Plot(M, temp)
    M.Vmod = [];
    Freqs = unique(M.Frequency);
    tol = 0.35;
    num = 170;
    M.ModVOsc_Times_Frequency = []


 


    for i = 1:length(Freqs)
        selection = M.VOsc_Times_Frequency((M.PointsPerFrequency*i-M.PointsPerFrequency*tol:M.PointsPerFrequency*i-M.PointsPerFrequency*tol+num));
        wholeselec = selection;
        for i = length(selection):-1:1
            if abs(selection(i) - mean(wholeselec)) > 0.1*std(wholeselec)
                selection(i) = [];
            end
        end
        M.Vmod = [M.Vmod, mean(selection)];
    end
    %for i=1:length(M.Vmod)
        for i = length(Freqs):-1:1
            if Freqs(i) == 90 || Freqs(i) == 95
                Freqs(i) = [];
                M.Vmod(i) = [];
            end
        end
        if temp == 30
            M.Vmod(27) = [];
            Freqs(27) = [];
        end
        M.VmodUn = M.Vmod;
        
        for j = 1:0
        for i = length(M.Vmod)-1:-1:2
            M.Vmod(i) = (M.Vmod(i+1)+M.Vmod(i)+M.Vmod(i-1))/3;
        end

        end 
        temp
        M.Vmod./mean(M.Vmod)
    plot(Freqs, M.Vmod./mean(M.Vmod), 'Color', ColormapInterpolate_InBounds(temp, [0,70]), 'LineWidth', 3, 'DisplayName', [int2str(temp),'K']); %ColormapInterpolate_InBounds(temp, [0,70]))
    %plot(Freqs, M.VmodUn./mean(M.VmodUn), 'Color', ColormapInterpolate_InBounds(temp, [0,70]), 'LineWidth', 3, 'DisplayName', [int2str(temp),'K']); %ColormapInterpolate_InBounds(temp, [0,70]))

end


% Freq = [175, 162, 158, 153, 150, 145, ]
%Temper = [5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 26, 30, 34, 38, 42, 46, 50]