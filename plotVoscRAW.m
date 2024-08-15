function plotVoscRAW()
    temp = [4,5:1:10, 12:2:20, 22:4:50];
    hold off;
    for i=1:length(temp)
        i;
        M(i) = LoadData(temp(i));
        Plot(M(i),temp(i)) ;
    end
 hold off;

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
    
    for i = 1:length(Freqs)
        M.Vmod = [M.Vmod, mean(M.VOsc_Times_Frequency((M.PointsPerFrequency*i-200:M.PointsPerFrequency*i)))];
    end
    hold on;
    %for i=1:length(M.Vmod)
        if length(Freqs) == 29 && Freqs(27) == 90
            Freqs;
            M.Vmod(27) = [];
            Freqs(27) = [];
        end
        for i = length(M.Frequency):-1:1
            if M.Frequency(i) == 90 || M.Frequency(i) == 95
                M.Frequency(i) = [];
                M.VOsc_Times_Frequency(i) = [];
            end
        end
    
    unique(M.Frequency)           
    plot(M.Frequency, M.VOsc_Times_Frequency, 'Color', ColormapInterpolate_InBounds(temp, [0,70]), 'LineWidth', 3, 'DisplayName', [int2str(temp),'K']); %ColormapInterpolate_InBounds(temp, [0,70]))
end