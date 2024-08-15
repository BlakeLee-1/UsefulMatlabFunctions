function [M] = BswpAnalysis()
    
    M5 = LoadData(5);
    M10 = LoadData(10);
    M15 = LoadData(15);
    M20 = LoadData(20);
    process(M5, 5);
    process(M10, 10);
    process(M15, 15);
    process(M20, 20);
    

    

end





function [M] = LoadData(temp) %Function to create a structure MS


    fileroot = '/home/blake/Documents/MATLAB/Bswp/'; %Parent file location for our data
    cd(fileroot); %We need to get there
    if temp == 5
        filedata = dir('MgCrO110_ML_CD1_5K__2022_5_18_13_51_48.mat');
    elseif temp == 10
        filedata = dir('MgCrO110_ML_CD1_10p0002K__2022_5_18_14_55_13.mat');
    elseif temp ==15
        filedata = dir('MgCrO110_ML_CD1_15p0012K__2022_5_18_16_21_52.mat');
    elseif temp ==20
        filedata = dir('MgCrO110_ML_CD1_20p0006K__2022_5_18_18_22_18.mat');
    end
    for i = 1:length(filedata)  %for all of these elements
       M = load(filedata(i).name);  % Load all data first before adding fields so they are all the same length
    end
    for i = 1:length(filedata) %straight up addressing all of those named elements from last loop and giving them these attributes
       M.name = filedata(i).name;
       M.Fileroot = fileroot;
       M.DateStr = filedata(i).date;
       M.Datenum = filedata(i).datenum;
    end
end

function [M,delT] = process(M, temp)
        
    HTRres = max(M.HTRres);
    HTRres
    current = max(M.HtrCurrent);
    powers = HTRres*current^2;
    kappas = [];
    BathTemp = [];
        Tolerance=0.5;
        n = floor(length(M.HotTemp)/M.SamplesPerStep);
        a = logical([zeros(1,floor(M.SamplesPerStep*(1-Tolerance))), ones(1,floor(M.SamplesPerStep*(Tolerance)))]);
        b = logical([repmat(a,1,n)]); 


        CutFields = {'Bfield','HtrCurrent','Time','VTITemp','ProbeTemp','HotRes', ...
            'ColdRes','BathRes','HotTemp','ColdTemp','BathTemp'};
        for q = 1:length(CutFields)
            M.(['cut',CutFields{q}]) = M.(CutFields{q})(b);
        end
        
        M.delT = M.cutHotTemp - M.cutColdTemp;
        delT = [];
        kappaB = [];
        for i = 1:n
            kappaB = [kappaB, mean(M.cutBfield(((i-1)*M.SamplesPerStep*Tolerance + 1):i*M.SamplesPerStep*Tolerance))];
            delT = [delT, mean(M.delT(((i-1)*M.SamplesPerStep*Tolerance + 1):i*M.SamplesPerStep*Tolerance))];
        end

        xvals = [];
        for i = 1:n
            xvals = [xvals, i*M.SamplesPerStep*Tolerance];
        end
        
        deldelT = [];

        if mod(length(delT),2) == 1
            for i = 1:2:(length(delT)-1)
                deldelT = [deldelT, delT(i+1)-delT(i)];
            end
        elseif mod(length(delT),2) == 0
            for i = 1:2:length(delT)
                deldelT = [deldelT, delT(i+1)-delT(i)];
            end  
        end

        area = 5.25e-8;
        Length = 1.875e-4;
        alpha = Length/area;
        kappas = (deldelT/(powers*alpha)).^(-1);
        kappaBhalf = [];

        for i = 1:2:(length(kappaB))-1
            kappaBhalf = [kappaBhalf, (kappaB(i)+kappaB(i+1))/2];
        end

        if temp == 5
            hold on;
            plot(kappaBhalf, kappas, 'or', 'MarkerFaceColor', 'r', 'DisplayName', '5')
        elseif temp == 10
            hold on;
            plot(kappaBhalf, kappas, 'og', 'MarkerFaceColor', 'g', 'DisplayName', '10')
        elseif temp ==15
            hold on;
            plot(kappaBhalf, kappas, 'ob', 'MarkerFaceColor', 'b', 'DisplayName', '15')
        elseif temp ==20
            hold on;
            plot(kappaBhalf, kappas, 'ok', 'MarkerFaceColor', 'k', 'DisplayName', '20')
        end
        
end

