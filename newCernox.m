function [f,M] = newCernox()
    
    M0 = LoadData('0');
    M18 = LoadData('18');
    %M12 = LoadData('12');
    kappa (M18, 'r', '18T')
    M18 = LoadData('18T Perpendicular');
    %kappa(M12,'g')
    %ProcessData(M);
    kappa(M0,'b', '0T')
    %[f,~] = resistanceCalc(M);
    
end


function [M] = LoadData(temp) %Function to create a structure M
    fileroot = '/home/blake/Documents/MATLAB/Tswp/'; %Parent file location for our data
    cd(fileroot); %We need to get there
    filedata = dir(['*',temp,'T*.mat']); %define structure filedata as all files containing string
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
                
function [M] = ProcessData(M) %Creates a structure M taking in Tolerance and a different M structure
    for i = 1:length(M)
        firstStepRes(i) = mean(M(i).BathRes(1:M(i).SampPerStep));
        firstStepTemp(i) = mean(M(i).BathTemp(1:M(i).SampPerStep));
        firstStepResHot(i) = mean(M(i).HotRes(1:M(i).SampPerStep));
        firstStepResCold(i) = mean(M(i).ColdRes(1:M(i).SampPerStep));
    end
    
    %hold on;plot(firstStepTemp,firstStepResHot,'or', 'MarkerFaceColor', 'r');
    %plot(firstStepTemp,firstStepResCold,'ob', 'MarkerFaceColor', 'b');
    Hotstuff = readmatrix("/home/blake/Documents/MATLAB/Cali/CX_M_CrI3_Hot.csv");
    Coldstuff = readmatrix("/home/blake/Documents/MATLAB/Cali/CX_N_CrI3_Trans.csv");

    %plot(Hotstuff(:,2), Hotstuff(:,1), 'og', 'MarkerFaceColor', 'g');
    %plot(Coldstuff(:,2), Coldstuff(:,1), 'ok', 'MarkerFaceColor', 'k');

    %HotFit = RtoT_cal_InputRT(firstStepTemp,firstStepResHot,'fitRvT',7,'Yes','Polynomial');
    %ColdFit = RtoT_cal_InputRT(firstStepTemp,firstStepResCold,'fitRvT',7,'Yes','Polynomial');;
   
end

function [M] = kappa(M, color, name)
    kappas = [];
    BathTemp = [];
    for i = 1:length(M)
        Tolerance=0.4;
        
        a = logical([zeros(1,floor(M(i).SampPerStep*(1-Tolerance))), ones(1,floor(M(i).SampPerStep*(Tolerance)))]);
        b = logical([repmat(a,1,5)]); %Cuts all the non-flat bits


        M(i).cutHotTemp = M(i).HotTemp(b);    %cut M temperatures
        M(i).cutColdTemp = M(i).ColdTemp(b);  %Cut N temperatures
        M(i).cutCurrent = M(i).HtrCurrent(b); %Cut Htr Current
        M(i).cutBathTemp = M(i).BathTemp(b);

        M(i).delT = M(i).cutHotTemp - M(i).cutColdTemp;  %Convert to this run's Temp difference
        M(i).powers = (M(i).cutCurrent.*M(i).cutCurrent).*M(i).HTRres;    %Square this run's currents
        [M(i).fitdata, M(i).error] = polyfit(M(i).powers,M(i).delT, 1);
        M(i).rsq = 1 - (M(i).error.normr/norm(M(i).delT - mean(M(i).delT)))^2;
        M(i).kappa = M(i).fitdata(1);
        M(i).intercept = M(i).fitdata(2);
        kappas = [kappas, M(i).kappa];
        BathTemp = [BathTemp, (mean(M(i).cutHotTemp)+mean(M(i).cutColdTemp))/2];



    end
    
    area = 5.25e-8;
    Length = 1.875e-4;
    alpha = Length/area;
    deleted = [];
    kappapure = kappas;
   for i = 1:length(kappas)
       if (kappas(i)/alpha)^(-1) < -0.1
           deleted = [deleted, i];
       end
   end
   

   for i = length(deleted):-1:1
      kappas(deleted(i)) = [];
      BathTemp(deleted(i)) = [];
   end
    hold on;
    plot(BathTemp, (kappas/alpha).^(-1), ['o',color], 'MarkerFaceColor', color, 'DisplayName', name)
    xlabel('Sample Temperature (K)')
    ylabel('\kappa (W/Km)')
end

function [f,M] = resistanceCalc(M)

    heat = [];
    res = [];
    heat2 = [];

    for i = 1:length(M)
        heat = [heat,mean(M(i).BathTemp)];
        heat2 = [heat2,(mean(M(i).HotTemp)+mean(M(i).ColdTemp))/2];
        res = [res, M(i).HTRres];
        
    end

   deleted = [];
   cutoff = 1500;
   for i = 1:(length(res)-3)
       if (res(i)-res(i+1)) > cutoff 
            deleted = [deleted, i];    
       elseif (res(i)-res(i+2)) > cutoff
            deleted = [deleted, i];
       elseif (res(i)-res(i+3)) > cutoff
            deleted = [deleted, i];
       end
   end
   

   for i = length(deleted):-1:1
      res(deleted(i)) = [];
      heat(deleted(i)) = [];
      heat2(deleted(i)) = [];

   end

   logres = log(res);
   logheat = log(heat);
% 
    f = fit(heat2',res','power2');
    hold off;
    %plot(heat,res,'or')
    %hold on;
    %plot(heat2,res,'or')

    plot(f, heat2', res')
    %plot(heat,f(heat)'-res, 'ok')
end