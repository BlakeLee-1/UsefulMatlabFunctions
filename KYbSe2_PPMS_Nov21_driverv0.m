function output = KYbSe2_PPMS_Nov21_driverv0(option, CXfitcell, HtrRvTfit, varargin)

numvarargs = length(varargin);
if numvarargs > 2
    error('TooManyInputs', ...
        'Anarequires at most 2 optional inputs');
end

optargs = {0,0};
optargs(1:numvarargs) = varargin;
[recalcPrimary, recalcSecondary] = optargs{:};


%% Set Geometric factor
w = .98e-3;
h = 45e-6;
l = 1.03e-3;

AoL = (w*h)./l;


%% Specify the number of distinct cooldowns (CD)
fr = {'Tswp_Hllab','Tswp_Hllc'};
%fr = {'Tswp_90deg_outofplane'};
numCD = length(fr);
shape = {'s','d','o'};

%% Cell name for primary calculations and Fields measured (to index over)
cn = 'PPMS_May2019';
Fields = {'0T','18T'};
FieldNum = [0, 18, 0, 18];
leg = {'0 T', '18 T (H||ab)', '18 T (H||c)'};

if recalcPrimary
    i = 0;
    for m = 1:numCD
        %% Fileroot for loading data
        Fileroot1 = ['C:\Users\LeeLabLaptop\Documents\ML_Nov21_Transport\KYbSe2\' fr{m} '\'];
        %Fileroot1 = '/Users/iCAPOS/Desktop/Lab Work/';
        
        %Fileroot2 = '/Users/iCAPOS/Google Drive/Ba4Ir3O10/Tests_O2/';
       % Fileroot1 = [Fileroot1 '/'];
        
        cd(Fileroot1)
        D = dir;
        D = D(4:end-1);
        
        [~,ind] = sort([D.datenum]);
        D = D(ind);
        
        for l = 1:length(dir)-4
            Data{l} = strrep(D(l).name, '.mat', '') ;
        end
        
        FieldIndexing{m} = cell(1,length(Fields));
        
        
        %% Calculations
        for l = 1:length(Data)
            i = i+1;
            filename = [Fileroot1 Data{l}];
            
            %% Sort the index i of file by Field (Using string in filename)
            for j = 1:length(Fields)
                test = strfind(filename,[Fields{j} '_']);
                if not(isempty(test))
                    Bcell(i) = FieldNum(j);
                    criterion = true;
                    FieldIndexing{m}{j} = [ FieldIndexing{m}{j} i];
                end
                
            end
            
            %% Load the file data into a cell
            celld = load(filename);
            %display(filename)
            
            %% Fill datacell with all the raw data & primary calculated values (T1, T2, ... misc Temps)
            datacell(i).SampPerStep = celld.SampPerStep;
            try 
                datacell(i).HtrCurrent = celld.HtrCurrent;
            catch
               datacell(i).HtrCurrent = nan; display(filename) 
            end
                
                datacell(i).Time = celld.Time;
            %datacell(i).BathRes = celld.BathRes;
            %datacell(i).BathTemp = mean(celld.BathTemp);
%            datacell(i).ProbeTemp = mean(celld.ProbeTemp);
            datacell(i).filename = filename;
            
           
                datacell(i).CX2Res = celld.ColdRes;
                datacell(i).CX1Res = celld.HotRes;
%                datacell(i).CX3Res = celld.TransverseRes;
           
            B = Bcell(i);
            
            datacell(i).CX2Temp = celld.ColdTemp;
            datacell(i).CX1Temp = celld.HotTemp;
            
            %h = round(B);
            datacell(i).Bval = B;
            %fitInd = find(CXfitcell{5} == B);
            %display(h); display(fitInd)
            
            %datacell(i).CX2Temp = CXfitcell{2}(datacell(i).CX2Res);
            %datacell(i).CX1Temp = CXfitcell{1}(datacell(i).CX1Res);
            %datacell(i).CX2Temp = CXfitcell{2}{fitInd}(datacell(i).CX2Res.*fac_H);
            %datacell(i).CX1Temp = CXfitcell{1}{fitInd}(datacell(i).CX1Res);
           % datacell(i).CX1_Exc = celld.InputA(5);
           % datacell(i).CX2_Exc = celld.InputB(5);
        end
        clear('Data');
    end
    FieldIndexing{1}{1} = setdiff(FieldIndexing{1}{1}, FieldIndexing{1}{2});
    save(['C:\Users\LeeLabLaptop\Documents\ML_Nov21_Transport\KYbSe2\' 'Primary_datacell.mat'],'datacell');
    save(['C:\Users\LeeLabLaptop\Documents\ML_Nov21_Transport\KYbSe2\' 'Primary_FieldIndexing.mat'],'FieldIndexing');
    
else % (Load the primary calculations)
    d1 = load(['C:\Users\LeeLabLaptop\Documents\ML_Nov21_Transport\KYbSe2\' 'Primary_datacell.mat']);
    datacell = d1.datacell;
    
    d2 = load(['C:\Users\LeeLabLaptop\Documents\ML_Nov21_Transport\KYbSe2\' 'Primary_FieldIndexing.mat']);
    FieldIndexing = d2.FieldIndexing;
    
    end




%% Secondary Calcualtions
secondaryVar = {'K', 'T_CXavg', 'R_HTR', 'MaxDelT'};

if recalcSecondary
    L = length(datacell);
    
    %% Initialize Data Arrays
    for q = 1:length(secondaryVar)
        eval([secondaryVar{q} '=zeros(L,1);']);
    end 
    
    for i = 1:L
        
        samp = datacell(i).SampPerStep;
        %fac_C = mean(datacell(i).CX1Res(1:samp)./datacell(i).CX2Res(1:samp));
        fac_C = 1;
        
        datacell(i).CX2Temp = CXfitcell{2}(datacell(i).CX2Res*fac_C);
        datacell(i).CX1Temp = CXfitcell{1}(datacell(i).CX1Res);
%        datacell(i).CX3Temp = CXfitcell{3}(datacell(i).CX3Res);
        
        % Temperatures
        T1 = datacell(i).CX1Temp;
        T2 = datacell(i).CX2Temp;
        %T3 = datacell(i).CX3Temp;
        CX_R1(i) = mean(datacell(i).CX1Res);
        CX_R2(i) = mean(datacell(i).CX2Res);
        %CX_Tbath = mean(datacell(i).CX2Res);
        
        delT = T1 - T2;
        %T_bath(i) = datacell(i).BathTemp;
%        T_probe(i) = datacell(i).ProbeTemp;
       
        fac = 1;
        % Heater Power
        
        
        % Expt. Time
        t = datacell(i).Time;
        
        % Choosing indices such that points in the vicinity of HTR current changes are exlcuded
        
        
        exclude = round(.5*samp);
        %indices = [1:2*samp];
        %indices = [1:samp (2*samp+1):3*samp];%length(delT)];
        indices = [1:length(delT)];
        
        T1 = datacell(i).CX1Temp(indices);
        T2 = datacell(i).CX2Temp(indices);
        
        T_CXavg(i) = mean([T1; T2]);
        
        R_HTR(i) = HtrRvTfit(T_CXavg(i));
        Q = fac.*R_HTR(i).*power(datacell(i).HtrCurrent, 2);
        
         MaxI(i) = max(datacell(i).HtrCurrent);
        
        a={};
        for ii=1:exclude
            a{ii} = indices(ii:samp:length(indices));
        end
        b =[];
        for ii=1:exclude
            b = [b a{ii}];
        end
        %indices = [exclude:samp (samp+exclude):(2*samp)];% 
        indices = setdiff(indices, b);
        
        % Extracting the drifting nominal Delta T baseline
        indBase = [exclude:samp 4*samp+exclude:5*samp];
        
        %display(datacell(i).filename);
        %try            
        base = polyfit(t(indBase),delT(indBase)',1);
        %display(base);
       % display(datacell(i).filename);
        % Fitting Q vs Delta T to extract Kappa, MaxDelT and Rsq
        delT = delT - polyval(base, t);
        delT = delT(indices);
        %[c,s] = polyfit(Q(indices),delT,1);
        %K(i) = (1/c(1))/AoL;
        [c,s] = polyfit(Q(indices),delT,1);
        K(i) = (1/c(1))/AoL;
        
        b = Q(indices); yCalc = polyval(c,delT);
        Rsq(i) = 1 - sum((b - yCalc).^2)/sum((b - mean(b)).^2);
        
        MaxDelT(i) = max(delT);
        %end
    end
    
    for q = 1:length(secondaryVar)
        eval(['datacell2.' secondaryVar{q} '=' secondaryVar{q} ';']);
    end
    save(['C:\Users\LeeLabLaptop\Documents\ML_Nov21_Transport\KYbSe2\' 'Secondary_datacell.mat'],'datacell2');
    
else % (Load the secondary calculations)
    
    d3 = load(['C:\Users\LeeLabLaptop\Documents\ML_Nov21_Transport\KYbSe2\' 'Secondary_datacell.mat']);
    datacell2 = d3.datacell2;
    for q = 1:length(secondaryVar)
        eval([secondaryVar{q} '=' 'datacell2.' secondaryVar{q} ';']);
    end
    
end

%% Plotting options, variations on Thermal Coductance I^2R/DeltaT versus the Temperature
colorz = varycolor(length(FieldNum));



switch option
    
    case 'Plot_KvT_all'
        figure;
        
        
        
        for m = 1:numCD
            for k = 1:length(Fields)
                display([m,k])
 
                X = T_CXavg(FieldIndexing{m}{k});
                Y = K(FieldIndexing{m}{k});
                %loglog(X,Y,'sq','Color',colorz(k,:), 'Marker', shape{m}); hold on
                plot(X,Y,'sq', 'Marker', shape{m}); hold on
            end
            grid on;
            
        end
        %xlim([0,150]) ;
        legend(leg);
        xlabel('T[K]'); ylabel('\kappa [W/m K]')
        
        
    case 'Plot KvT clean'
        figure;
        
        m = 1; k = 1;
        X = T_CXavg(FieldIndexing{m}{k});
        Y = K(FieldIndexing{m}{k});
        excind = [1,77];
        ind = setdiff(1:length(X),excind);     
        plot(X(ind),Y(ind),'sq', 'Marker', shape{m}); hold on
        
        
        m = 1; k = 2;
        X = T_CXavg(FieldIndexing{m}{k});
        Y = K(FieldIndexing{m}{k}); ind2 = 196:203; ind3 = 249:258;
        excind = [67,192:195, ind2, 244:248];
        ind = setdiff(1:length(X),excind);
        plot([X(ind)' X(ind2)'],[Y(ind)' Y(ind2)'/1.0263],'sq', 'Marker', shape{m}); hold on
       
        
        m = 2; k = 2;
        X1 = T_CXavg(FieldIndexing{m}{k});
        Y1 = K(FieldIndexing{m}{k});
        excind = [62 206:228];
        ind = setdiff(1:length(X1),excind);     
        plot([X1(ind)' X(ind3)'], [Y1(ind)' Y(ind3)'],'sq', 'Marker', shape{m}); hold on
        

        %xlim([0,150]) ;
        legend(leg);
        xlabel('T[K]'); ylabel('\kappa [W/m K]')
        
    case 'Plot_KvT'
        figure;
        
        ExcInd{1}{1} = [];%[1:112 461:528];
        ExcInd{1}{2} = [];
        ExcInd{1}{3} = [];
        ExcInd{2}{1} = [];
        ExcInd{2}{2} = [];
        ExcInd{2}{3} = [];
        ExcInd{3}{1} = [];
        ExcInd{3}{2} = [];
        ExcInd{3}{3} = [];
        
        
%        ExcInd{1}{2} = [269:315, 385:428, 432:445, 61, 87, 93, 99, 53];

        
        
        for m = 1:numCD
            for k = 1:length(Fields)
                %display([m,k])
                ind1 = FieldIndexing{m}{k};
                X = T_CXavg(ind1);
                Y = K(ind1)/10;
                
                ind2 = setdiff(1:length(X), ExcInd{m}{k});
                X = X(ind2);
                Y = Y(ind2);
                
                %loglog(X,Y,'sq','Color',colorz(k,:), 'Marker', shape{m}); hold on
                plot(X,Y,'sq', 'Marker', shape{m}); hold on
            end
            grid on;
            
        end
        xlim([0,150]) ;
        legend(leg);
        xlabel('T[K]'); ylabel('\kappa [W/m K]')
            
    case '18TMR'
        
%        a =       121.6;
%        c =       4.451;
%        d =      0.3939;
%        f =     0.05391;
       a =       129.8;
       c =       4.469;
       d =      0.3924;
       f =     0.05044;

       F = @(x) a.*cos(f.*x+c).*exp(-d.*x);
        
        k = length(Fields);
        tol = .001;
        figure;
        for m = 1:numCD
            ind = FieldIndexing{m}{k};
            for i = ind
                 R1B = mean(datacell(i).CX1Res)  ;  R2B = mean(datacell(i).CX2Res); 
                 R1ZFprev = 0    ;    R2ZFprev = 0;
                 R1ZFcurr = R1B  ;  R2ZFcurr = R2B;
                % display(i)
                while (abs(R1ZFprev-R1ZFcurr)>tol) && (abs(R2ZFprev-R2ZFcurr)>tol)
                    %\
                    temp = R1ZFcurr;
                    f = F(CXfitcell{1}(temp))./100;
                    R1ZFcurr = R1B/(1+f);
                    R1ZFprev = temp;
                    %
                    temp = R2ZFcurr;
                    f = F(CXfitcell{2}(temp))./100;
                    R2ZFcurr = R2B/(1+f);
                    R2ZFprev = temp;
                end
                T1(i) = CXfitcell{1}(R1ZFcurr);
                T2(i) = CXfitcell{2}(R2ZFcurr); 
            end
            X = (T1(ind)+T2(ind))./2;
            Y = K(ind);
            loglog(X,Y,'sq'); hold on
        end    
        
    case 'RCXvT'
        %figure; %loglog(nan,nan);
        %hold on
        Tinterp = [5:1:90];
        X = [];
            Y1 = [];
            Y2 = [];
        for k = 1:length(Fields)
            
            for m = 1:numCD
                %X = [X T(FieldIndexing{m}{k})];
                X  = [X T_CXavg(FieldIndexing{m}{k})'];
                Y1 = [Y1 CX_R1(FieldIndexing{m}{k})];
                Y2 = [Y2 CX_R2(FieldIndexing{m}{k})];
                %Y = [Y MaxI(FieldIndexing{m}{k})];
            end
            %plot(X,Y,'sq','Color',colorz(k,:))
            %plot(X,Y,'sq','Color','r')
            grid on;
            
            
        end
        output = [X; Y1; Y2];
        xlim([0,110])
        %sylim([0,50])
        legend(Fields)
        xlabel('T[K]'); ylabel('\kappa W/mK') 
        
    otherwise
        display('invalid option')
        
end

end