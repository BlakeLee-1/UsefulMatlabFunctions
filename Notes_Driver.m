%% Fit code copypaste

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


figure('name','kappatry2')
plot(MS(1).HeaterPowerStepsCut,MS(1).DeltaTStepsCut,'or'); 
    hold on; plot(MS(500).HeaterPowerStepsCut,MS(500).DeltaTStepsCut,'og');
    hold on; plot(MS(1000).HeaterPowerStepsCut,MS(1000).DeltaTStepsCut,'ob');
    Label_Plot('Power','W','Delta T','K');  
    

%% Checking MS(850) data delta T vs time, delta T vs current

figure('name','MS850 time check');
   
            subplot(1,2,1); 
            plot(MS(850).Time,MS(850).DeltaT,'or');
            Label_Plot('Time','s','Delta T','K');
            
            subplot(1,2,2); 
            plot(MS(850).HeaterPower,MS(850).DeltaT,'ob');
            Label_Plot('Power','W','Delta T','K');
            
            drawnow
            
%% Checking MS(850) data 
figure('name','MS850');
   
            subplot(1,3,1); 
            plot(MS(850).HeaterPower,MS(850).DeltaT,'or');
            Label_Plot('Power','W','Delta T','K');
            
            subplot(1,3,2); 
            plot(MS(850).HeaterPowerStepsOnly,MS(850).DeltaTStepsOnly,'ob');
            Label_Plot('Power','W','Delta T','K');
            
            subplot(1,3,3); 
            plot(MS(850).HeaterPowerStepsCut,MS(850).DeltaTStepsCut,'og');
            Label_Plot('Power','W','Delta T','K');
 
            drawnow
            
 
%%  Comparing different MS(i) data
figure('name','Sorted Comparing Currents');
   
            plot(MS(1).HeaterPowerStepsCut,MS(1).DeltaTStepsCut,'or'); hold on;
            plot(MS(150).HeaterPowerStepsCut,MS(150).DeltaTStepsCut,'om'); hold on;
            plot(MS(300).HeaterPowerStepsCut,MS(300).DeltaTStepsCut,'oy'); hold on;
            plot(MS(450).HeaterPowerStepsCut,MS(450).DeltaTStepsCut,'og'); hold on;
            plot(MS(600).HeaterPowerStepsCut,MS(600).DeltaTStepsCut,'oc'); hold on;
            plot(MS(850).HeaterPowerStepsCut,MS(850).DeltaTStepsCut,'ob'); hold on;
            Label_Plot('Power','W','Delta T','K');
 
            drawnow
            
%% Getting sample avg temp and sorting it
for i = length(MS)
    MS(i).TotalAvgerageTemp = mean(MS(i).SampleAverageTemp);     % Create new field for single Total Average Temp per .mat data file
end
TotalAverageTemp = vertcat(MS.TotalAvgerageTemp);                % Create Total Average Temp vector from structure field
[Temp,Index] = sort(TotalAverageTemp);                           % Sort Total Average Temp vector and index values 
MS = MS(Index)

%% Use this to save figures, Saahhh
SavePngFig(figfr,'KVO_300K_Ex100_Sens100_400SPS')

%% Recalibrating Cernox's
% Load the SD cernox calibrations. 
% BathCell = {SD1, SD2}; 

% make new calibrations using x = bath temp, y = hot/cold cernox temps
BathCell = load('C:\Users\LeeLabLaptop\Documents\Minhyea Lee Research\Maglab_March2019\CXfits_March2019\SDCXcell.mat');
BathFit = BathCell.CXcell{1};

deg = 7;
CX_fr = 'C:\Users\LeeLabLaptop\Documents\Sarah\CXData\MN_3_2022\'; 
cd(CX_fr);
TLim = 70; 
Dathot = load('CX_M_CrI3_Hot.csv');  % Hot Cx Fit
x = Dathot(:,1); y = Dathot(:,2);
KeepInds = y>=TLim; 
x = x(KeepInds); y = y(KeepInds);
HotFit = RtoT_cal_inputRT(x,y,'getfit',deg,'Yes','Polynomial',[]);
Datcold = load('CX_N_CrI3_Trans.csv');   % Cold Cx Fit
x = Datcold(:,1); y = Datcold(:,2);
KeepInds = y>=TLim; 
x = x(KeepInds); y = y(KeepInds);
ColdFit = RtoT_cal_inputRT(x,y,'plotRvT',deg,'Yes','Polynomial',[]);