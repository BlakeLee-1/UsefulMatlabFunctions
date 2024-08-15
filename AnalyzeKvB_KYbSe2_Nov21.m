function output = AnalyzeKvB_KYbSe2_Nov21(option,CXcell,HtrFit,set)

%% Load Data

if set
    folders = 'C:\Users\LeeLabLaptop\Documents\ML_Nov21_Transport\KYbSe2\Bswp\';
    
    fl =   	  {...
        'KYbSe2_1p769K__2021_11_21_11_52_39.mat',...
        'KYbSe2_3p9971K__2021_11_21_14_31_1'...
        'KYbSe2_4p9938K__2021_11_21_21_31_2.mat',...
        'KYbSe2_5p9953K__2021_11_21_15_17_8.mat',...
        'KYbSe2_8p0039K__2021_11_21_16_2_4.mat',...
        'KYbSe2_10p0019K__2021_11_21_16_41_5.mat',...
        'KYbSe2_11p9973K__2021_11_21_17_17_12.mat',...
        'KYbSe2_15p0939K__2021_11_21_18_15_50.mat',...
        'KYbSe2_19p9989K__2021_11_21_18_56_4.mat',...
        'KYbSe2_29p9968K__2021_11_21_19_42_32.mat',...
        'KYbSe2_49p9765K__2021_11_21_20_40_2.mat',...
        };
    
    leg = {'1.6 K (H || ab)', '4 K', '5 K', '6 K', '8 K', '10 K', '12 K', '15 K', '20 K', '30 K', '50 K'};
    
else
    folders = 'C:\Users\LeeLabLaptop\Documents\ML_Nov21_Transport\KYbSe2\Bswp\';
    
    fl =   	  {...
    'KYbSe2_Hllc_1p773K__2021_11_22_9_20_27.mat',...
    'KYbSe2_Hllc_3p9983K__2021_11_22_10_36_5.mat',...
    'KYbSe2_4p9993K__2021_11_21_22_14_41.mat',...
    ...'KYbSe2_Hllc_5p0023K__2021_11_22_12_18_.mat',...
    'KYbSe2_Hllc_5p9923K__2021_11_22_8_17_8.mat',...
    'KYbSe2_Hllc_8p004K__2021_11_22_7_40_4.mat',...
    'KYbSe2_Hllc_9p9963K__2021_11_22_7_2_10.mat',...
    'KYbSe2_Hllc_11p9996K__2021_11_22_6_23_59.mat',...
    'KYbSe2_Hllc_14p9993K__2021_11_22_5_46_50.mat',...
    'KYbSe2_Hllc_20p0011K__2021_11_22_5_6_2.mat',...
    'KYbSe2_Hllc_30p0171K__2021_11_22_4_25_9.mat',...
    'KYbSe2_Hllc_50p0624K__2021_11_22_3_37_20.mat',...
        };
    
    leg = {'1.6 K (H || ab)', '4 K', '5 K', '6 K', '8 K', '10 K', '12 K', '15 K', '20 K', '30 K', '50 K'};
    

    
end


fieldX = [0,6,9,12,15,18];

for i = 1:length(fl)
    datacell{i} = load([folders fl{i}]);%fl1;
end

L = length(fl);

factor = 1.0;


w = .98e-3;
h = 45e-6;
l = 1.03e-3;

AoL = (w*h)./l;


ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.TolFun = 1e-09;
opts.TolX = 1e-09;

colz = varycolor(L);

%% Process Raw data
for i = 1:L
    N = datacell{i}.SamplesPerStep;
    Current = datacell{i}.HtrCurrent;
    Hot = CXcell{1}(datacell{i}.HotRes);
    Cold = CXcell{2}(datacell{i}.ColdRes);
%    Transverse = CXcell{3}(datacell{i}.TransverseRes);
    RawDelTxx = Hot - Cold;
  %  RawDelTxy = Transverse - Cold;
%    HtrResRaw = datacell{i}.HTRres;
    Bval = datacell{i}.Bfield;
    ProbeT = datacell{i}.ProbeTemp;
    %    BathT = datacell{i}.BathTemp;
    %display(mean(ProbeT))
    
    M = length(Bval);
    M = M - mod(M,2*N);
    
    %     temp = isnan(HtrResRaw);
    %     HtrResRaw(temp)=0;
    
    Ex = round(.3*N);
    
    %% Index according to Heater-On/Off and Get Heater-Off values of B
    for j = 1:M/(2*N)
        ind1{j} = (j-1)*2*N+1+Ex:j*2*N-N ;
        ind2{j} = (j-1)*2*N+N+1+Ex:j*2*N ;
        ind3{j} = (j-1)*2*N+1+1:j*2*N    ;
        cell.B_off{i}(j) = mean(Bval(ind1{j}));
    end
    
    %% Find ZF resistance values of each cernox (for a simple MR correction)
    tmp = abs(cell.B_off{i});
    ZF_ind = find(tmp == min(tmp));
    cell.R1_0{i} = mean(datacell{i}.HotRes(ind1{ZF_ind}));
    cell.R2_0{i} = mean(datacell{i}.ColdRes(ind1{ZF_ind}));
%    cell.R3_0{i} = mean(datacell{i}.TransverseRes(ind1{ZF_ind}));
    
    %% Compute Temperature Gradients
    for j = 1:M/(2*N)
        
        I1 = ind1{j};
        I2 = ind2{j};
        I3 = [I1 I2];
        
        %         I1 = [(ind1{j}-2*N) ind1{j} (ind1{j}+2*N)];
        %         I2 = [(ind2{j}-2*N) ind2{j} (ind2{j}+2*N)];
        %         I3 = [(ind3{j}-2*N) ind3{j} (ind3{j}+2*N)];
        %
        %
        %          mfDTxx = medfilt1(RawDelTxx);
        %          mfDTxy = medfilt1(RawDelTxy);
        %
        %          xx_on  = @(b) polyval(polyfit(Bval(I2), mfDTxx(I2)', 2),b);
        %          xx_off = @(b) polyval(polyfit(Bval(I1), mfDTxx(I1)', 1),b);
        %
        %          xy_on  = @(b) polyval(polyfit(Bval(I2), mfDTxy(I2)', 2),b);
        %          xy_off = @(b) polyval(polyfit(Bval(I1), mfDTxy(I1)', 1),b);
        
        h  = mean(Bval(I3));
        cell.B{i}(j) = h;
        %
        %          cell.DelTxx{i}(j) = (xx_on(h) - xx_off(h));
        %          cell.DelTxy{i}(j) = (xy_on(h) - xy_off(h));
        
        cell.R1_off{i}(j) = mean(datacell{i}.HotRes(I1));
        cell.R2_off{i}(j) = mean(datacell{i}.ColdRes(I1));
        %cell.R3_off{i}(j) = mean(datacell{i}.TransverseRes(I1));
        
        MRc1 = cell.R1_0{i}./cell.R1_off{i}(j);
        MRc2 = cell.R2_0{i}./cell.R2_off{i}(j);
%        MRc3 = cell.R3_0{i}./cell.R3_off{i}(j);
        
        cell.R1_on{i}(j) = mean(datacell{i}.HotRes(I2));
        cell.R2_on{i}(j) = mean(datacell{i}.ColdRes(I2));
        %cell.R3_on{i}(j) = mean(datacell{i}.TransverseRes(I2));
        
        cell.HotTemp{i}(j)   = mean([CXcell{1}(MRc1*cell.R1_off{i}(j)) CXcell{1}(MRc1*cell.R1_on{i}(j))]);
        cell.ColdTemp{i}(j)  = mean([CXcell{2}(MRc2*cell.R2_off{i}(j)) CXcell{2}(MRc2*cell.R2_on{i}(j))]);
       % cell.TransTemp{i}(j) = mean([CXcell{3}(MRc3*cell.R3_off{i}(j)) CXcell{3}(MRc3*cell.R3_on{i}(j))]);
        
        
        cell.DelTxx{i}(j) =  CXcell{1}(cell.R1_on{i}(j))  - CXcell{2}(cell.R2_on{i}(j)) - ...
            (CXcell{1}(cell.R1_off{i}(j)) - CXcell{2}(cell.R2_off{i}(j)) );
        
        %cell.DelTxx{i}(j) =  CXcell{1}(MRc1*cell.R1_on{i}(j))  - CXcell{2}(MRc2*cell.R2_on{i}(j)) - ...
        %                    (CXcell{1}(MRc1*cell.R1_off{i}(j)) - CXcell{2}(MRc2*cell.R2_off{i}(j)) );
        
        %cell.DelTxy{i}(j) =  CXcell{3}(MRc3*cell.R3_on{i}(j))  - CXcell{2}(MRc2*cell.R2_on{i}(j)) - ...
       %     (CXcell{3}(MRc3*cell.R3_off{i}(j)) - CXcell{2}(MRc2*cell.R2_off{i}(j)) );
        
        
        
        %         cell.T{i}(j) = mean(Hot([ind1 ind2])+Cold([ind1 ind2]))/2;
        cell.R{i}(j) = 1000;%mean(HtrResRaw(I2));
        cell.I{i}(j) = mean(Current(I2));
        %display(cell.I{i}(j))
        %
        %          = mean(Bval(ind1));
        %         cell.Tbath{i}(j) = mean(BathT([ind1 ind2]));
        %         cell.Tprobe{i}(j) = mean(ProbeT([ind1 ind2]));
        
        
    end
%     if cell.B{i}(end) > cell.B{i}(1)
%         meanTemp{i} = mean(cell.HotTemp{i}(1:4)+cell.ColdTemp{i}(1:4))./2;
%     else
%         meanTemp{i} = mean(cell.HotTemp{i}(end-3:end)+cell.ColdTemp{i}(end-3:end))./2;
%     end
    meanTemp(i) = mean([cell.HotTemp{i} cell.ColdTemp{i}]);
    display(meanTemp(i));
    
    
    for m = 1:length(fieldX)
        
        X = abs(cell.B{i}-fieldX(m));
        ind = find(X==min(X));
        k = power(cell.I{i},2).*cell.R{i}./(cell.DelTxx{i}.*AoL);
        
        KX(i,m) = k(ind);
    end
    %display(mean(cell.HotTemp{i}+cell.ColdTemp{i})./2);
    %     try
    %         TempAvg{i} = mean(cell.T{i}(:));
    %     end
end

%output = TempAvg;


%% Meta Processing



% for n = 1:L%length(PatchSweeps)
%     H{n} = []; H_off{n} = [];
%     Q{n} = []; Rfull{n} = [];
%     DelTxx{n} = []; Kxx{n} = [];
%     DelTxy{n} = []; Kxy{n} = [];
%     T_Probe{n} = []; T_Bath{n} = [];
%     R1{n} = []; R2{n} = []; R3{n} = [];
%    % figure; hold on
%     for i = RunInd{n}
%
%         Htemp = cell.B{i}(:);
%         Htemp_off = cell.B_off{i}(:);
%
%         %plot(Htemp)
%         H{n} = [H{n}(:)' Htemp'];
%         H_off{n} = [H_off{n}(:)' Htemp_off'];
%         ind = find(not(isnan(cell.R{i})));
%         R{n} = 1000;%mean(cell.R{i}(ind));
%         Rfull{n} = [Rfull{n} cell.R{i}];
%
%         Qtemp = Rfull{n}.*max(cell.I{i}(:)).^2;
%         Q{n} = [Q{n} Qtemp];
%
%         R1{n} = [R1{n} cell.R1{i}];
%         R2{n} = [R2{n} cell.R2{i}];
%         R3{n} = [R3{n} cell.R3{i}];
%
%         T_Probe{n} = [T_Probe{n} cell.Tprobe{i}];
%         T_Bath{n} = [T_Bath{n} cell.Tbath{i}];
%
%         %display(mean(cell.T{i}))
%
%         DelTxx{n} = [DelTxx{n}(:)' cell.DelTxx{i}(:)'];
%         DelTxy{n} = [DelTxy{n}(:)' cell.DelTxy{i}(:)'];
%
%         Kxxtemp = R{n}.*cell.I{i}(:).^2./((w*t/l)*cell.DelTxx{i}(:));
%
%
%         Kxx{n} = [Kxx{n}(:)' Kxxtemp'];
%
%     end
%     %ind = find(and(not(isnan(K{n})),not(K{n}==0)));
%
%     ind = find(not(isnan(Kxx{n})));
%     Kxx{n} = Kxx{n}(ind); H{n} = H{n}(ind);
%     DelTxx{n} = DelTxx{n}(ind); DelTxy{n}(ind);
%
%     [H{n} ind] = sort(H{n}); Kxx{n} = Kxx{n}(ind);
%
%     DelTxx{n} = DelTxx{n}(ind);  DelTxy{n} = DelTxy{n}(ind);
%     T_Probe{n} = T_Probe{n}(ind); T_Bath{n} = T_Bath{n}(ind);
%     Rfull{n} = Rfull{n}(ind);
%     %Q{n} = Q{n}(ind);
%
%     R1{n} = R1{n}(ind);
%     R2{n} = R2{n}(ind);
%     R3{n} = R3{n}(ind);
%
%     ZFP = find(abs(H{n})==min(abs(H{n})));
%     AvgSampT{n} = mean([CXcell{1}(R1{n}(ZFP)) CXcell{2}(R2{n}(ZFP)) CXcell{3}(R3{n}(ZFP))]);
%     %display(['Avg CX read sample temp = ' num2str(AvgSampT{n})])
%
%     HTRRes{n} = HtrFit(AvgSampT{n});
%     HTRI{n}   = max(cell.I{n});
%
%     Q{n} = power(HTRI{n},2).*HTRRes{n};
%     %display(['R = ' num2str(HTRI{n})])
% end

switch option
    
    case 'plot KX'
        figure; hold on;
        for m = 1:length(fieldX)
            plot(meanTemp, KX(:,m), 'sq-','displayname',num2str(fieldX(m)))
        end
        legend show
    
    case 'Plot KxxvB'
        
        figure; hold on;
        for n = 1:L
            K{n} = power(cell.I{n},2).*cell.R{n}./(cell.DelTxx{n}.*AoL);
            [X ind] = sort(abs(cell.B{n}));
            
            Y = smooth(K{n}(ind),1);
            plot(X, Y,'sq','color',colz(n,:))
        end
        grid on; box on; xlabel('\mu_0H [T]'); ylabel('\kappa [W/Km]')
        legend(leg,'location','eastoutside'); xlim([0,14]);
              
    case 'KxxvB out'
        
        %figure; hold on;
        for n = 1:L
            K{n} = power(cell.I{n},2).*cell.R{n}./(cell.DelTxx{n}.*AoL);
            [X ind] = sort(abs(cell.B{n}));
            
            Y = smooth(K{n}(ind),1);
            output.H{n} = X/1e4;
            output.K{n} = Y;
            output.desc{n} = leg{n};
            %plot(X, Y,'sq','color',colz(n,:))
        end
%         grid on; box on; xlabel('\mu_0H [T]'); ylabel('\kappa [W/Km]')
%         legend(leg,'location','eastoutside'); xlim([0,18]);
              
    case 'Plot KxxvB MR Cal'
        
        figure; hold on;
        
        
        for n = 1:L
            K{n} = power(cell.I{n},2).*cell.R{n}./(cell.DelTxx{n}.*AoL);
            [X ind] = sort(abs(cell.B{n}));
            
            Y = smooth(K{n}(ind),1);
            plot(X, Y,'sq','color',colz(n,:))
        end
        
        legend(leg)
        
    case 'Plot K/K0'
        
        figure;
        hold on;
        for n = 1:L
            K{n} = power(cell.I{n},2).*cell.R{n}./(cell.DelTxx{n}.*AoL);
            
            K{n} = power(cell.I{n},2).*cell.R{n}./(cell.DelTxx{n}.*AoL);
            [X ind] = sort(abs(cell.B{n}));
            
            Y = smooth(K{n}(ind),1);
            %plot(X, Y,'sq')
            
            ind = find(X == min(X));
            K0 = Y(ind(1));
            plot(X, Y/K0,'sq-','color',colz(n,:),'displayname',leg{n});
            
            if n >= 5
                ind = find(Y==min(Y));
                Hc(n) = X(ind);
                plot(Hc(n),Y(ind)/K0,'ksq','markersize',10)
            elseif n ==4
                Hc(n) = 6.692;
                plot(Hc(n),1.0302,'ksq','markersize',10)
            end
        end
        
        
        %legend(leg,'location','eastoutside'); xlim([0,18]);
        grid on; box on; xlabel('\mu_0H [T]'); ylabel('\kappa/\kappa_0')
        
%         figure; hold on;
%         plot(meanTemp,Hc)
        
    case 'Sym/AntiSym DT'
        
        figure; hold on;
        h =  -18:.01:18;
        for n = 4
            DTxx = reshape(Interp1NonUnique(cell.B{n}, cell.DelTxx{n},h),[1,length(h)]);
            DTxy = reshape(Interp1NonUnique(cell.B{n}, cell.DelTxy{n},h),[1,length(h)]);
            
            DTxxS =  (DTxx + fliplr(DTxx))/2;
            DTxxAS =  (DTxx - fliplr(DTxx))/2;
            
            DTxyS = (DTxy + fliplr(DTxy))/2;
            DTxyAS = (DTxy - fliplr(DTxy))/2;
            
            plot(h*cosd(30), DTxxS, h*cosd(30), DTxxAS, h*cosd(30), DTxyS, h*cosd(30), DTxyAS)
        end
        legend({'\DeltaT_{xx} (Symm) (slow)', '\DeltaT_{xx} (AntiSymm)', '\DeltaT_{xy} (Symm)', '\DeltaT_{xy} (AntiSymm)'})
        grid on; box on;
        xlabel('\mu_0H [T]')
        ylabel('\Delta T [K]')
        
    case 'MR Corr Temps'
        
        figure; hold on;
        
        for n = 1:L
            plot(cell.B{n}, cell.HotTemp{n},'r')
            plot(cell.B{n}, cell.ColdTemp{n},'b')
            plot(cell.B{n}, cell.TransTemp{n},'g')
            
        end
        
        
    otherwise
        display('Invlaid option')
end


end

function datacell = loadfolder(fd)
cd(fd);
D = dir; D = D(3:end);
for i = 1:length(D)
    datacell{i} = load(D(i).name);
end
end
