function [F,bvector, temps, f0] = RTMdriver(plotoption,session, fitt, desire, evenplot)
%fallsmall
switch session
    case 'fall KYbSePhi1'
        temps = [1.5,4,10,20,25,30,35,40,45,50,55,60,70,100];
    case 'fall KYbSePhi2'
        temps = [10,20,30,40,45,50,70];
    case 'fall CsErSeS1'
        temps = [2,15,35,70,140];
    case 'fall CsErSeS2'
        temps = [2,7,15,35,45,55,70,90,110,140];
    case 'spring CsPrSe'
        temps = [2,5,20,60,120];
    case 'fall KYbSeSmall'
        temps = [1.5,4,10,20,30,40,45,50,60,80];
end
desire = str2num(desire);
[M, titleboi,fr] = LoadData(temps, session);

[F,bvector, f0] = Plot(M, plotoption, temps, titleboi,fr, fitt, desire,evenplot);
end


colorbar
function [outcell,titleboi,frt] = LoadData(whatTemp,session)
for j = 1:length(whatTemp)
loadopt = 'Phi1';
switch session
    case 'fall CsErSeS1'
        frt = '/home/blake/Documents/MATLAB/fallRTM/CsErSe2S1/';
        fr = [frt,num2str(whatTemp(j)),'K/'];
        headerlines = 8;
        titleboi = 'CsErSe_{2} Sample 1';
    case 'fall CsErSeS2'
        frt = '/home/blake/Documents/MATLAB/fallRTM/CsErSe2S2/';
        fr = [frt,num2str(whatTemp(j)),'K/'];
        headerlines = 9;
        titleboi = 'CsErSe_{2} Sample 2';
    case 'fall KYbSePhi1'
        frt = '/home/blake/Documents/MATLAB/fallRTM/KYbSePhi1/';
        fr = ['/home/blake/Documents/MATLAB/fallRTM/KYbSePhi1/',num2str(whatTemp(j)),'K/'];
        Phi_Angle = 0;
        headerlines = 8;
        titleboi = 'KYbSe_{2} Phi1';
    case 'spring CsPrSe'
        frt = '/home/blake/Documents/MATLAB/SCM2/CsPrSe2Data/';
         fr = ['/home/blake/Documents/MATLAB/SCM2/CsPrSe2Data/',num2str(whatTemp(j)),'K/'];
        Phi_Angle = 0;
        headerlines = 11;
        titleboi = 'CsPrSe_{2}';

    case 'fall KYbSePhi2'
        frt = '/home/blake/Documents/MATLAB/fallRTM/KYbSePhi2/';
        fr = ['/home/blake/Documents/MATLAB/fallRTM/KYbSePhi2/',num2str(whatTemp(j)),'K/'];
        Phi_Angle = 0;
        headerlines = 9;
        titleboi = 'KYbSe_{2} Phi2';

    case 'fall KYbSeSmall'
        frt = '/home/blake/Documents/MATLAB/fallRTM/KYbSeSmall/';
        fr = ['/home/blake/Documents/MATLAB/fallRTM/KYbSeSmall/',num2str(whatTemp(j)),'K/'];
        Phi_Angle = 0;
        headerlines = 9;
        titleboi = 'KYbSe_{2} Small';

    otherwise
        disp('Invalid loadopt in RTM_DC_ML_Loader!');
end

cd(fr);
DStore = dir('*.txt');
switch session
    case 'fall CsErSeS1'
        for i=1:length(DStore)
            temp = readmatrix(DStore(i).name);
            outcell(j,i).Field = temp(headerlines:end,1);
            outcell(j,i).Frequency = temp(headerlines:end,2);
            outcell(j,i).Gamma = temp(headerlines:end,3);
            outcell(j,i).Rotator = temp(headerlines:end,4);
            outcell(j,i).T_Probe = temp(headerlines:end,5);
            outcell(j,i).T_VTI = temp(headerlines:end,6);
            outcell(j,i).Time = temp(headerlines:end,7);
            outcell(j,i).Filename = DStore(i).name;
            outcell(j,i).meanTemp = mean(outcell(j,i).T_Probe);
            outcell(j,i).meanField = round(2*mean(outcell(j,i).Field))/2;
            outcell(j,i).rawmeanField = mean(outcell(j,i).Field);
            outcell(j,i).roundmeanTemp = round(2*outcell(j,i).meanTemp)/2;
            closest = 1;
            for k = 1:length(whatTemp)
                if abs(outcell(j,i).meanTemp - whatTemp(k)) < abs(outcell(j,i).meanTemp - whatTemp(closest))
                    closest = k;
                end
            end
            outcell(j,i).closestTemp = whatTemp(closest);

        end
    case 'fall CsErSeS2'
        for i=1:length(DStore)
            temp = readmatrix(DStore(i).name);
            outcell(j,i).Field = temp(headerlines:end,1);
            outcell(j,i).Frequency = temp(headerlines:end,2);
            outcell(j,i).Gamma = temp(headerlines:end,3);
            outcell(j,i).Rotator = temp(headerlines:end,4);
            outcell(j,i).T_Probe = temp(headerlines:end,5);
            outcell(j,i).T_VTI = temp(headerlines:end,6);
            outcell(j,i).Time = temp(headerlines:end,7);
            outcell(j,i).Filename = DStore(i).name;
            outcell(j,i).meanTemp = mean(outcell(j,i).T_Probe);
            outcell(j,i).meanField = round(2*mean(outcell(j,i).Field))/2;
            outcell(j,i).rawmeanField = mean(outcell(j,i).Field);
            outcell(j,i).roundmeanTemp = round(2*outcell(j,i).meanTemp)/2;
            closest = 1;
            for k = 1:length(whatTemp)
                if abs(outcell(j,i).meanTemp - whatTemp(k)) < abs(outcell(j,i).meanTemp - whatTemp(closest))
                    closest = k;
                end
            end
            outcell(j,i).closestTemp = whatTemp(closest);
        end
    case 'spring CsPrSe'
        for i=1:length(DStore)
            temp = readmatrix(DStore(i).name);
            outcell(j,i).Field = temp(headerlines:end,1);
            outcell(j,i).Hallx = temp(headerlines:end,2);
            outcell(j,i).Hally = temp(headerlines:end,3);
            outcell(j,i).T_VTI = temp(headerlines:end,4);
            outcell(j,i).T_Probe = temp(headerlines:end,5);
            outcell(j,i).Frequency = temp(headerlines:end,6);
            outcell(j,i).Gamma = temp(headerlines:end,7);
            outcell(j,i).Rotator = temp(headerlines:end,8);  
            outcell(j,i).Time = temp(headerlines:end,9);
            outcell(j,i).Filename = DStore(i).name;
            outcell(j,i).meanTemp = mean(outcell(j,i).T_Probe);
            outcell(j,i).meanField = round(2*mean(outcell(j,i).Field))/2;
            outcell(j,i).rawmeanField = mean(outcell(j,i).Field);
            outcell(j,i).roundmeanTemp = round(2*outcell(j,i).meanTemp)/2;
            closest = 1;
            for k = 1:length(whatTemp)
                if abs(outcell(j,i).meanTemp - whatTemp(k)) < abs(outcell(j,i).meanTemp - whatTemp(closest))
                    closest = k;
                end
            end
            outcell(j,i).closestTemp = whatTemp(closest);
        end
    case 'fall KYbSePhi1'
        for i=1:length(DStore)
            temp = readmatrix(DStore(i).name);
            outcell(j,i).Field = temp(headerlines:end,1);
            outcell(j,i).Frequency = temp(headerlines:end,2);
            outcell(j,i).Gamma = temp(headerlines:end,3);
            outcell(j,i).Rotator = temp(headerlines:end,4);
            outcell(j,i).T_Probe = temp(headerlines:end,5);
            outcell(j,i).T_VTI = temp(headerlines:end,6);
            outcell(j,i).Time = temp(headerlines:end,7);
            outcell(j,i).Filename = DStore(i).name;
            outcell(j,i).meanTemp = mean(outcell(j,i).T_Probe);
            outcell(j,i).meanField = round(2*mean(outcell(j,i).Field))/2;
            outcell(j,i).rawmeanField = mean(outcell(j,i).Field);
            outcell(j,i).roundmeanTemp = round(2*outcell(j,i).meanTemp)/2;
            closest = 1;
            for k = 1:length(whatTemp)
                if abs(outcell(j,i).meanTemp - whatTemp(k)) < abs(outcell(j,i).meanTemp - whatTemp(closest))
                    closest = k;
                end
            end
            outcell(j,i).closestTemp = whatTemp(closest);
        end
    case 'fall KYbSePhi2'
        for i=1:length(DStore)
            temp = readmatrix(DStore(i).name);
            outcell(j,i).Field = temp(headerlines:end,1);
            outcell(j,i).Frequency = temp(headerlines:end,2);
            outcell(j,i).Gamma = temp(headerlines:end,3);
            outcell(j,i).Rotator = temp(headerlines:end,4);
            outcell(j,i).T_Probe = temp(headerlines:end,5);
            outcell(j,i).T_VTI = temp(headerlines:end,6);
            outcell(j,i).Time = temp(headerlines:end,7);
            outcell(j,i).Filename = DStore(i).name;
            outcell(j,i).meanTemp = mean(outcell(j,i).T_Probe);
            outcell(j,i).meanField = round(2*mean(outcell(j,i).Field))/2;
            outcell(j,i).rawmeanField = mean(outcell(j,i).Field);
            outcell(j,i).roundmeanTemp = round(2*outcell(j,i).meanTemp)/2;
            closest = 1;
            for k = 1:length(whatTemp)
                if abs(outcell(j,i).meanTemp - whatTemp(k)) < abs(outcell(j,i).meanTemp - whatTemp(closest))
                    closest = k;
                end
            end
            outcell(j,i).closestTemp = whatTemp(closest);
            
        end
    case 'fall KYbSeSmall'
        for i=1:length(DStore)
            temp = readmatrix(DStore(i).name);
            outcell(j,i).Field = temp(headerlines:end,1);
            outcell(j,i).Frequency = temp(headerlines:end,2);
            outcell(j,i).Gamma = temp(headerlines:end,3);
            outcell(j,i).Rotator = temp(headerlines:end,4);
            outcell(j,i).T_Probe = temp(headerlines:end,5);
            outcell(j,i).T_VTI = temp(headerlines:end,6);
            outcell(j,i).Time = temp(headerlines:end,7);
            outcell(j,i).Filename = DStore(i).name;
            outcell(j,i).meanTemp = mean(outcell(j,i).T_Probe);
            outcell(j,i).meanField = round(2*mean(outcell(j,i).Field))/2;
            outcell(j,i).rawmeanField = mean(outcell(j,i).Field);
            outcell(j,i).roundmeanTemp = round(2*outcell(j,i).meanTemp)/2;
            closest = 1;
            for k = 1:length(whatTemp)
                if abs(outcell(j,i).meanTemp - whatTemp(k)) < abs(outcell(j,i).meanTemp - whatTemp(closest))
                    closest = k;
                end
            end
            outcell(j,i).closestTemp = whatTemp(closest);
        end
end
end
end

function [F,bvector, f0] = Plot(M, po,temps, titleboi, fr, fitt, desire,evenplot)
bvector = [];
f0 = [];
grid on;
            lensiz = size(M);
            len = lensiz(2);

colorz = lines(36);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Sort M to make it more manageable
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
allfields = [];
        for j = 1:length(temps)
            for i = 1:len
                if isempty(M(j,i).Rotator) == 0
                    
                    allfields = [allfields,M(j,i).meanField];
                end
            end
        end
        uniquefields = sort(unique(allfields));
        for i = length(uniquefields):-1:1
            if uniquefields(i) == 8.5
                uniquefields(i) = [];
            end
        end

       
        for j = 1:length(temps)
            for i = 1:len
                if isempty(M(j,i).Rotator) == 0
                    leng(j) = i;
                end
            end
        end
        for j = 1:length(temps)
            for i = 1:leng(j)
                for k = 1:length(temps)
                    for l = 1:length(uniquefields)
                        if (M(j,i).closestTemp) == temps(k)
                            if round(2*M(j,i).meanField)/2 == uniquefields(l)
                                Msort(k,l) = M(j,i);
                            elseif (floor(M(j,i).meanField) == 8) && uniquefields(l) == 8  
                                Msort(k,l) = M(j,i);  
                            end
                        end
                    end
                    
                end
                
            end
        end

    lowestindex = [];
    for j = 1:length(temps)
    i = 1;
        while (isempty(Msort(j,i).Field)) == 1
            i = i + 1;
        end
    lowestindex = [lowestindex,i];
    end

    lensiz = size(Msort);
    len = lensiz(2);

 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%End of Sort M
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




switch po
    
    case 'delf v theta'
        for j = 1:length(temps)

            figure;
            
            for i = lowestindex(j):len
                if isempty(Msort(j,i).Field) == 0
                hold on;
                clz = ColormapInterpolate_InBounds(1:36,[0,36]);
switch evenplot
    case 'y'
                plot(Msort(j,i).Rotator./103, Msort(j,i).Frequency - (max(Msort(j,lowestindex(j)).Frequency)+min(Msort(j,lowestindex(j)).Frequency))/2, 'color', ...
                    clz(round(2*mean(Msort(j,i).Field)),:), 'MarkerFaceColor', clz(round(2*mean(Msort(j,i).Field)),:), ...
                    'DisplayName', [num2str(0.5*round(2*mean(Msort(j,i).Field))),'T   ', num2str(temps(j)),'K'], ...
                    'LineWidth', 3)
    case 'n'
end
                end
                
                
                
            end
               % xlim([0,180]);
                xlabel('Angle (\circ)'); ylabel('\Deltaf (Hz)'); title([titleboi,' RTM Response v Angle (',num2str(temps(j)),'K)']);
                hold on;
                %                 If Saving Plots
                saveas(gcf,[fr,'plots/',num2str(temps(j)),'K.jpg'])
                saveas(gcf,[fr,'plots/',num2str(temps(j)),'K.fig'])
        end

  

   
    case 'temp v theta'
        for j = 1:length(temps)
            colorz = lines(len); 
            if M(j,1).Rotator(500) < M(j,1).Rotator(1)
                    colorz = flip(colorz);
            end
            for i = 1:len
            figure(j)
            switch evenplot
                case 'y'
            plot(M(i).Rotator./103, M(i).T_Probe, 'color', ...
                colorz(i,:), 'MarkerFaceColor', colorz(i,:), ...
                'DisplayName', [num2str(round(mean(M(i).Field))),'T   ', num2str(mean(M(i).T_Probe)),'K'], ...
                'LineWidth', 1.5)
                case 'n'
            end
            xlim([0,200]);
            xlabel('Angle (\circ)'); ylabel('Temperature (K)'); title([titleboi,' RTM Response v Angle (',num2str(temps(j)),'K)']);
            hold on;
            end
        end

    case 'res vs B'
        coloring = lines(length(temps));


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        viewAngle = findviewangle(desire, Msort, temps,lowestindex);
        %viewAngle = 7
        tol = 5;
        f0 = [];
        bvector = [];
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
        for j = 1:length(temps)
            F(j).temp = temps(j);
            resonants = [];
            temperatures = [];
            fields = [];
            colorz = lines(36); 
            f0 = [f0, (max(Msort(j,lowestindex(j)).Frequency)+min(Msort(j,lowestindex(j)).Frequency))/2];
            
            f0(j);
            
            for i = 1:len
                jj = length(Msort(j,i).Rotator);
                for k = 1:length(Msort(j,i).Rotator)
                    if abs(Msort(j,i).Rotator(k) - viewAngle*103) < abs(abs(Msort(j,i).Rotator(jj) - viewAngle*103))
                        jj = k;
                    end
                end
                lengthofrot = length(Msort(j,i).Rotator);
                lengthofrot;
                jj;


                if isempty(Msort(j,i).Rotator) == 0
                    length(Msort(j,i).Rotator);
                    f0;
                    resonants = [resonants, mean(Msort(j,i).Frequency(jj-tol:jj+tol))-(max(Msort(j,lowestindex(j)).Frequency)+min(Msort(j,lowestindex(j)).Frequency))/2];
                    fields = [fields, mean(Msort(j,i).Field)];
                    mean(Msort(j,i).Field);
                    fields;
                end
            end
            
            F(j).C_frequency = resonants;
            F(j).C_Field = fields;
            F(j).zf_frequency = f0(j);
            F(j).Temperature = temps(j);
            
            

        hold on;
        length(fields);
        length(resonants);
        
        
        clz = ColormapInterpolate_InBounds(temps,[0,max(temps)]);
        white = [1,1,1];
        F(j).UColor = clz(j,:);
        switch evenplot
            case 'y'

        if fitt == 'y'
            curveFitter(fields,resonants)
        elseif fitt == 'n'
        plot(fields, resonants,'-sq', 'Color', ...
            clz(j,:), 'MarkerFaceColor', white', ...
            'DisplayName', [num2str(round(mean(M(j,1).T_Probe))),'K'], ...
            'LineWidth', 2.5);
        xlabel('Field (T)'); ylabel('\Deltaf (Hz)'); title([titleboi,' \Deltaf v Field at ',num2str(viewAngle),' Degrees']);
        hold on;

        bvector = [bvector,findcoeff(fields,resonants)];
        
        end
        case 'n'
        bvector = [bvector,findcoeff(fields,resonants)];
        end
        end
    

    case 'res vs T'
        figure;
        allfields = [];
        for j = 1:length(temps)
            for i = 1:len
                if isempty(M(j,i).Rotator) == 0
                    
                    allfields = [allfields,M(j,i).meanField];
                end
            end
        end
        uniquefields = sort(unique(allfields));
        for i = length(uniquefields):-1:1
            if uniquefields(i) == 8.5
                uniquefields(i) = [];
            end
        end

       
        for j = 1:length(temps)
            for i = 1:len
                if isempty(M(j,i).Rotator) == 0
                    leng(j) = i;
                end
            end
        end
        
        for j = 1:length(temps)
            for i = 1:leng(j)
                for k = 1:length(temps)
                    for l = 1:length(uniquefields)
                        if round(M(j,i).meanTemp) == temps(k)
                            if round(2*M(j,i).meanField)/2 == uniquefields(l)
                                Msort(k,l) = M(j,i);
                            elseif floor(M(j,i).meanField) == 8
                                if uniquefields(l) == 8
                                    Msort(k,l) = M(j,i);
                                end
                            end
                        end
                    end
                    
                end
                
            end
        end

        for i = 1:len
            resonants = [];
            temperatures = [];
            fields = [];
            colorz = lines(36); 
            
            
            for j = 1:length(temps)

                if isempty(Msort(j,i).Rotator) == 0
                    resonants = [resonants, (max(Msort(j,i).Frequency)+min(Msort(j,i).Frequency))/2];
                    temperatures = [temperatures, mean(Msort(j,i).T_Probe)];
                    fields = [fields, Msort(j,i).meanField];
                end
                avgF = mean(fields);
    
            end
        
        hold on;
switch evenplot
    case 'y'
        plot(temperatures, resonants,'o', 'color', ...
            colorz(round(2*mean(M(j,i).Field)),:), 'MarkerFaceColor', colorz(round(2*mean(M(j,i).Field)),:), ...
            'DisplayName', [num2str(round(2*avgF)/2),'T'], ...
            'LineWidth', 1.5)
    case 'n'
end

        xlabel('Temperature (K)'); ylabel('Resonant Frequency (Hz)'); title([titleboi,' RTM Response v Temp']);
        hold on;
        
        end



    
end

   
        
    
    grid on;

hold off;
end

function viewangle = findviewangle(desired,M,temps,lowfield)
viewangle = 0;
shifts = [];
verticals = [];
for i = 1:2
    for j = 1:length(temps)
        if isempty(M(j,max(lowfield)-1+i).Rotator) == 0
            x = M(j,max(lowfield)-1+i).Rotator./103;
            y = M(j,max(lowfield)-1+i).Frequency - (max(M(j,max(lowfield)-1+i).Frequency) + min(M(j,max(lowfield)-1+i).Frequency))/2;
            % Define Start points, fit-function and fit curve
            x0 = [1 1 1 1 1]; 
            fitfunc = fittype( @(a,b,c,d,f, x) a*cosd(2*(x+b))+c + d*cosd(4*(x+f)));
            [fitted_curve,~] = fit(x,y,fitfunc,'StartPoint',x0);
            % Save the coeffiecient values for a,b,c and d in a vector
            coeffvals = coeffvalues(fitted_curve);
%            Plot results
%             figure;
%             scatter(x, y, 'r+')
%             hold on
%             plot(x,fitted_curve(x))
%             title(temps(j))

            hold off;
            shifts = [shifts,coeffvals(2)];
            verticals = [verticals, coeffvals(3)];
            
        end
    end
end
    shifts2keep = [shifts(2:4),shifts(length(temps)+2:length(temps)+4)];
    viewangle = mean(abs(shifts2keep))+desired;

end

function bvectorcomp = findcoeff(fields,resonants)
            whereis6 = 1;
            for i = 1:length(fields)
                if fields(i) < 7
                    whereis6 = i;
                end
            end
            x = fields(1:whereis6)';
            y = resonants(1:whereis6)';         % Define Start points, fit-function and fit curve
            x0 = [1 1]; 
            fitfunc = fittype( @(a,b,x) a+b*x.^2);
            [fitted_curve,~] = fit(x,y,fitfunc,'StartPoint',x0);
            % Save the coeffiecient values for a,b,c and d in a vector
            coeffvals = coeffvalues(fitted_curve);
%            Plot results
%             figure;
%             scatter(x, y, 'r+')
%             hold on
%             plot(x,fitted_curve(x))
%             title(temps(j))

            hold off;
            bvectorcomp = coeffvals(2);
            


end
