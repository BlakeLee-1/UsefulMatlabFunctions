function outcell = RTMfall(plotoption)
M2 = LoadData('1.5');
M10 = LoadData ('10');
M20 = LoadData ('20');
M30 = LoadData ('30');
M40 = LoadData ('40');
M55 = LoadData ('55');
M60 = LoadData ('60');
M70 = LoadData ('70');
M100 = LoadData ('100');

colorz = lines(11);

Plot(M2, 1, colorz, plotoption)
Plot(M10, 2, flip(colorz), plotoption)
Plot(M20, 3, colorz, plotoption)
Plot(M30, 4, flip(colorz), plotoption)
Plot(M40, 5, colorz, plotoption)
Plot(M55, 6, colorz, plotoption)
Plot(M60, 7, colorz, plotoption)
Plot(M70, 8, colorz, plotoption)
Plot(M100, 9, colorz, plotoption)

end



function [outcell] = LoadData(whatTemp)
loadopt = 'Phi1';
switch loadopt
    case 'Phi1'
        fr = ['/home/blake/Documents/MATLAB/fallRTM/',whatTemp,'K/'];
        Phi_Angle = 0;
        headerlines = 8;
    case 'Phi2'
        fr = '/home/blake/Documents/MATLAB/SCM2/';
        Phi_Angle = 0;
        headerlines = 11;
    case 'S2 Phi1'
        fr = '/home/blake/Documents/MATLAB/SCM2/';
        Phi_Angle = 0;
        headerlines = 11;

    otherwise
        disp('Invalid loadopt in RTM_DC_ML_Loader!');
end

cd(fr);
DStore = dir('*.txt');
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
end
end

function [M] = Plot(M, p2,colorz, po)
% colorz = lines(length(M));
grid on;
length(M())

switch po
    case 'delf v theta'
        for i = 1:length(M)
            figure(p2)
            plot(M(i).Rotator./108, M(i).Frequency - (max(M(i).Frequency)+min(M(i).Frequency))/2, 'color', ...
                colorz(i,:), 'MarkerFaceColor', colorz(i,:), ...
                'DisplayName', [num2str(round(mean(M(i).Field))),'T   ', num2str(mean(M(i).T_Probe)),'K'], ...
                'LineWidth', 1.5)
            xlim([0,200]);
      %      ylim([-30,40]);
            xlabel('Angle (\circ)'); ylabel('\Deltaf (Hz)'); title(['CsPrSe_{2} RTM Response v Angle (',num2str(round(mean(M(i).T_Probe))),'K)']);
            hold on;
        end
    case 'f v theta'
        for i = 1:length(M)
            figure(p2)
            plot(M(i).Rotator./108, M(i).Frequency, 'color', ...
                colorz(i,:), 'MarkerFaceColor', colorz(i,:), ...
                'DisplayName', [num2str(round(mean(M(i).Field))),'T   ', num2str(mean(M(i).T_Probe)),'K'], ...
                'LineWidth', 1.5)
           
            xlim([0,200]);
            xlabel('Angle (\circ)'); ylabel('f (Hz)'); title(['CsPrSe_{2} RTM Response v Angle (',num2str(round(mean(M(i).T_Probe))),'K)']);
            hold on;
        end
    case 'temp v theta'
        for i = 1:length(M)
            figure(p2)
            plot(M(i).Rotator./108, M(i).T_Probe, 'color', ...
                colorz(i,:), 'MarkerFaceColor', colorz(i,:), ...
                'DisplayName', [num2str(round(mean(M(i).Field))),'T   ', num2str(mean(M(i).T_Probe)),'K'], ...
                'LineWidth', 1.5)
            xlim([0,200]);
            xlabel('Angle (\circ)'); ylabel('Temperature (K)'); title(['CsPrSe_{2} Angle v Temp (',num2str(round(mean(M(1).T_Probe))),'K)']);
            hold on;
        end

    case 'all delf v theta'
        hold on;
        for i = 1:length(M)
        
        figure(6)
        plot(M(i).Rotator./108, M(i).Frequency - mean(M(i).Frequency), 'color', ...
            colorz(i,:), 'MarkerFaceColor', colorz(i,:), ...
            'DisplayName', [num2str(round(mean(M(i).Field))),'T   ', num2str(mean(M(i).T_Probe)),'K'], ...
            'LineWidth', 1.5)
        xlabel('Angle (\circ)'); ylabel('\Deltaf (Hz)'); title(['CsPrSe_{2} RTM Response v Angle (',num2str(round(mean(M(1).T_Probe))),'K)']);
        hold on;
        end

    case 'res freqs harder'
        figure(6)
        coloring = lines(5);
        cutM = cutter(M);

    case 'res freqs'
        figure(6)
        coloring = lines(5);
        resonants = [];
        fields = [];
        for i = 1:length(M)
            resonants = [resonants, (max(M(i).Frequency)+min(M(i).Frequency))/2];
            fields = [fields, mean(M(i).Field)];

        end
        hold on;
        plot(fields, resonants,'o', 'color', ...
            coloring(p2,:), 'MarkerFaceColor', coloring(p2,:), ...
            'DisplayName', [num2str(mean(M(1).T_Probe)),'K'], ...
            'LineWidth', 1.5)
        xlabel('Field (T)'); ylabel('Resonant Frequency (Hz)'); title(['CsPrSe_{2} RTM Resonant Frequency v Field (',num2str(round(mean(M(i).T_Probe))),'K)']);
        hold on;
    case 'fit'
        for i = 1:length(M)
            curveFitter(M(i).Rotator.*pi/(180*108), M(i).Frequency - (max(M(i).Frequency)+min(M(i).Frequency))/2)
            mean(M(i).Field)
        end
    
    
    end

   
        
    
    grid on;

hold off;
end


function M = cutter(Min)

for i = length(Min.Fre)

end
end

