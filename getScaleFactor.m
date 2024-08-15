function scalefactor = getScaleFactor()

tempmin = 4;
tempmax = 80;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Load Susceptibility Data and Extract
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
intx = linspace(tempmin,tempmax,(tempmax-tempmin)*1+1);

[~,bvectorperp, ~, f0perp] = RTMdriver('res vs B','fall KYbSeSmall', 'n', '90', 'n');
[~,bvectorpara, temps, f0para] = RTMdriver('res vs B','fall KYbSeSmall', 'n', '0', 'n');
perpvec = bvectorperp.*f0perp./10^6;
%perpvec = perpvec - perpvec(length(perpvec));
paravec = -1.*bvectorpara.*f0para./10^6;
%paravec = paravec - paravec(length(paravec));
intperpvec = Interp1NonUnique(temps,perpvec,intx);
intparavec = Interp1NonUnique(temps,paravec,intx);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Load Susceptibility Data and Extract
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
axisperp = 1;
axispara = 3;
open '/home/blake/Documents/MATLAB/KYbSuscep/KYbSe2_InvChi_molOe_per_emu_12-Jul-2022.fig'
temp = gca;
perpx = temp.Children(axisperp).XData;
perpy = temp.Children(axisperp).YData;
parax = temp.Children(axispara).XData;
paray = temp.Children(axispara).YData;
intparay = Interp1NonUnique(parax,paray,intx);
intperpy = Interp1NonUnique(perpx,perpy,intx);
delchi = intparay.^(-1)-intperpy.^(-1);
%delchi = delchi - delchi(length(delchi));



perpcoeff = delchi./intperpvec;
paracoeff = delchi./intparavec;


% scalepara = mean([paracoeff(1:40-tempmin+1),paracoeff(47-tempmin+1:length(paracoeff))]);
% scaleperp = mean([perpcoeff(1:40-tempmin+1),perpcoeff(47-tempmin+1:length(perpcoeff))]);

scalepara = mean(paracoeff(1:25-tempmin+1));
scaleperp = mean(perpcoeff(1:25-tempmin+1));


scalefactor = (scalepara + scaleperp) / 2;
close(gcf)
figure;
plot(intx,delchi, 'Color', 'k')
hold on;
plot(intx,scaleperp.*intperpvec, 'Color','r', 'DisplayName', 'H||ab')
%curveFitter(intx,scalefactor.*intperpvec)
hold on;
plot(intx,scalepara.*intparavec,'Color','b', 'DisplayName', 'H||c')
%curveFitter(intx,scalefactor.*intparavec)

xlabel('Temperature (K)'); ylabel('\Delta\chi^{-1} (emu/mol*Oe)'); title(['\Delta\chi Scaled Coefficients', num2str(scalefactor)]);
x

end