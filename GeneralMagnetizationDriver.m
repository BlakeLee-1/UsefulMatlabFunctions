function GeneralMagnetizationDriver(SN,plotoption,newfigure,interpoption, titleoption)

sampleopt = SN;
switch SN
    case 'MgCr2O4' % Just an example! Change these fields!
        fr = '/home/blake/Documents/MATLAB/MPMS Data/';
        SN = 'MgCr2O4';
        tstr = 'MgCr2O4, c-axis';
        PP.mass = 1.33e-3;
        PP.Density = 1;
        PP.MolarMass = CalculateMolarMass({'Mg','Cr','O'},[1,2,4]);
        convertMolCr = PP.mass/PP.MolarMass*2;
    otherwise
        disp('Invalid SN!');
end
BohrMagInEMU=9.27400968e-21;
PP.NumMoles=PP.mass./PP.MolarMass;
PP.NumFU=PP.NumMoles.*6.022e23;
PP.ConvertToMuB_Per_FU=1./(PP.NumFU.*BohrMagInEMU);

Mastercell_MvT=DirectorySnatchMPMSData([fr,'T/110/'],'MvsT');
Mastercell_MvH=DirectorySnatchMPMSData([fr,'H/110/'],'MvsH');

Mastercell_MvT = AddColumns_MagData(Mastercell_MvT,PP);
Mastercell_MvH = AddColumns_MagData(Mastercell_MvH,PP);

Mastercell_MvT=Sort_MagDat_Cell(Mastercell_MvT,'Field');
Mastercell_MvH=Sort_MagDat_Cell(Mastercell_MvH,'Temperature');

Mastercell_MvT = GenUColor_MvsT(Mastercell_MvT);
Mastercell_MvH = GenUColor_MvsH(Mastercell_MvH);

h2 = newfigcheck(newfigure);
switch interpoption
    case 'Interp'
        IFstring='IF_';
        ITstring='IT_';
    case 'No'
        IFstring='';
        ITstring='';
    otherwise
        disp('Default chosen...not plotting interpolated data!');
        interpstring='';
        IFstring='';
        ITstring='';
end

switch plotoption
    case 'Skip'
        disp('No plotting done, exiting...');
        
    case 'MvT All' %Plot All MvT Runs, raw data
        xarg='TemperatureK';
        yarg='LongMomentemu';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        xlabel('T [K]'); ylabel('M [emu]'); title([SN,': MvT ', titleoption]);
        
        
    case 'MvT All Error' %Plot All MvT Runs, raw data
        xarg='TemperatureK';
        yarg='LongPercentError';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        xlabel('T [K]'); ylabel('Error [%]'); title([SN,': MvT ', titleoption]);
        
    case 'M over H v T All Cutoff' %Plot M/H per molar mass. Cutoff noisy data -- over 5% error.
        cutoff=2.5;
        cutoffarg='LongPercentError';
        xarg='TemperatureK';
        yarg='Moment_over_FieldT_PerMole';
        legendarg='LegendMvT';
        plot_GiveArgs_NewCutoff_WithError(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg,cutoffarg,cutoff,SN);
        xlabel('T [K]'); ylabel('M/H [emu/(mol\cdotT)]'); title([SN,': M/H v T ', titleoption]);
        
    case 'MvH All StdDev' %Plot All MvT Runs, raw data
        xarg='FieldT';
        yarg='LongScanStdDev';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        xlabel('T [K]'); ylabel('\sigma_M [emu]'); title([SN,': MvH ', titleoption]);
        
    case 'M over H vs T' %Plot All MvT Runs, raw data
        xarg='TemperatureK';
        yarg='Moment_over_FieldT';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        xlabel('T [K]'); ylabel('M/H [emu/T]'); title([SN,': M/H v T ', titleoption]);
        xlim([0 302]);

    case 'M over H vs T per mol Cr' %Plot All MvT Runs, raw data
        xarg='TemperatureK';
        yarg='Moment_over_FieldT_PMCr';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        xlabel('T [K]'); ylabel('M/H [emu/mol_{Cr}\cdotT]'); title([SN,': M/H v T ', titleoption]);
        xlim([0 302]);
        
    case 'MvT All Per FU' %Plot All MvT Runs, raw data
        xarg='TemperatureK';
        yarg='Moment_MuB_Per_FU';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        xlabel('T [K]'); ylabel('M [\mu_B/FU]'); title([SN,': MvT ', titleoption]);
        
    case 'ChivT All Per FU' %Plot All MvT Runs, raw data
        xarg='TemperatureK';
        yarg='Chi_MuB_Per_FU';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        xlabel('T [K]'); ylabel('\chi [\mu_B/(FU*T)]'); title([SN,': \chi v T ', titleoption]);
        xlim([0 305]);
        
    case 'ChivT Volume Units' %Plot All MvT Runs, raw data
        xarg='TemperatureK';
        yarg='Chi_Volume';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        xlabel('T [K]'); ylabel('\chi [emu/(cm^3*Oe)]'); title([SN,': \chi v T ', titleoption]);
        
        
    case 'ChivT Dimensionless' %Plot All MvT Runs, raw data
        xarg='TemperatureK';
        yarg='Chi_Dimensionless';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        xlabel('T [K]'); ylabel('\chi [Dimensionless]'); title([SN,': \chi v T ', titleoption]);
        
    case 'ChivT Mass' %Plot All MvT Runs, raw data
        xarg='TemperatureK';
        yarg='Chi_Mass';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        xlabel('T [K]'); ylabel('\chi_{Mass} [emu/(g*Oe)]'); title([SN,': \chi_{Mass} v T ', titleoption]);
        
    case 'ChivT Molar' %Plot All MvT Runs, raw data
        xarg='TemperatureK';
        yarg='Chi_Molar';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        xlabel('T [K]'); ylabel('\chi_{Mass} [emu/(mol*Oe)]'); title([SN,': \chi_{Molar} v T ', titleoption]);
        
        
        
    case 'ChivT All Per FU Cutoff' %Plot All MvT Runs, raw data
        cutoff=2.5;
        cutoffarg='LongPercentError';
        xarg='TemperatureK';
        yarg='Chi_MuB_Per_FU';
        legendarg='LegendMvT';
        plot_GiveArgs_Cutoff(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg,cutoffarg,cutoff);
        xlabel('T [K]'); ylabel('\chi [\mu_B/FU*T]'); title([SN,': \chi v T ', titleoption]);
        
    case 'ChivT All Per FU Cutoff 1T' %Plot All MvT Runs, raw data
        cutoff=2.5;
        cutoffarg='LongPercentError';
        xarg='TemperatureK';
        yarg='Chi_MuB_Per_FU';
        legendarg='LegendMvT';
        plot_GiveArgs_Cutoff_1T(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg,cutoffarg,cutoff,ProcessParameters.SN);
        xlabel('T [K]'); ylabel('\chi [\mu_B/FU*T]'); title([SN,': \chi v T ', titleoption]);
        
        
    case 'MvT Specific' %Plot Specific MvT Runs
        i=input('What run would you like to see? ');
        xarg='TemperatureK';
        yarg='LongMomentemu';
        legendarg='LegendMvT';
        plot_GiveArgs({Mastercell_MvT{i}},[ITstring,xarg],[ITstring,yarg],legendarg);
        title(Mastercell_MvT{i}.LegendMvT);
        
        
    case 'MvH All' %Plot All MvH Runs
        xarg='FieldT';
        yarg='LongMomentemu';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        
    case 'MvH All MassNormal' %Plot All MvH Runs
        xarg='FieldT';
        yarg='MassNormal_Moment';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        title([SN,': M v T ', titleoption]);
        xlabel('\mu_0H [T]'); ylabel('M [emu/g]');
        xlim([-7.05 7.05])

    case 'MvH All MolNormal' %Plot All MvH Runs
        xarg='FieldT';
        yarg='MassNormal_Moment';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        title([SN,': M v T ', titleoption]);
        xlabel('\mu_0H [T]'); ylabel('M [emu/g]');
        xlim([-7.05 7.05])
        
    case 'MvH All MuBPerFU' %Plot All MvH Runs
        xarg='FieldT';
        yarg='Moment_MuB_Per_FU';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        title([SN,': M v T ', titleoption]);
        xlabel('\mu_0H [T]'); ylabel('M [\mu_B/FU]');
        xlim([-7.05 7.05])
        
    case 'dM/dH Manual' %Plot All MvH Runs
        %         xarg='FieldT';
        %         yarg='LongMomentemu';
        %         legendarg='LegendMvH';
        %         plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        for i=1:length(Mastercell_MvH)
            plot(Mastercell_MvH{i}.FieldT(1:end-1),diff(Mastercell_MvH{i}.LongMomentemu)./...
                diff(Mastercell_MvH{i}.FieldT),'-o','Color',Mastercell_MvH{i}.UColor,...
                'DisplayName',Mastercell_MvH{i}.LegendMvH)
        end
        title(tstr);
        xlabel('\mu_0H [T]'); ylabel('dM/dH [emu/T]');
        
        
    case 'Sub HF Slope'
        for i=1:length(Mastercell_MvH)
            %             cutfield = 2;
            %             kinds = Mastercell_MvH{i}.FieldT>cutfield;
            cutfield = -2;
            %             cutfield = -0.8;
            kinds = and(Mastercell_MvH{i}.FieldT<=cutfield,Mastercell_MvH{i}.FieldT>-7);
            %             kinds = and(Mastercell_MvH{i}.FieldT<=cutfield,Mastercell_MvH{i}.FieldT>-1.5);
            x = Mastercell_MvH{i}.FieldT(kinds);
            y = Mastercell_MvH{i}.LongMomentemu(kinds);
            [xData, yData] = prepareCurveData( x, y );
            ft = fittype( 'poly1' );
            opts = fitoptions( 'Method', 'LinearLeastSquares' );
            [fitresult, ~] = fit( xData, yData, ft, opts );
            cvals = coeffvalues(fitresult);
            Mastercell_MvH{i}.LegendMvH
            cvals(1)
            hold on; plot(Mastercell_MvH{i}.FieldT,Mastercell_MvH{i}.LongMomentemu-cvals(1).*Mastercell_MvH{i}.FieldT,...
                '-o','Color',Mastercell_MvH{i}.UColor,...
                'DisplayName',Mastercell_MvH{i}.LegendMvH);
        end
        xlabel('\mu_0H [T]'); ylabel('M-HF Slope [emu]');
    case 'Sub HF Slope Mass Normalized'
        for i=1:length(Mastercell_MvH)
            %             cutfield = 2;
            %             kinds = Mastercell_MvH{i}.FieldT>cutfield;
            cutfield = -2;
            kinds = and(Mastercell_MvH{i}.FieldT<=cutfield,Mastercell_MvH{i}.FieldT>-6.5);
            x = Mastercell_MvH{i}.FieldT(kinds);
            y = Mastercell_MvH{i}.MassNormal_Moment(kinds);
            [xData, yData] = prepareCurveData( x, y );
            ft = fittype( 'poly1' );
            opts = fitoptions( 'Method', 'LinearLeastSquares' );
            [fitresult, ~] = fit( xData, yData, ft, opts );
            cvals = coeffvalues(fitresult);
            Mastercell_MvH{i}.LegendMvH
            cvals(1)
            hold on; plot(Mastercell_MvH{i}.FieldT,Mastercell_MvH{i}.MassNormal_Moment-cvals(1).*Mastercell_MvH{i}.FieldT,...
                '-o','Color',Mastercell_MvH{i}.UColor,...
                'DisplayName',Mastercell_MvH{i}.LegendMvH);
        end
        xlabel('\mu_0H [T]'); ylabel('M-HF Slope [emu/g]');
        xlim([-7.05 7.05]);
        title([SN,': M vs H']);
        
    case 'Sub HF Slope Mass Normalized Data and Fits'
        for i=1:length(Mastercell_MvH)
            cutfield = 1.8;
            kinds = Mastercell_MvH{i}.FieldT>cutfield;
            x = Mastercell_MvH{i}.FieldT(kinds);
            y = Mastercell_MvH{i}.MassNormal_Moment(kinds);
            [xData, yData] = prepareCurveData( x, y );
            ft = fittype( 'poly1' );
            opts = fitoptions( 'Method', 'LinearLeastSquares' );
            [fitresult, ~] = fit( xData, yData, ft, opts );
            cvals = coeffvalues(fitresult);
            %             Mastercell_MvH{i}.LegendMvH;
            %             cvals(1);
            %             hold on; plot(Mastercell_MvH{i}.FieldT,Mastercell_MvH{i}.MassNormal_Moment-cvals(1).*Mastercell_MvH{i}.FieldT,...
            %                 '-o','Color',Mastercell_MvH{i}.UColor,...
            %                 'DisplayName',Mastercell_MvH{i}.LegendMvH);
            hold on; plot(Mastercell_MvH{i}.FieldT,Mastercell_MvH{i}.MassNormal_Moment,...
                '-o','Color',Mastercell_MvH{i}.UColor,...
                'DisplayName',Mastercell_MvH{i}.LegendMvH);
            hold on; plot(Mastercell_MvH{i}.FieldT,cvals(1).*Mastercell_MvH{i}.FieldT,...
                '--k','DisplayName',[Mastercell_MvH{i}.LegendMvH,' Fit']);
        end
        xlabel('\mu_0H [T]'); ylabel('M [emu/g]');
        xlim([-7.05 7.05]);
        title([SN,': M vs H']);
        legend show;
        
    case 'MoverHvH All' %Plot All MvH Runs
        xarg='FieldT';
        yarg='Moment_over_FieldT';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        title(tstr);
        xlabel('\mu_0H [T]'); ylabel('M/H [emu/T]');
        
        
    case 'MoverHvH All MassNormalized' %Plot All MvH Runs
        xarg='FieldT';
        yarg='MassNormal_Moment_over_FieldT';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        title(tstr);
        xlabel('\mu_0H [T]'); ylabel('M/H [emu/(g*T)]');
        xlim([-7.05 7.05]);

    case 'dMdHvH All MassNormalized' %Plot All MvH Runs
        xarg='FieldT';
        yarg='MassNormal_Moment_over_FieldT';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        title(tstr);
        xlabel('\mu_0H [T]'); ylabel('M/H [emu/(g*T)]');
        xlim([-7.05 7.05]);
    case 'MvH All Error' %Plot All MvH Runs
        xarg='FieldT';
        yarg='LongPercentError';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        title(tstr);
        xlabel('\mu_0H [T]'); ylabel('Moment Error [%]');
        xlim([-7.05 7.05]);
        
    case 'MoverHvH All MuBPerFU' %Plot All MvH Runs
        xarg='FieldT';
        yarg='Moment_MuB_Per_FU_over_FieldT';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        title(tstr);
        xlabel('\mu_0H [T]'); ylabel('M/H [\mu_B/(FU*T)]');
        xlim([-7.05 7.05]);
        
    case 'MoverHvH All Dimensionless' %Plot All MvH Runs
        xarg='FieldT';
        yarg='Chi_Dimensionless';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        title(tstr);
        xlabel('\mu_0H [T]'); ylabel('M/H [Dimensionless]');
        xlim([-7.05 7.05]);
        
        
    case 'MoverH v 1/H All' %Plot All MvH Runs
        xarg='InverseFieldT';
        yarg='Moment_over_FieldT';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        title(tstr);
        xlabel('1/\mu_0H [T^{-1}]'); ylabel('M/H [emu/T]');
        
    case 'MvH Specific' %Plot Specific MvH Runs
        i=input('What run would you like to see? ');
        xarg='FieldT';
        yarg='LongMomentemu';
        legendarg='LegendMvH';
        plot_GiveArgs({Mastercell_MvH{i}},[IFstring,xarg],[IFstring,yarg],legendarg);
        title(Mastercell_MvH{i}.LegendMvH);
        
    case 'Chi v T' %Plot All ChivT
        xarg='TemperatureK';
        yarg='Chi';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        
    case 'dMdH All' %Plot All dM/dH vs H
        xarg='FieldT';
        yarg='dMdH';
        legendarg='LegendMvH';
        plot_GiveArgs(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],legendarg);
        xlabel('\mu_0H [T]'); ylabel('M/H [emu/T]');

        
    case 'dMdH All Smooth' %Plot All dM/dH vs H Runs with smoothing
        xarg='FieldT';
        yarg='dMdH';
        legendarg='LegendMvH';
        smoothwindow=5;
        plot_GiveArgs_smooth(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],...
            legendarg,smoothwindow);
        title(tstr);
        xlabel('\mu_0H [T]'); ylabel('dM/dH [emu/T]');

    case 'dMdH fit' %Plot All dM/dH vs H Runs with smoothing
    
        xarg='FieldT';
        yarg='FitdMdH';
        legendarg='LegendMvH';
        [IFstring,yarg]
        smoothwindow=5;
        figure(10)
        plot_GiveArgs_smooth(Mastercell_MvH,[IFstring,xarg],[IFstring,yarg],...
           legendarg,smoothwindow);
        title(tstr);
        xlabel('\mu_0H [T]'); ylabel('dM/dH [emu/T]');
        
    case 'dMdT All' %Plot All dM/dT vs T
        xarg='TemperatureK';
        yarg='dMdT';
        legendarg='LegendMvT';
        plot_GiveArgs(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],legendarg);
        
    case 'dMdT All Smooth' %Plot All dM/dT vs T Runs with smoothing
        xarg='TemperatureK';
        yarg='dMdT';
        legendarg='LegendMvT';
        smoothwindow=5;
        plot_GiveArgs_smooth(Mastercell_MvT,[ITstring,xarg],[ITstring,yarg],...
            legendarg,smoothwindow);
        
    case 'Fit M vs H Mass Normalized'
        
        for i=1:length(Mastercell_MvH)
            x2 = Mastercell_MvH{i}.FieldT;
            y2 = Mastercell_MvH{i}.MassNormal_Moment;
            
            [xData, yData] = prepareCurveData( x2, y2 );
            
            % Set up fittype and options.
            ft = fittype( 'a.*(x-x0)', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.StartPoint = [0.849129305868777 0.933993247757551];
            
            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            Mastercell_MvH{i}.LegendMvH
            coeffvalues(fitresult)
        end
        
    case 'Fit M/H vs H Mass Normalized'
        
        for i=1:length(Mastercell_MvH)
            x2 = Mastercell_MvH{i}.FieldT;
            y2 = Mastercell_MvH{i}.MassNormal_Moment./x2;
            kinds = abs(x2) > 0.1;
            x = x2(kinds);
            y = y2(kinds);
            
            [xData, yData] = prepareCurveData( x, y);
            
            % Set up fittype and options.
            ft = fittype( 'a.*(abs(x)-x0)./abs(x)', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.StartPoint = [0.849129305868777 0.933993247757551];
            
            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            Mastercell_MvH{i}.LegendMvH
            coeffvalues(fitresult)
            gof.rsquare
        end
    case 'Fit M/H vs H Mass Normalized variable exponent'
        
        for i=1:length(Mastercell_MvH)
            x2 = Mastercell_MvH{i}.FieldT;
            y2 = Mastercell_MvH{i}.MassNormal_Moment./x2;
            kinds = abs(x2) > 0.4;
            x = x2(kinds);
            y = y2(kinds);
            
            [xData, yData] = prepareCurveData( x, y);
            
            % Set up fittype and options.
            ft = fittype( 'a.*(abs(x)-x0)./abs(x).^b', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.StartPoint = [0.849129305868777 0.933993247757551 .5390];
            
            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            Mastercell_MvH{i}.LegendMvH
            coeffvalues(fitresult)
            gof.rsquare
        end
        
end

end


% Internal Processing Functions
function outcell=AddColumns_MagData(incell,ProcessParameters)
PP.mass = 1.33e-3;
PP.Density = 1;
PP.MolarMass = CalculateMolarMass({'Mg','Cr','O'},[1,2,4]);
convertMolCr = PP.mass/PP.MolarMass*2;
molmass = ProcessParameters.MolarMass;
outcell=incell;
interpdensity=1;
SampleDensity = ProcessParameters.Density; % g per cm^3
for i=1:length(outcell)
    outcell{i}.FieldT=outcell{i}.FieldOe/10000;
    outcell{i}.InverseFieldT = 1./outcell{i}.FieldT;
    outcell{i}.InterpField=linspace(min(outcell{i}.FieldT),...
        max(outcell{i}.FieldT),interpdensity*length(outcell{i}.FieldT));
    outcell{i}.InterpTemperature=linspace(min(outcell{i}.TemperatureK),...
        max(outcell{i}.TemperatureK),interpdensity*length(outcell{i}.TemperatureK));
    outcell{i}.Chi=outcell{i}.LongMomentemu./outcell{i}.FieldT;
    outcell{i}.ChiUnits='emu/T';
    outcell{i}.Chi_Mass=outcell{i}.LongMomentemu./(outcell{i}.FieldOe.*ProcessParameters.mass);
    outcell{i}.Chi_Mass_Units = 'emu/(g*Oe)';
    outcell{i}.Chi_Volume = outcell{i}.Chi_Mass.*SampleDensity;
    outcell{i}.Chi_Volume_Units = 'emu/(cm^3*Oe)';
    outcell{i}.Chi_Dimensionless = outcell{i}.Chi_Volume.*(4*pi);
    outcell{i}.Chi_Molar_Units = 'emu/mol*Oe';
    outcell{i}.Chi_Molar = outcell{i}.Chi_Mass.*molmass;
    outcell{i}.MassNormal_Moment=outcell{i}.LongMomentemu./ProcessParameters.mass;
    outcell{i}.Moment_over_FieldT = outcell{i}.LongMomentemu./outcell{i}.FieldT;
    outcell{i}.MassNormal_Moment_over_FieldT = outcell{i}.MassNormal_Moment./outcell{i}.FieldT;
    
    outcell{i}.Moment_over_FieldT_PerMole = outcell{i}.LongMomentemu.*molmass./...
        (outcell{i}.FieldT.*ProcessParameters.mass);
    outcell{i}.Moment_over_FieldT_PerMole_Error = outcell{i}.LongScanStdDev.*molmass./...
        (outcell{i}.FieldT.*ProcessParameters.mass);
    if ~isfield(outcell{i},'LongPercentError')
        outcell{i}.LongPercentError = abs(100.*outcell{i}.LongScanStdDev./outcell{i}.LongMomentemu);
    end
    outcell{i}.Moment_MuB_Per_FU=outcell{i}.LongMomentemu.*ProcessParameters.ConvertToMuB_Per_FU;
    outcell{i}.Moment_MuB_Per_FU_over_FieldT=outcell{i}.LongMomentemu.*ProcessParameters.ConvertToMuB_Per_FU./outcell{i}.FieldT;
    outcell{i}.Chi_MuB_Per_FU=outcell{i}.Moment_MuB_Per_FU./outcell{i}.FieldT;
    %     cutinds=outcell{i}.LongPercentError<=2.5;
    cutinds = outcell{i}.TemperatureK > 0;
    outcell{i}.Chi_MuB_cutoff=outcell{i}.Chi_MuB_Per_FU(cutinds);
    outcell{i}.TemperatureK_cutoff=outcell{i}.TemperatureK(cutinds);
    
    storefieldnames=fieldnames(outcell{i});

%     for j = length(outcell{i}.LongPercentError):-1:2
% 
%         if outcell{i}.FieldT(j) == outcell{i}.FieldT(j-1)
%             if outcell{i}.LongPercentError(j) < outcell{i}.LongPercentError(j-1)
%             %outcell{i}.LongMomentemu(j-1) =(outcell{i}.LongMomentemu(j)+outcell{i}.LongMomentemu(j-1))/2;
%             outcell{i}.LongMomentemu(j-1) = outcell{i}.LongMomentemu(j);
% 
%             end
%             outcell{i}.LongMomentemu(j) = []
%             outcell{i}.FieldT(j) = []
%             outcell{i}.TemperatureK(j) = []
%             outcell{i}.LongPercentError(j) = []
%         
% 
%         elseif outcell{i}.LongPercentError(j) > 0.5
%             outcell{i}.LongMomentemu(j) = []
%             outcell{i}.FieldT(j) = []
%             outcell{i}.TemperatureK(j) = []
%             outcell{i}.LongPercentError(j) = []
% 
%         end
%     end


    for j=1:length(storefieldnames)
        try
            outcell{i}.(['IF_',storefieldnames{j}])=Interp1NonUnique(outcell{i}.FieldT,...
                outcell{i}.(storefieldnames{j}),outcell{i}.InterpField);
        catch
            outcell{i}.(['IF_',storefieldnames{j}])=zeros(1,length(outcell{i}.InterpField));
        end
        try
            outcell{i}.(['IT_',storefieldnames{j}])=Interp1NonUnique(outcell{i}.TemperatureK,...
                outcell{i}.(storefieldnames{j}),outcell{i}.InterpTemperature);
        catch
            outcell{i}.(['IT_',storefieldnames{j}])=zeros(1,length(outcell{i}.InterpTemperature));
        end
    end
    
    outcell{i}.LegendMvH=[num2str(mean(outcell{i}.TemperatureK),'%0.f'),' K'];
    outcell{i}.LegendMvT=[num2str(mean(outcell{i}.FieldT)),' T'];
    outcell{i}.dMdH=diff(outcell{i}.LongMomentemu)./diff(outcell{i}.FieldT);
    outcell{i}.dMdT=diff(outcell{i}.LongMomentemu)./diff(outcell{i}.TemperatureK);
    outcell{i}.IF_dMdH=diff(outcell{i}.IF_LongMomentemu)./diff(outcell{i}.InterpField);
    outcell{i}.IT_dMdT=diff(outcell{i}.IT_LongMomentemu)./diff(outcell{i}.InterpTemperature);
    outcell{i}.cutinds=cutinds;
    outcell{i}.cutoff=2.5;
    outcell{i}.Moment_over_FieldT_PMCr = outcell{i}.Moment_over_FieldT./convertMolCr;
    if length(unique(outcell{i}.InterpField)) == length(outcell{i}.InterpField)
        
    [xData, yData] = prepareCurveData(outcell{i}.InterpField, outcell{i}.IF_LongMomentemu );

    % Set up fittype and options.
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline' );
    opts.SmoothingParam = 0.999;

% Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
%     figure(i)
%     plot(outcell{i}.FieldT,outcell{i}.LongMomentemu)
%     hold on;
%     length(outcell{i}.FieldT)
%     length(outcell{i}.LongPercentError)
%     yyaxis right;
%     plot(outcell{i}.FieldT,outcell{i}.LongPercentError)

    outcell{i}.IF_FitM = fitresult(outcell{i}.InterpField);
    

    outcell{i}.IF_FitdMdH = diff(outcell{i}.IF_FitM)./diff(outcell{i}.InterpField)';

    length(outcell{i}.InterpField)
    length(outcell{i}.IF_FitdMdH)
    
    end
end
end

function outcell = GenUColor_MvsT(incell)
outcell = incell;
for i=1:length(outcell)
    outcell{i}.UColor =  ColormapInterpolate_InBounds(abs(mean(outcell{i}.FieldT)),[0,7.25]);
end
end

function outcell = GenUColor_MvsH(incell)
outcell = incell;
for i=1:length(outcell)
    outcell{i}.UColor =  ColormapInterpolate_InBounds(mean(outcell{i}.TemperatureK),[0,310]);
end

end

% Sortby can be 'Temperature' or 'Field'
function sortedoutcell=Sort_MagDat_Cell(incell,sortby)
sortedcell=incell;
for i=1:length(sortedcell)
    Temperatures(i)=mean(sortedcell{i}.TemperatureK);
    FieldVals(i)=mean(sortedcell{i}.FieldT);
end
switch sortby
    case 'Temperature'
        [~,sortinds]=sort(Temperatures);
    case 'Field'
        [~,sortinds]=sort(FieldVals);
    otherwise
        disp('No sorting done! Please sort by a valid parameter!');
end
for i=1:length(sortinds)
    sortedoutcell{i}=incell{sortinds(i)};
end
end

% Plot functions.
function plot_GiveArgs(incell,arg1,arg2,legendarg)
%colorz=distinguishable_colors(length(incell));
colorz=varycolor(length(incell));
%colorz = ['r','b','g','k'];
for i=1:length(incell)
    clz = incell{i}.UColor;
    %clz = colorz(i);
    %colors(i)
    try
%         plot(incell{i}.(arg1),incell{i}.(arg2),...
%             'Color',clz,'MarkerFaceColor',clz); hold on;
        curveFitter(incell{i}.(arg1),(incell{i}.(arg2)).^(-1))
    catch
        disp(['Plot arguments not the same length, shortening x var length to y var length: Run ',num2str(i)]);
%         incell{i}.(arg1)
%         incell{i}.(arg2)
        plot(incell{i}.(arg1)(1:length(incell{i}.(arg2))),incell{i}.(arg2),...
            'Color',clz); hold on;
    end
    %legendstore{i}=[int2str(round(mean(incell{i}.TemperatureK))),'K'];
    legendstore{i}=incell{i}.(legendarg);
end
legend(legendstore{1:end});

end
function HarryPlotter(incell)
for i = length(incell)
plot(incell{i}.IF_FitdMdH)
end

end
function plot_GiveArgs_Cutoff_1T(incell,arg1,arg2,legendarg,cutoffarg,cutoff,SN)
% colorz=distinguishable_colors(length(incell));
colorz=varycolor(length(incell));
for i=length(incell)
    cutinds=incell{i}.(cutoffarg)<=cutoff;
    [markerstr,col]=MarkerOpt_MPMS(SN);
    try
        plot(incell{i}.(arg1)(cutinds),incell{i}.(arg2)(cutinds),markerstr,...
            'Color',col,'MarkerFaceColor',brightcolor(col,.5),'MarkerSize',8,...
            'Linewidth',1.5,'DisplayName',SN); hold on;
    catch
        disp(['Plot arguments not the same length, shortening x var length to y var length: Run ',num2str(i)]);
        plot(incell{i}.(arg1)(1:length(incell{i}.(arg2))),incell{i}.(arg2),...
            markerstr,'Color',col,'MarkerFaceColor',brightcolor(col,.5),'MarkerSize',8,...
            'DisplayName',SN); hold on;
    end
    %     legendstore{i}=incell{i}.(legendarg);
end


    function [markerstring,colorval]=MarkerOpt_MPMS(SN)
        %         ShapeOpts={'o','s','^','d'};
        ShapeOpts={'o','s'};
        base='';
        colore=distinguishable_colors(2);
        switch SN
            case 'NbP'
                markerstring=[base,ShapeOpts{1}];
                colorval=colore(2,:);
            case 'TaP'
                markerstring=[base,ShapeOpts{1}];
                colorval=colore(2,:);
            case 'NbSb2'
                markerstring=[base,ShapeOpts{2}];
                colorval=colore(1,:);
            case 'TaSb2'
                markerstring=[base,ShapeOpts{2}];
                colorval=colore(1,:);
        end
    end


end

function plot_GiveArgs_smooth(incell,arg1,arg2,legendarg,window)
% colorz=distinguishable_colors(length(incell));
colorz=varycolor(length(incell));
for i=1:length(incell)
    try
        plot(incell{i}.(arg1),smooth(incell{i}.(arg2)',window),...
            'Color',colorz(i,:)); hold on;
    catch
        length(incell{i}.(arg1))
        length(incell{i}.(arg2))
        disp(['Plot arguments not the same length, shortening x var length to y var length: Run ',num2str(i)]);
        plot(incell{i}.(arg1)(1:length(incell{i}.(arg2))),smooth(incell{i}.(arg2),window),...
            'Color',colorz(i,:)); hold on;
    end
    legendstore{i}=incell{i}.(legendarg);
end
legend(legendstore{1:end});

end

function plot_GiveArgs_Cutoff(incell,arg1,arg2,legendarg,cutoffarg,cutoff)
% colorz=distinguishable_colors(length(incell));
colorz=varycolor(length(incell));
for i=1:length(incell)
    cutinds=incell{i}.(cutoffarg)<=cutoff;
    try
        plot(incell{i}.(arg1)(cutinds),incell{i}.(arg2)(cutinds),'-o',...
            'Color',colorz(i,:)); hold on;
    catch
        disp(['Plot arguments not the same length, shortening x var length to y var length: Run ',num2str(i)]);
        plot(incell{i}.(arg1)(1:length(incell{i}.(arg2))),incell{i}.(arg2),...
            'Color',colorz(i,:)); hold on;
    end
    legendstore{i}=incell{i}.(legendarg);
end
legend(legendstore{1:end});

end

function h2 = newfigcheck(newfigopt)
if strcmp(newfigopt,'New')
    h2 = figure('Position',[264 282 819 668]); hold on;
else
    h2=[];
end
end

function plot_GiveArgs_NewCutoff_WithError(incell,arg1,arg2,legendarg,cutoffarg,cutoff,sample)
% colorz=distinguishable_colors(length(incell));
for i=1:length(incell)
    cutinds=incell{i}.(cutoffarg)<=cutoff;
    markerstr = 'o';
    %     col = colorz(i,:);
    col = incell{i}.UColor;
    try
        plot(incell{i}.(arg1)(cutinds),incell{i}.(arg2)(cutinds),markerstr,...
            'Color',col,'MarkerFaceColor',brightcolor(col,.5),'MarkerSize',8,...
            'Linewidth',1.5,'DisplayName',incell{i}.(legendarg)); hold on;
        %         errorbar(incell{i}.(arg1)(cutinds),incell{i}.(arg2)(cutinds),...
        %             incell{i}.([arg2,'_Error'])(cutinds));
        %         lineProps.col = col;
        %         mseb(incell{i}.(arg1)(cutinds),incell{i}.(arg2)(cutinds)',...
        %             incell{i}.([arg2,'_Error'])(cutinds)',[],1);
    catch
        disp(['Plot arguments not the same length, shortening x var length to y var length: Run ',num2str(i)]);
        plot(incell{i}.(arg1)(1:length(incell{i}.(arg2))),incell{i}.(arg2),...
            markerstr,'Color',col,'MarkerFaceColor',brightcolor(col,.5),'MarkerSize',8,...
            'DisplayName',incell{i}.(legendarg)); hold on;
    end
    %     legendstore{i}=incell{i}.(legendarg);
end
end


