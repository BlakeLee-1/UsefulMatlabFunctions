% Ian Leahy
% 1/25/2022
% CEF Magnetotropic Modeling User Friendly Rewrite
% ======================================================================
% =============================Inputs===================================
% ======================================================================
%  SampleOption (string): Currently configured for two samples
% 'CsYbSe2' and 'KYbSe2', new sample options can be added in the
% 'LoadSampleData' internal function.
% ======================================================================
% AnalysisOption (string): Current options are 'FitChi, Magnetotropic',
% 'Fit Chi, Magnetotropic, with Bound Eigenenergies', 'Display Fit Results
% k', and 'Find Scaling Amplitude'.
% ======================================================================
% InputParameters (1 x 9 vector): Input vector for parameters. Should be
% shaped as [Bmn.*FAC Jxx Jzz Chi0], where Bmn are the steven's operator
% coefficients in meV, FAC is given as [1, 10, 1, 1e3, 1e2, 1e2], Jxx and
% Jzz are calculated from a low temperature Curie-Weiss fit and Chi0 is an
% optional diamagnetic offset in the susceptibility (often not necessary
% and set to zero). If the Analysis Option does not require a set of Input
% Parameters, the InputParameters can be set to a string with any contents
% and it will be replaced with a vector of zeros.
% ======================================================================
% ==============================Outputs=================================
% ======================================================================
% Outcell -- Contingent on AnalysisOption, usually returns fit parameter
% structure, if not, empty.
%
% Output fitting Parameters have FAC within them.
%
% =============================ChangeLog================================
% v0  - 1/25/2022 -- Only includes option 'Fit Chi, Magnetotropic' and
% 'Display Fit Results'. Created and tested all internal functions
% necessary for the calculations of k. This method does not use the
% self-consistent mean field theory for calculation.
%     -  2/4/2022 -- Added a new cumulative error function
%     (CumulativeError_BoundEigenvalues) that properly weights the
%     residuals for the magnetotropic and chi data and offers the
%     additional functionality of bounding the fit result to have
%     eigenenergies that match a known set of eigenenergies. These
%     bounding energies are set in 'LoadSampleData' as a ProcessParameters
%     entry. Added new ProcessParameters functionality, now contains fields
%     that can bias fitting to prioritize energy eigenvalues, chi, or
%     magnetotropic data. Added an AnalysisOption, 'Find Scaling Amplitude'
%     that helps in calculating the scale factor necessary to scale the
%     magnetotropic data into real units. Once found, the amplitude can be
%     set as Amp in ProcessParameters.
%     - 2/8/2022 -- Added an AnalysisOption,'Find B20 from Susceptibility'
%     for calculating the B20 Steven's operator coefficient from
%     susceptibility data.

function [outcell] = CEF_MagnetotropicFittingAndCalculation(SampleOption,AnalysisOption,InputParameters)

% Run internal functions to generate matrices and constants for
% calculations.
[ID,Jx,~,Jz,~,~] = GenerateSpinMatrices();
[muB,kB,gJ,A1,C0,~,~] = GenerateConstants();
[MagnetotropicData_AB,MagnetotropicData_C,SusceptibilityData,ProcessParameters]=LoadSampleData(SampleOption);
[MagnetotropicData_AB]=AddColumnsMagnetotropicData(MagnetotropicData_AB,ProcessParameters);
[MagnetotropicData_C]=AddColumnsMagnetotropicData(MagnetotropicData_C,ProcessParameters);
[Bmn_Parameters,Jxx_Parameter,Jzz_Parameter,chi0_Parameter] = GenerateParameters(InputParameters);
[CEF_Hamiltonian] = GenerateCEFHamiltonian();


switch AnalysisOption
    case 'Fit Chi, Magnetotropic'

        % First, plot the susceptibility and magnetotropic data to be fit.
        h2 = figure;
        [axlist,~] = tight_subplot(1,3,[.01,.04],[.10 .10],[0.075 .075]);

        % Plot susceptibility in subplot 1.
        figure(h2);
        axes(axlist(1)); hold on;
        plot(SusceptibilityData.Temperature,SusceptibilityData.Inverse_ChiAB,...
            'bsq','MarkerFaceColor',brc([0 0 1],.5),'DisplayName','H||ab');
        plot(SusceptibilityData.Temperature,SusceptibilityData.Inverse_ChiC,...
            'rsq','MarkerFaceColor',brc([1 0 0],.5),'DisplayName','H||c');
        xlabel('T [K]'); ylabel('\chi^{-1} [mol Oe/emu]');
        axes(axlist(2)); hold on;
        axes(axlist(3)); hold on;
        for q=1:length(MagnetotropicData_AB)
            % Plot AB data in subplot 2.
            plot(axlist(2),MagnetotropicData_AB(q).InterpFields,...
                MagnetotropicData_AB(q).k_interp,'sq',...
                'Color', MagnetotropicData_AB(q).UColor,...
                'MarkerFaceColor',brc(MagnetotropicData_AB(q).UColor,.5),...
                'DisplayName',num2str(MagnetotropicData_AB(q).Temperature));
            % Plot C data in subplot 3.
            plot(axlist(3),MagnetotropicData_C(q).InterpFields,...
                MagnetotropicData_C(q).k_interp,'sq',...
                'Color', MagnetotropicData_C(q).UColor,...
                'MarkerFaceColor',brc(MagnetotropicData_C(q).UColor,.5),...
                'DisplayName',num2str(MagnetotropicData_C(q).Temperature));
        end
        axes(axlist(2)); xlabel('\mu_0H [T]'); ylabel('k_{ab}');
        axes(axlist(3)); xlabel('\mu_0H [T]'); ylabel('k_c');


        % Define the cumulative error function. First, need to define all
        % relevant inputs to the Cumulative Error function.
        MinimizationFunction = @(Bmn) CumulativeError(Bmn,...
            Jxx_Parameter,Jzz_Parameter,chi0_Parameter,...
            CEF_Hamiltonian,ID,Jx,Jz,{MagnetotropicData_AB,MagnetotropicData_C},...
            SusceptibilityData,muB,kB,gJ,A1,C0,ProcessParameters);

        % fminsearchbnd needs the minimization function, an initial guess
        % at the parameters, lower bounds, and upper bounds.
        InitialGuess = Bmn_Parameters;
        LowerBounds = -1.*ones(1,length(Bmn_Parameters));
        UpperBounds = ones(1,length(Bmn_Parameters));
        [FitParameters,EndFunctionValue] = fminsearchbnd(MinimizationFunction,...
            InitialGuess,LowerBounds,UpperBounds);

        % Display the results.
        display(FitParameters);
        display(EndFunctionValue);

        % Calculate the relevant magnetotropic traces at different
        % temperatures and plot them on top of the data. Need to first
        % calculate the Hamiltonian as it is used to calculate the k.
        CEF_Hamiltonian_MatrixFit = CEF_Hamiltonian(FitParameters);
        [~,EigenValues_Fit] = eig(CEF_Hamiltonian_MatrixFit);
        [EnergyValues_Fit,~] = sort(diag(EigenValues_Fit));
        ZeroEnergy_Fit = EnergyValues_Fit(1);

        % We will plot the resulting fit at fields between 0 and 60T.
        InterpField = 0:0.1:20;
        LengthFitRange = length(ProcessParameters.FitRange) ;
        FitRange = ProcessParameters.FitRange;
        for q=1:LengthFitRange
            k = FitRange(q);
            CurrentTemperature = MagnetotropicData_AB(k).Temperature;

            % Plot the fits.
            kab_Calc = Magnetotropic_AB_vs_H(CEF_Hamiltonian_MatrixFit,...
                ZeroEnergy_Fit,ID,Jx,Jz,Jxx_Parameter,Jzz_Parameter,...
                CurrentTemperature,InterpField,muB,kB,gJ,A1);
            plot(axlist(2),InterpField,kab_Calc,'-','Color',...
                brc([0 0 0],.5),'DisplayName',[num2str(CurrentTemperature),' K Fit']);

            kc_Calc = Magnetotropic_C_vs_H(CEF_Hamiltonian_MatrixFit,...
                ZeroEnergy_Fit,ID,Jx,Jz,Jxx_Parameter,Jzz_Parameter,...
                CurrentTemperature,InterpField,muB,kB,gJ,A1);
            plot(axlist(3),InterpField,kc_Calc,'-','Color',...
                brc([0 0 0],.5),'DisplayName',[num2str(CurrentTemperature),' K Fit']);
        end

        % Plot the susceptibility as well.
        axes(axlist(1)); hold on;
        TemperaturesPlot = 1:300;
        Inverse_ChiAB_Calc = InverseSusceptibility_AB_vs_T(CEF_Hamiltonian_MatrixFit,...
            ZeroEnergy_Fit,ID,Jx,Jxx_Parameter,TemperaturesPlot,muB,kB,gJ,C0);
        Inverse_ChiC_Calc = InverseSusceptibility_C_vs_T(CEF_Hamiltonian_MatrixFit,...
            ZeroEnergy_Fit,ID,Jz,Jzz_Parameter,TemperaturesPlot,muB,kB,gJ,C0);
        plot(TemperaturesPlot,Inverse_ChiAB_Calc,'b');
        plot(TemperaturesPlot,Inverse_ChiC_Calc,'r');
    case 'Fit Chi, Magnetotropic, with Bound Eigenenergies'

        % First, plot the susceptibility and magnetotropic data to be fit.
        h2 = figure;
        [axlist,~] = tight_subplot(1,3,[.01,.04],[.10 .10],[0.075 .075]);

        % Plot susceptibility in subplot 1.
        figure(h2);
        axes(axlist(1)); hold on;
        plot(SusceptibilityData.Temperature,SusceptibilityData.Inverse_ChiAB,...
            'bsq','MarkerFaceColor',brc([0 0 1],.5),'DisplayName','H||ab');
        plot(SusceptibilityData.Temperature,SusceptibilityData.Inverse_ChiC,...
            'rsq','MarkerFaceColor',brc([1 0 0],.5),'DisplayName','H||c');
        xlabel('T [K]'); ylabel('\chi^{-1} [mol Oe/emu]');
        set(gca,'XTickLabelMode','auto'); set(gca,'YTickLabelMode','auto');
        axes(axlist(2)); hold on;
        set(gca,'XTickLabelMode','auto'); set(gca,'YTickLabelMode','auto');
        axes(axlist(3)); hold on;
        set(gca,'XTickLabelMode','auto'); set(gca,'YTickLabelMode','auto');

        for q=ProcessParameters.FitRange
            % Plot AB data in subplot 2.
            plot(axlist(2),MagnetotropicData_AB(q).InterpFields,...
                MagnetotropicData_AB(q).k_interp,'sq',...
                'Color', MagnetotropicData_AB(q).UColor,...
                'MarkerFaceColor',brc(MagnetotropicData_AB(q).UColor,.5),...
                'DisplayName',num2str(MagnetotropicData_AB(q).Temperature));
            % Plot C data in subplot 3.
            plot(axlist(3),MagnetotropicData_C(q).InterpFields,...
                MagnetotropicData_C(q).k_interp,'sq',...
                'Color', MagnetotropicData_C(q).UColor,...
                'MarkerFaceColor',brc(MagnetotropicData_C(q).UColor,.5),...
                'DisplayName',num2str(MagnetotropicData_C(q).Temperature));
        end
        axes(axlist(2)); xlabel('\mu_0H [T]'); ylabel('k_{ab}');
        axes(axlist(3)); xlabel('\mu_0H [T]'); ylabel('k_c');


        % Define the cumulative error function. First, need to define all
        % relevant inputs to the Cumulative Error function.
        MinimizationFunction = @(Bmn) CumulativeError_BoundEigenvalues(Bmn,...
            Jxx_Parameter,Jzz_Parameter,chi0_Parameter,...
            CEF_Hamiltonian,ID,Jx,Jz,{MagnetotropicData_AB,MagnetotropicData_C},...
            SusceptibilityData,muB,kB,gJ,A1,C0,ProcessParameters);

        % fminsearchbnd needs the minimization function, an initial guess
        % at the parameters, lower bounds, and upper bounds. Define
        % FitStruct, which holds information about the fitting routine.
        InitialGuess = Bmn_Parameters;


         LowerBounds = -100.*ones(1,length(Bmn_Parameters));
         UpperBounds = 100.*ones(1,length(Bmn_Parameters));



%         LowerBounds = [-100   -0.7747    0.5631   -100   -0.3641   -1.7954];
%         UpperBounds = [100   -0.7747    0.5631   100   -0.3641   -1.7954];
%                 LowerBounds = [-6.2,-100.*ones(1,length(Bmn_Parameters)-1)];
%                 UpperBounds = [-5.5,100.*ones(1,length(Bmn_Parameters)-1)];

        % %         fminoptions = optimset('PlotFcns',@optimplotfval);
        %         [FitStruct.FitParameters(1,:),FitStruct.EndFunctionValue(1)] = fminsearchbnd(MinimizationFunction,...
        %             InitialGuess,LowerBounds,UpperBounds,fminoptions);
        [FitStruct.FitParameters(1,:),FitStruct.EndFunctionValue(1)] = fminsearchbnd(MinimizationFunction,...
            InitialGuess,LowerBounds,UpperBounds);
        FitStruct.LowerBounds = LowerBounds;
        FitStruct.UpperBounds = UpperBounds;
        FitStruct.InitialGuess = InitialGuess;
        FitStruct.ChiFittingWeight = ProcessParameters.ChiFittingWeight;
        FitStruct.RamanEnergyFittingWeight = ProcessParameters.RamanEnergyFittingWeight;

        CEF_Hamiltonian_Initial = CEF_Hamiltonian(InitialGuess);
        [~,EigenValues_Initial] = eig(CEF_Hamiltonian_Initial);
        EigenValues_Initial = sort(EigenValues_Initial);
        FitStruct.EigenValues_Initial = EigenValues_Initial - EigenValues_Initial(1);


        disp('Iteration 1')
        disp('=================================================');
        disp('Fit Parameters:');
        disp(FitStruct.FitParameters(1,:));
        disp('End Function Value:');
        disp(FitStruct.EndFunctionValue(1));
        disp('=================================================');

        FitIterations = ProcessParameters.FitIterations;
        if FitIterations>0
            for z = 2:(FitIterations+1)
                [FitStruct.FitParameters(z,:),FitStruct.EndFunctionValue(z)] = fminsearchbnd(MinimizationFunction,...
                    FitStruct.FitParameters(z-1,:),LowerBounds,UpperBounds);
                disp(['Iteration ',num2str(z)]);
                disp('=================================================');
                disp('Fit Parameters:');
                disp(FitStruct.FitParameters(z,:));
                disp('End Function Value:');
                disp(FitStruct.EndFunctionValue(z));
                disp('=================================================');
            end
        end

        % Display the results.


        % Calculate the relevant magnetotropic traces at different
        % temperatures and plot them on top of the data. Need to first
        % calculate the Hamiltonian as it is used to calculate the k.
        FitParameters = FitStruct.FitParameters(FitIterations+1,:);
        CEF_Hamiltonian_MatrixFit = CEF_Hamiltonian(FitParameters);
        [~,EigenValues_Fit] = eig(CEF_Hamiltonian_MatrixFit);
        [EnergyValues_Fit,~] = sort(diag(EigenValues_Fit));
        ZeroEnergy_Fit = EnergyValues_Fit(1);

        disp('Eigenvalue differences: ');
        disp((EnergyValues_Fit(1:2:end)' - ZeroEnergy_Fit) - ProcessParameters.Eigenenergies);
        FitStruct.EigenValues_Final = EnergyValues_Fit-ZeroEnergy_Fit;
        FitStruct.Eigenvalue_Differences = (EnergyValues_Fit(1:2:end)' - ZeroEnergy_Fit) - ProcessParameters.Eigenenergies;

        % We will plot the resulting fit at fields between 0 and 60T.
        InterpField = 0:0.1:20;
        LengthFitRange = length(ProcessParameters.FitRange) ;
        FitRange = ProcessParameters.FitRange;
        for q=1:LengthFitRange
            k = FitRange(q);
            CurrentTemperature = MagnetotropicData_AB(k).Temperature;

            % Plot the fits.
            kab_Calc = Magnetotropic_AB_vs_H(CEF_Hamiltonian_MatrixFit,...
                ZeroEnergy_Fit,ID,Jx,Jz,Jxx_Parameter,Jzz_Parameter,...
                CurrentTemperature,InterpField,muB,kB,gJ,A1);
            plot(axlist(2),InterpField,kab_Calc,'-','Color',...
                brc([0 0 0],.5),'DisplayName',[num2str(CurrentTemperature),' K Fit']);

            kc_Calc = Magnetotropic_C_vs_H(CEF_Hamiltonian_MatrixFit,...
                ZeroEnergy_Fit,ID,Jx,Jz,Jxx_Parameter,Jzz_Parameter,...
                CurrentTemperature,InterpField,muB,kB,gJ,A1);
            plot(axlist(3),InterpField,kc_Calc,'-','Color',...
                brc([0 0 0],.5),'DisplayName',[num2str(CurrentTemperature),' K Fit']);
        end

        % Plot the susceptibility as well.
        axes(axlist(1)); hold on;
        TemperaturesPlot = 1:300;
        Inverse_ChiAB_Calc = InverseSusceptibility_AB_vs_T(CEF_Hamiltonian_MatrixFit,...
            ZeroEnergy_Fit,ID,Jx,Jxx_Parameter,TemperaturesPlot,muB,kB,gJ,C0);
        Inverse_ChiC_Calc = InverseSusceptibility_C_vs_T(CEF_Hamiltonian_MatrixFit,...
            ZeroEnergy_Fit,ID,Jz,Jzz_Parameter,TemperaturesPlot,muB,kB,gJ,C0);
        plot(TemperaturesPlot,Inverse_ChiAB_Calc,'b');
        plot(TemperaturesPlot,Inverse_ChiC_Calc,'r');
        outcell = FitStruct;

    case 'Fit Chi, Magnetotropic, with Bound Eigenenergies, with Bound Chi Params'

        % First, plot the susceptibility and magnetotropic data to be fit.
        h2 = figure;
        [axlist,~] = tight_subplot(1,3,[.01,.04],[.10 .10],[0.075 .075]);

        % Plot susceptibility in subplot 1.
        figure(h2);
        axes(axlist(1)); hold on;
        plot(SusceptibilityData.Temperature,SusceptibilityData.Inverse_ChiAB,...
            'bsq','MarkerFaceColor',brc([0 0 1],.5),'DisplayName','H||ab');
        plot(SusceptibilityData.Temperature,SusceptibilityData.Inverse_ChiC,...
            'rsq','MarkerFaceColor',brc([1 0 0],.5),'DisplayName','H||c');
        xlabel('T [K]'); ylabel('\chi^{-1} [mol Oe/emu]');
        axes(axlist(2)); hold on;
        axes(axlist(3)); hold on;
        for q=1:length(MagnetotropicData_AB)
            % Plot AB data in subplot 2.
            plot(axlist(2),MagnetotropicData_AB(q).InterpFields,...
                MagnetotropicData_AB(q).k_interp,'sq',...
                'Color', MagnetotropicData_AB(q).UColor,...
                'MarkerFaceColor',brc(MagnetotropicData_AB(q).UColor,.5),...
                'DisplayName',num2str(MagnetotropicData_AB(q).Temperature));
            % Plot C data in subplot 3.
            plot(axlist(3),MagnetotropicData_C(q).InterpFields,...
                MagnetotropicData_C(q).k_interp,'sq',...
                'Color', MagnetotropicData_C(q).UColor,...
                'MarkerFaceColor',brc(MagnetotropicData_C(q).UColor,.5),...
                'DisplayName',num2str(MagnetotropicData_C(q).Temperature));
        end
        axes(axlist(2)); xlabel('\mu_0H [T]'); ylabel('k_{ab}');
        axes(axlist(3)); xlabel('\mu_0H [T]'); ylabel('k_c');


        % Define the cumulative error function. First, need to define all
        % relevant inputs to the Cumulative Error function.
        MinimizationFunction = @(Bmn) CumulativeError_BoundEigenvalues_AndChiParams(Bmn,...
            Jxx_Parameter,Jzz_Parameter,chi0_Parameter,...
            CEF_Hamiltonian,ID,Jx,Jz,{MagnetotropicData_AB,MagnetotropicData_C},...
            SusceptibilityData,muB,kB,gJ,A1,C0,ProcessParameters);

        % fminsearchbnd needs the minimization function, an initial guess
        % at the parameters, lower bounds, and upper bounds. Define
        % FitStruct, which holds information about the fitting routine.
        InitialGuess = Bmn_Parameters;
        LowerBounds = -100.*ones(1,length(Bmn_Parameters));
        UpperBounds = 100.*ones(1,length(Bmn_Parameters));
%         LowerBounds = InitialGuess;
%         UpperBounds = InitialGuess;
%         
        %         LowerBounds = [-.17,-100.*ones(1,length(Bmn_Parameters)-1)];
        %         UpperBounds = [-.17,100.*ones(1,length(Bmn_Parameters)-1)];

        % %         fminoptions = optimset('PlotFcns',@optimplotfval);
        %         [FitStruct.FitParameters(1,:),FitStruct.EndFunctionValue(1)] = fminsearchbnd(MinimizationFunction,...
        %             InitialGuess,LowerBounds,UpperBounds,fminoptions);
        [FitStruct.FitParameters(1,:),FitStruct.EndFunctionValue(1)] = fminsearchbnd(MinimizationFunction,...
            InitialGuess,LowerBounds,UpperBounds);
        FitStruct.LowerBounds = LowerBounds;
        FitStruct.UpperBounds = UpperBounds;
        FitStruct.InitialGuess = InitialGuess;
        FitStruct.ChiFittingWeight = ProcessParameters.ChiFittingWeight;
        FitStruct.RamanEnergyFittingWeight = ProcessParameters.RamanEnergyFittingWeight;

        CEF_Hamiltonian_Initial = CEF_Hamiltonian(InitialGuess);
        [~,EigenValues_Initial] = eig(CEF_Hamiltonian_Initial);
        EigenValues_Initial = sort(EigenValues_Initial);
        FitStruct.EigenValues_Initial = EigenValues_Initial - EigenValues_Initial(1);


        disp('Iteration 1')
        disp('=================================================');
        disp('Fit Parameters:');
        disp(FitStruct.FitParameters(1,:));
        disp('End Function Value:');
        disp(FitStruct.EndFunctionValue(1));
        disp('=================================================');

        FitIterations = ProcessParameters.FitIterations;
        if FitIterations>0
            for z = 2:(FitIterations+1)
                [FitStruct.FitParameters(z,:),FitStruct.EndFunctionValue(z)] = fminsearchbnd(MinimizationFunction,...
                    FitStruct.FitParameters(z-1,:),LowerBounds,UpperBounds);
                disp(['Iteration ',num2str(z)]);
                disp('=================================================');
                disp('Fit Parameters:');
                disp(FitStruct.FitParameters(z,:));
                disp('End Function Value:');
                disp(FitStruct.EndFunctionValue(z));
                disp('=================================================');
            end
        end

        % Display the results.


        % Calculate the relevant magnetotropic traces at different
        % temperatures and plot them on top of the data. Need to first
        % calculate the Hamiltonian as it is used to calculate the k.
        FitParameters = FitStruct.FitParameters(FitIterations+1,:);
        CEF_Hamiltonian_MatrixFit = CEF_Hamiltonian(FitParameters);
        [~,EigenValues_Fit] = eig(CEF_Hamiltonian_MatrixFit);
        [EnergyValues_Fit,~] = sort(diag(EigenValues_Fit));
        ZeroEnergy_Fit = EnergyValues_Fit(1);

        disp('Eigenvalue differences: ');
        disp((EnergyValues_Fit(1:2:end)' - ZeroEnergy_Fit) - ProcessParameters.Eigenenergies);
        FitStruct.EigenValues_Final = EnergyValues_Fit-ZeroEnergy_Fit;
        FitStruct.Eigenvalue_Differences = (EnergyValues_Fit(1:2:end)' - ZeroEnergy_Fit) - ProcessParameters.Eigenenergies;

        % We will plot the resulting fit at fields between 0 and 60T.
        InterpField = 0:0.1:20;
        LengthFitRange = length(ProcessParameters.FitRange) ;
        FitRange = ProcessParameters.FitRange;
        for q=1:LengthFitRange
            k = FitRange(q);
            CurrentTemperature = MagnetotropicData_AB(k).Temperature;

            % Plot the fits.
            kab_Calc = Magnetotropic_AB_vs_H(CEF_Hamiltonian_MatrixFit,...
                ZeroEnergy_Fit,ID,Jx,Jz,Jxx_Parameter,Jzz_Parameter,...
                CurrentTemperature,InterpField,muB,kB,gJ,A1);
            plot(axlist(2),InterpField,kab_Calc,'-','Color',...
                brc([0 0 0],.5),'DisplayName',[num2str(CurrentTemperature),' K Fit']);

            kc_Calc = Magnetotropic_C_vs_H(CEF_Hamiltonian_MatrixFit,...
                ZeroEnergy_Fit,ID,Jx,Jz,Jxx_Parameter,Jzz_Parameter,...
                CurrentTemperature,InterpField,muB,kB,gJ,A1);
            plot(axlist(3),InterpField,kc_Calc,'-','Color',...
                brc([0 0 0],.5),'DisplayName',[num2str(CurrentTemperature),' K Fit']);
        end

        % Plot the susceptibility as well.
        axes(axlist(1)); hold on;
        TemperaturesPlot = 1:300;
        Inverse_ChiAB_Calc = InverseSusceptibility_AB_vs_T(CEF_Hamiltonian_MatrixFit,...
            ZeroEnergy_Fit,ID,Jx,Jxx_Parameter,TemperaturesPlot,muB,kB,gJ,C0);
        Inverse_ChiC_Calc = InverseSusceptibility_C_vs_T(CEF_Hamiltonian_MatrixFit,...
            ZeroEnergy_Fit,ID,Jz,Jzz_Parameter,TemperaturesPlot,muB,kB,gJ,C0);
        plot(TemperaturesPlot,Inverse_ChiAB_Calc,'b');
        plot(TemperaturesPlot,Inverse_ChiC_Calc,'r');
        outcell = FitStruct;

    case 'Display Fit Results k'

        % First, plot the susceptibility and magnetotropic data to be fit.
        h2 = figure;
        [axlist,~] = tight_subplot(1,3,[.01,.04],[.10 .10],[0.075 .075]);
        for i=1:length(axlist)
            set(axlist(i),'YTickLabelMode','auto');
            set(axlist(i),'XTickLabelMode','auto');
        end
        % Plot susceptibility in subplot 1.
        figure(h2);
        axes(axlist(1)); hold on;
        plot(SusceptibilityData.Temperature,SusceptibilityData.Inverse_ChiAB,...
            'bsq','MarkerFaceColor',brc([0 0 1],.5),'DisplayName','H||ab');
        plot(SusceptibilityData.Temperature,SusceptibilityData.Inverse_ChiC,...
            'rsq','MarkerFaceColor',brc([1 0 0],.5),'DisplayName','H||c');
        xlabel('T [K]'); ylabel('\chi^{-1} [mol Oe/emu]');
        axes(axlist(2)); hold on;
        axes(axlist(3)); hold on;
        for q=1:length(MagnetotropicData_AB)
            % Plot AB data in subplot 2.
            plot(axlist(2),MagnetotropicData_AB(q).InterpFields,...
                MagnetotropicData_AB(q).k_interp,'sq',...
                'Color', MagnetotropicData_AB(q).UColor,...
                'MarkerFaceColor',brc(MagnetotropicData_AB(q).UColor,.5),...
                'DisplayName',num2str(MagnetotropicData_AB(q).Temperature));
            % Plot C data in subplot 3.
            plot(axlist(3),MagnetotropicData_C(q).InterpFields,...
                MagnetotropicData_C(q).k_interp,'sq',...
                'Color', MagnetotropicData_C(q).UColor,...
                'MarkerFaceColor',brc(MagnetotropicData_C(q).UColor,.5),...
                'DisplayName',num2str(MagnetotropicData_C(q).Temperature));
        end
        axes(axlist(2)); xlabel('\mu_0H [T]'); ylabel('k_{ab}');
        axes(axlist(3)); xlabel('\mu_0H [T]'); ylabel('k_c');


        % Define the cumulative error function. First, need to define all
        % relevant inputs to the Cumulative Error function.
        %         MinimizationFunction = @(Bmn) CumulativeError(Bmn,...
        %             Jxx_Parameter,Jzz_Parameter,chi0_Parameter,...
        %             CEF_Hamiltonian,ID,Jx,Jz,{MagnetotropicData_AB,MagnetotropicData_C},...
        %             SusceptibilityData,muB,kB,gJ,A1,C0,ProcessParameters);
        %
        % fminsearchbnd needs the minimization function, an initial guess
        % at the parameters, lower bounds, and upper bounds.
        %         InitialGuess = Bmn_Parameters;
        %         LowerBounds = repmat(-1,1,length(Bmn_Parameters));
        %         UpperBounds = repmat(1,1,length(Bmn_Parameters));
        %         [FitParameters,EndFunctionValue] = fminsearchbnd(MinimizationFunction,...
        %             InitialGuess,LowerBounds,UpperBounds);
        %
        % Display the results.
        %         display(FitParameters);
        %         display(EndFunctionValue);

        FitParameters = InputParameters(1:6);
        % Calculate the relevant magnetotropic traces at different
        % temperatures and plot them on top of the data. Need to first
        % calculate the Hamiltonian as it is used to calculate the k.
        CEF_Hamiltonian_MatrixFit = CEF_Hamiltonian(FitParameters);
        [~,EigenValues_Fit] = eig(CEF_Hamiltonian_MatrixFit);
        [EnergyValues_Fit,~] = sort(diag(EigenValues_Fit));
        ZeroEnergy_Fit = EnergyValues_Fit(1);
        disp(EnergyValues_Fit-ZeroEnergy_Fit);

        % We will plot the resulting fit at fields between 0 and 60T.
        InterpField = 0:0.1:100;
        LengthFitRange = length(ProcessParameters.FitRange) ;
        FitRange = ProcessParameters.FitRange;
        for q=1:LengthFitRange
            k = FitRange(q);
            CurrentTemperature = MagnetotropicData_AB(k).Temperature;

            % Plot the fits.
            kab_Calc = Magnetotropic_AB_vs_H(CEF_Hamiltonian_MatrixFit,...
                ZeroEnergy_Fit,ID,Jx,Jz,Jxx_Parameter,Jzz_Parameter,...
                CurrentTemperature,InterpField,muB,kB,gJ,A1);
            plot(axlist(2),InterpField,kab_Calc,'-','Color',...
                brc([0 0 0],.5),'DisplayName',[num2str(CurrentTemperature),' K Fit']);
            

            kc_Calc = Magnetotropic_C_vs_H(CEF_Hamiltonian_MatrixFit,...
                ZeroEnergy_Fit,ID,Jx,Jz,Jxx_Parameter,Jzz_Parameter,...
                CurrentTemperature,InterpField,muB,kB,gJ,A1);
            plot(axlist(3),InterpField,kc_Calc,'-','Color',...
                brc([0 0 0],.5),'DisplayName',[num2str(CurrentTemperature),' K Fit']);
            
        end
        axes(axlist(1)); hold on;
        TemperaturesPlot = 1:300;
        Inverse_ChiAB_Calc = InverseSusceptibility_AB_vs_T(CEF_Hamiltonian_MatrixFit,...
            ZeroEnergy_Fit,ID,Jx,Jxx_Parameter,TemperaturesPlot,muB,kB,gJ,C0);
        Inverse_ChiC_Calc = InverseSusceptibility_C_vs_T(CEF_Hamiltonian_MatrixFit,...
            ZeroEnergy_Fit,ID,Jz,Jzz_Parameter,TemperaturesPlot,muB,kB,gJ,C0);
        plot(TemperaturesPlot,Inverse_ChiAB_Calc,'b');
        plot(TemperaturesPlot,Inverse_ChiC_Calc,'r');

        [AX0,AZ0,A0_AB,A0_C] = Calculate_AZ(Bmn_Parameters,...
            Jxx_Parameter,Jzz_Parameter,chi0_Parameter,...
            CEF_Hamiltonian,ID,Jx,Jz,muB,gJ)
%         xlim([-5,5]);
        ylim([0,600]);
        axes(axlist(2));
        ylim([-5,5]);
        xlim([0,100]);
        axes(axlist(3));
        ylim([-5,5]);
        xlim([0,100]);

    case 'Find Scaling Amplitude'
        % Use this case to calibrate magnetotropic data to measured
        % susceptibility.
        h2 = figure;
        [axlist,~] = tight_subplot(1,3,[.01,.04],[.10 .10],[0.075 .075]);

        % Plot susceptibility in subplot 1.
        figure(h2);
        axes(axlist(1)); hold on;
        plot(SusceptibilityData.Temperature,...
            SusceptibilityData.ChiAB_meV_Per_T2 - SusceptibilityData.ChiC_meV_Per_T2,...
            '-ok','MarkerFaceColor',brc([0 0 0],.5),'DisplayName','Delta \chi [meV/T^2]');
        xlabel('T [K]'); ylabel('\Delta \chi [meV/T^2]');
        set(gca,'YTickLabelMode','auto');
        set(gca,'XTickLabelMode','auto');
        axes(axlist(2)); hold on;
        axes(axlist(3)); hold on;
        CutoffField = 15;
        for q=1:length(MagnetotropicData_AB)
            % Plot AB data in subplot 2. Prepare Fitinds for fitting for
            % delta chi.
            x = MagnetotropicData_AB(q).C_Field;
            y = MagnetotropicData_AB(q).C_frequency;
            FitInds = x<=CutoffField;
            plot(axlist(2),x,y,'sq',...
                'Color', MagnetotropicData_AB(q).UColor,...
                'MarkerFaceColor',brc(MagnetotropicData_AB(q).UColor,.5),...
                'DisplayName',num2str(MagnetotropicData_AB(q).Temperature));

            % Fit the low field data to get the coefficients.
            [xData, yData] = prepareCurveData( x(FitInds), y(FitInds) );
            ft = fittype( 'a0+b0.*x.^2', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.DiffMaxChange = 1000;
            opts.DiffMinChange = 1e-10;
            opts.Display = 'Off';
            opts.MaxFunEvals = 1000;
            opts.MaxIter = 1000;
            opts.TolFun = 1e-10;
            opts.TolX = 1e-10;
            [FitResult, ~] = fit( xData, yData, ft, opts );
            CoefficientValues = coeffvalues(FitResult);
            FitData.Offset_AB(q) = CoefficientValues(1);
            FitData.Parabola_AB(q) = CoefficientValues(2);
            FitData.Temperature(q) = MagnetotropicData_AB(q).Temperature;
            hold on; plot(axlist(2),x,FitResult(x)','-k','Linewidth',1);
            FitData.Parabola_AB_Units(q) = FitData.Parabola_AB(q).*....
                MagnetotropicData_AB(q).zf_frequency./(1e6);
            % Plot C data in subplot 3.
            x = MagnetotropicData_C(q).C_Field;
            y = MagnetotropicData_C(q).C_frequency;
            FitInds = x<=CutoffField;
            plot(axlist(3),x,y,'sq',...
                'Color', MagnetotropicData_C(q).UColor,...
                'MarkerFaceColor',brc(MagnetotropicData_C(q).UColor,.5),...
                'DisplayName',num2str(MagnetotropicData_C(q).Temperature));

            % Fit the low field data to get the coefficients.
            [xData, yData] = prepareCurveData( x(FitInds), y(FitInds) );
            ft = fittype( 'a0+b0.*x.^2', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.DiffMaxChange = 1000;
            opts.DiffMinChange = 1e-10;
            opts.Display = 'Off';
            opts.MaxFunEvals = 1000;
            opts.MaxIter = 1000;
            opts.TolFun = 1e-10;
            opts.TolX = 1e-10;
            [FitResult, ~] = fit( xData, yData, ft, opts );
            CoefficientValues = coeffvalues(FitResult);
            FitData.Offset_C(q) = CoefficientValues(1);
            FitData.Parabola_C(q) = CoefficientValues(2);
            FitData.Parabola_C_Units(q) = FitData.Parabola_C(q).*....
                MagnetotropicData_C(q).zf_frequency./(1e6);
            plot(axlist(3),x,FitResult(x)','-k','Linewidth',1);
        end
        axes(axlist(2)); xlabel('\mu_0H [T]'); ylabel('\Deltaf_{ab}');
        set(gca,'YTickLabelMode','auto');
        set(gca,'XTickLabelMode','auto');
        axes(axlist(3)); xlabel('\mu_0H [T]'); ylabel('\Deltaf_c');
        set(gca,'YTickLabelMode','auto');
        set(gca,'XTickLabelMode','auto');
    case 'Find B20 From Susceptibility'
        % Plot susceptibility data in molar susceptibility units.
        figure; hold on;
        plot(SusceptibilityData.FullTemperature_AB,...
            1./SusceptibilityData.ChiAB_cm3_Per_mol,'bsq','MarkerFaceColor',...
            brc([0 0 1],.5),'DisplayName','H||ab');
        plot(SusceptibilityData.FullTemperature_C,...
            1./SusceptibilityData.ChiC_cm3_Per_mol,'rsq','MarkerFaceColor',...
            brc([1 0 0],.5),'DisplayName','H||c');
        %         plot(SusceptibilityData.FullTemperature_AB,...
        %             1./SusceptibilityData.Inverse_ChiAB_Full_molOe_Per_emu,'bsq','MarkerFaceColor',...
        %             brc([0 0 1],.5),'DisplayName','H||ab');
        %         plot(SusceptibilityData.FullTemperature_C,...
        %             1./SusceptibilityData.Inverse_ChiC_Full_molOe_Per_emu,'rsq','MarkerFaceColor',...
        %             brc([1 0 0],.5),'DisplayName','H||c');
        xlabel('T [K]'); ylabel('\chi^{-1} [cm^3/mol]^{-1}');
        %         ylabel('\chi^{-1} [mol\cdotOe/emu]');
        xlim([0 310]);
        %         set(gca,'Yscale','Log');set(gca,'Xscale','Log');

        CutTemp = 250;
        CutInds_AB = SusceptibilityData.FullTemperature_AB>=CutTemp;
        XDat_AB = SusceptibilityData.FullTemperature_AB(CutInds_AB);
        YDat_AB = 1./SusceptibilityData.ChiAB_cm3_Per_mol(CutInds_AB);
        CutInds_C = SusceptibilityData.FullTemperature_AB>=CutTemp;
        XDat_C = SusceptibilityData.FullTemperature_AB(CutInds_C);
        YDat_C = 1./SusceptibilityData.ChiC_cm3_Per_mol(CutInds_C);
        FitResult_AB = polyfit(XDat_AB,YDat_AB,1);
        FitResult_C = polyfit(XDat_C,YDat_C,1);

        PlotTemperatures = -100:.1:300;
        hold on; plot(PlotTemperatures,polyval(FitResult_AB,PlotTemperatures),...
            '-k','Linewidth',1.5,'DisplayName','ABFit');
        hold on; plot(PlotTemperatures,polyval(FitResult_C,PlotTemperatures),...
            '-k','Linewidth',1.5,'DisplayName','CFit');
        xlim([-60, 300])

        B20Func_Par = @(theta,Jfit,j) -1.25.*8.62e-2.*(theta + (j*(j+1))*2*Jfit)*(1./((j-.5)*(j+1.5)));
        B20Func_Perp = @(theta,Jfit,j) (5/2).*8.62e-2.*(theta + (j*(j+1))*2*Jfit)*(1./((j-.5)*(j+1.5)));


        text(0,10,...
            ['\Theta_{AB} = ',num2str(roots(FitResult_AB)),' K; B_2^0 = ',...
            num2str(B20Func_Perp(roots(FitResult_AB),Jxx_Parameter,7/2)),' meV'],...
            'FontSize',20);
        text(0,8,...
            ['\Theta_C = ',num2str(roots(FitResult_C)),' K; B_2^0 = ',...
            num2str(B20Func_Par(roots(FitResult_C),Jzz_Parameter,7/2)),' meV'],...
            'FontSize',20);
    case 'Find B20 From Susceptibility Higher Order'
        % Plot susceptibility data in molar susceptibility units.
        figure; hold on;
        plot(SusceptibilityData.FullTemperature_AB,...
            1./SusceptibilityData.ChiAB_cm3_Per_mol,'bsq','MarkerFaceColor',...
            brc([0 0 1],.5),'DisplayName','H||ab');
        plot(SusceptibilityData.FullTemperature_C,...
            1./SusceptibilityData.ChiC_cm3_Per_mol,'rsq','MarkerFaceColor',...
            brc([1 0 0],.5),'DisplayName','H||c');
        %         plot(SusceptibilityData.FullTemperature_AB,...
        %             1./SusceptibilityData.Inverse_ChiAB_Full_molOe_Per_emu,'bsq','MarkerFaceColor',...
        %             brc([0 0 1],.5),'DisplayName','H||ab');
        %         plot(SusceptibilityData.FullTemperature_C,...
        %             1./SusceptibilityData.Inverse_ChiC_Full_molOe_Per_emu,'rsq','MarkerFaceColor',...
        %             brc([1 0 0],.5),'DisplayName','H||c');
        xlabel('T [K]'); ylabel('\chi^{-1} [cm^3/mol]^{-1}');
        %         ylabel('\chi^{-1} [mol\cdotOe/emu]');
        xlim([0 310]);
        %         set(gca,'Yscale','Log');set(gca,'Xscale','Log');

        CutTemp = 250;
        CutInds_AB = SusceptibilityData.FullTemperature_AB>=CutTemp;
        XDat_AB = SusceptibilityData.FullTemperature_AB(CutInds_AB);
        YDat_AB = 1./SusceptibilityData.ChiAB_cm3_Per_mol(CutInds_AB);
        CutInds_C = SusceptibilityData.FullTemperature_AB>=CutTemp;
        XDat_C = SusceptibilityData.FullTemperature_AB(CutInds_C);
        YDat_C = 1./SusceptibilityData.ChiC_cm3_Per_mol(CutInds_C);

        [XDat_AB, YDat_AB] = prepareCurveData( XDat_AB, YDat_AB);
        ft = fittype( 'An1*8.62e-2.*x.^2./(x+TCW)', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        [FitResult_AB, ~] = fit( XDat_AB, YDat_AB, ft, opts );


        [XDat_C, YDat_C] = prepareCurveData( XDat_C, YDat_C);
        ft = fittype( 'An1*8.62e-2.*x.^2./(x+TCW)', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        [FitResult_C, ~] = fit( XDat_C, YDat_C, ft, opts );

        PlotTemperatures = 50:.1:300;
        hold on; plot(PlotTemperatures,FitResult_AB(PlotTemperatures)',...
            '-k','Linewidth',1.5,'DisplayName','ABFit');
        hold on; plot(PlotTemperatures,FitResult_C(PlotTemperatures)',...
            '-k','Linewidth',1.5,'DisplayName','CFit');
        xlim([-60, 300]);

        AB_Coefficients = coeffvalues(FitResult_AB);
        C_Coefficients  = coeffvalues(FitResult_C);

        B20Func_Par = @(theta,Jfit,j) -1.25.*8.62e-2.*(theta + (j*(j+1))*2*Jfit)*(1./((j-.5)*(j+1.5)));
        B20Func_Perp = @(theta,Jfit,j) (5/2).*8.62e-2.*(theta + (j*(j+1))*2*Jfit)*(1./((j-.5)*(j+1.5)));


        text(0,10,...
            ['\Theta_{AB} = ',num2str(AB_Coefficients(2)),' K; B_2^0 = ',...
            num2str(B20Func_Perp(AB_Coefficients(2),Jxx_Parameter,7/2)),' meV'],...
            'FontSize',20);
        text(0,8,...
            ['\Theta_C = ',num2str(C_Coefficients(2)),' K; B_2^0 = ',...
            num2str(B20Func_Par(C_Coefficients(2),Jzz_Parameter,7/2)),' meV'],...
            'FontSize',20);
    case 'Find B20 From Susceptibility Higher Order Variable Cutoff'
        % Plot susceptibility data in molar susceptibility units.
        figure;
        B20Func_Par = @(theta,Jfit,j) -1.25.*8.62e-2.*(theta + (j*(j+1))*2*Jfit)*(1./((j-.5)*(j+1.5)));
        B20Func_Perp = @(theta,Jfit,j) (5/2).*8.62e-2.*(theta + (j*(j+1))*2*Jfit)*(1./((j-.5)*(j+1.5)));

        CutoffTemperatures = 150:2:260;
        NumTemps = length(CutoffTemperatures);
        for q = 1:length(CutoffTemperatures)
            CutTemp = CutoffTemperatures(q);
            CutInds_AB = SusceptibilityData.FullTemperature_AB>=CutTemp;
            XDat_AB = SusceptibilityData.FullTemperature_AB(CutInds_AB);
            YDat_AB = 1./SusceptibilityData.ChiAB_cm3_Per_mol(CutInds_AB);
            CutInds_C = SusceptibilityData.FullTemperature_AB>=CutTemp;
            XDat_C = SusceptibilityData.FullTemperature_AB(CutInds_C);
            YDat_C = 1./SusceptibilityData.ChiC_cm3_Per_mol(CutInds_C);

            [XDat_AB, YDat_AB] = prepareCurveData( XDat_AB, YDat_AB);
            %             ft = fittype( '.345*8.62e-2.*x.^2./(x+TCW)', 'independent', 'x', 'dependent', 'y' );
            ft = fittype( 'An1*8.62e-2.*x.^2./(x+TCW)', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            [FitResult_AB, ~] = fit( XDat_AB, YDat_AB, ft, opts );


            [XDat_C, YDat_C] = prepareCurveData( XDat_C, YDat_C);
            ft = fittype( 'An1*8.62e-2.*x.^2./(x+TCW)', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            [FitResult_C, ~] = fit( XDat_C, YDat_C, ft, opts );


            AB_Coefficients = coeffvalues(FitResult_AB);
            C_Coefficients  = coeffvalues(FitResult_C);

            FitStruct.Theta_AB(q) = AB_Coefficients(2);
            FitStruct.Theta_C(q)  = C_Coefficients(2);
            FitStruct.B20_ABVals(q) = B20Func_Perp(FitStruct.Theta_AB(q),Jxx_Parameter,7/2);
            FitStruct.B20_CVals(q) = B20Func_Par(FitStruct.Theta_C(q),Jzz_Parameter,7/2);
        end

        subplot(1,2,1); hold on; plot(CutoffTemperatures,FitStruct.Theta_AB,'-b',...
            'DisplayName','AB');
        hold on; plot(CutoffTemperatures,FitStruct.Theta_C,'-r',...
            'DisplayName','C');
        xlabel('Fit Above T [K]');
        ylabel('\Theta_\alpha [K]');

        subplot(1,2,2); hold on; plot(CutoffTemperatures,FitStruct.B20_ABVals,'-b',...
            'DisplayName','AB');
        hold on; plot(CutoffTemperatures,FitStruct.B20_CVals,'-r',...
            'DisplayName','C');
        xlabel('Fit Above T [K]');
        ylabel('B_2^0 [meV]');


    case 'Find B20 From Susceptibility Simultaneous'
        % Plot susceptibility data in molar susceptibility units.
        B20Func_Par = @(theta,Jfit,j) -1.25.*8.62e-2.*(theta + (j*(j+1))*2*Jfit)*(1./((j-.5)*(j+1.5)));
        B20Func_Perp = @(theta,Jfit,j) (5/2).*8.62e-2.*(theta + (j*(j+1))*2*Jfit)*(1./((j-.5)*(j+1.5)));
        kb = 8.62e-2;
        ThetaC = @(B20,JFit,j) -(4./(5*kb)).*(j-.5).*(j+1.5).*B20 - j.*(j+1).*2.*JFit;
        ThetaAB = @(B20,JFit,j) +(2./(5*kb)).*(j-.5).*(j+1.5).*B20 - j.*(j+1).*2.*JFit;
        FitFunc_AB = @(params,XValues) params(1).*(kb).*XValues.^2./(XValues+ThetaAB(params(2),Jxx_Parameter,7/2));
        FitFunc_C = @(params,XValues) params(1).*(kb).*XValues.^2./(XValues+ThetaC(params(2),Jzz_Parameter,7/2));

        XData = SusceptibilityData.InterpTemp;
        YData_AB = SusceptibilityData.Interp_Inverse_ChiAB_cm3_Per_mol;
        YData_C  = SusceptibilityData.Interp_Inverse_ChiC_cm3_Per_mol;
        figure; hold on; plot(XData,YData_AB,'-b');
        hold on; plot(XData,YData_C,'-r');
        interptemp = 60:400;
        CutoffTemperatures = 150:2:292;
        for q=1:length(CutoffTemperatures)

            CutTemp = CutoffTemperatures(q);
            CutInds = XData>=CutTemp;
            XData = XData(CutInds);
            YData_AB = YData_AB(CutInds);
            YData_C  = YData_C(CutInds);

            DenomAB = sum((YData_AB-mean(YData_AB)).^2);
            DenomC = sum((YData_C-mean(YData_C)).^2);

            ErrorFunction = @(params) ( sum((YData_AB -...
                FitFunc_AB(params,XData)).^2)./DenomAB ) +...
                ( sum((YData_C -...
                FitFunc_C(params,XData)).^2)./DenomC );

            LowerBounds = [0,-1];
            UpperBounds = [1e6,1];
            InitialGuess= [1,+.01];
            [FitStruct.FitParameters(q,:),FitStruct.FValEnd(q)] = fminsearchbnd(ErrorFunction,...
                InitialGuess,LowerBounds,UpperBounds);
            FitStruct.B20_Result(q) = FitStruct.FitParameters(q,2);
            FitStruct.ThetaAB_Value(q) = ThetaAB(FitStruct.B20_Result(q),Jxx_Parameter,7/2);
            FitStruct.ThetaC_Value(q) = ThetaC(FitStruct.B20_Result(q),Jzz_Parameter,7/2);
            FitStruct.B20_Perp_Calc(q) = B20Func_Perp(FitStruct.ThetaAB_Value(q) ,Jxx_Parameter,7/2);
            FitStruct.B20_Par_Calc(q) = B20Func_Par(FitStruct.ThetaC_Value(q),Jzz_Parameter,7/2);
            hold on; plot(interptemp,FitFunc_AB(FitStruct.FitParameters(q,:),interptemp),'-k','Linewidth',1.25)
            hold on; plot(interptemp,FitFunc_C(FitStruct.FitParameters(q,:),interptemp),'-k','Linewidth',1.25)
        end
        FitStruct.Temperatures = CutoffTemperatures;


        figure; hold on;

        subplot(1,3,1);
        hold on; plot(FitStruct.Temperatures,FitStruct.B20_Result,'-ok',...
            'MarkerFaceColor',brc([0 0 0],.5));
        xlabel('Fit Above T [K]'); ylabel('B_2^0 [meV]');
        yline(-.17,'-k','Linewidth',1.25)

        subplot(1,3,2);
        hold on; plot(FitStruct.Temperatures,FitStruct.ThetaAB_Value,'-ob',...
            'MarkerFaceColor',brc([0 0 1],.5),'DisplayName','AB');
        hold on; plot(FitStruct.Temperatures,FitStruct.ThetaC_Value,'-or',...
            'MarkerFaceColor',brc([1 0 0],.5),'DisplayName','C');
        xlabel('Fit Above T [K]'); ylabel('\Theta_\alpha [K]');

        subplot(1,3,3);
        hold on; plot(FitStruct.Temperatures,FitStruct.FValEnd,'-ok',...
            'MarkerFaceColor',brc([0 0 0],.5));
        xlabel('Fit Above T [K]'); ylabel('EndFuncVal \propto Error');


    otherwise
        disp('Invalid AnalysisOption in CEF_MagnetotropicFittingAndCalculation!');
end


end



%% Internal Functions Section.


function [ID,Jx,Jy,Jz,Jp,Jm] = GenerateSpinMatrices()
% This internal function populates the matrices to be used in the
% calculations for a spin 7/2 system. The function has no inputs and
% outputs the relevant matrices [ID,Jz,Jp,Jm,Jx,Jy] = [Identity, Jx-Jz
% operators, Raising Operator, and Lowering Operator].

% 8x8 identity matrix. Locally defining the identity is faster than using
% eye to create the identity wherever needed.
ID = eye(8);

% Jz is diagonal in this {|7/2,7/2>.,...,|7/2,-7/2>| basis. Create the
% matrix and populate the diagonals with the Jz values (-7/2...7/2).
Jz = diag([7/2,5/2,3/2,1/2,-1/2,-3/2,-5/2,-7/2]);

% Ladder operator +- J matrices. They have the structure of connecting
% neighboring spin states so should have entries on the super diagonals
% (one diagonal above (+) or below (-) the center). We also make these with
% the diag() matlab function.
Jp = diag([sqrt(7),sqrt(12),sqrt(15),4,sqrt(15),sqrt(12),sqrt(7)],1);
Jm = diag([sqrt(7),sqrt(12),sqrt(15),4,sqrt(15),sqrt(12),sqrt(7)],-1);

% Construct the Jx and Jy operators to be used in the Hamiltonian from sums
% and differences of the raising and lowering operators.
Jx = (Jp+Jm)./2;
Jy = (Jp-Jm)/2i;
end

function [muB,kB,gJ,A1,C0,SecondDerivative_C1,SecondDerivative_C2] = GenerateConstants()
% This internal function generates a set of physical and mathematical
% constants used throughout the calculations. Units of physical constants
% are listed where relevant.

% Constants used throughout the code.
muB = 5.78828e-2; % [meV/T]; Bohr Magneton
kB  = 8.617e-2  ; % [meV/K]; Boltzmann constant.
gJ = 8/7; % L=3, S=1/2, J=7/2 g-lande factor;
A1 = 6*kB/(muB*gJ);
C0 = 2.0416; % Factor for scaling the susceptibility, in units of [mol Oe/(K emu)]

% These vectors are used in numerical calculation of second derivatives.
SecondDerivative_C1 = [1/12, 	-2/3, 	0,      2/3, 	-1/12];
SecondDerivative_C2 = [-1/12 	4/3 	-5/2 	4/3 	-1/12];

end

function [OutcellAB,OutcellC,SusceptibilityData,ProcessParameters]=LoadSampleData(SampleOption)
% This internal function populates the experimental magnetotropic data, the
% experimental susceptibility data, and a set of process parameters that
% are specific to each sample.

% Calling external drivers where possible, load magnetotropic and
% susceptibility data.
BohrMagInEMU=9.27400968e-21;
BohrMagInMevPerT = 5.788e-2;
AvogadrosNumber = 6.022e23;
switch SampleOption
    case 'CsYbSe2'
        [OutcellAB] = CsYbSe2_PaperData_LANL_Driver('AB T Dependence',...
            'Down Rough Delta rmoutliers then smooth no plot','Temperature','No');
        [OutcellC] = CsYbSe2_PaperData_LANL_Driver('C T Dependence',...
            'Down Rough Delta rmoutliers then smooth no plot','Temperature','No');

        % Define ProcessParameters.

        ProcessParameters.Amp = 9.95; % Amplitude for scaling the magnetotropic data.

        ProcessParameters.FitRange = 2:8; % Which datasets are fit? This was hardcoded previously, but it should be variable.
        ProcessParameters.FitIterations = 3; % In new fitting routine, how many additional iterations of fminsearchbnd?

        % Need a weighting factor for the fitting of magnetotropic data.
        % It was previously set to 10, but I wonder if this changes if the
        % residuals are normalized -- check that!
        ProcessParameters.ChiFittingWeight = 0;
        ProcessParameters.RamanEnergyFittingWeight = [100 100 50 10];

        % Eigenenergies determined from raman data. Bound eigenvalues to be
        % similar to these.
        ProcessParameters.Eigenenergies = [0 114.063 184.983 271]./8.06;

        % Chris previously hardcoded these values of the susceptibility in
        % -- I don't have access to the CsYbSe2 susceptibbility curves, so
        % I am keeping this the same for now. :[ I think this is in units
        % of inverse chi ~ [mol Oe/emu]
        SusceptibilityData.Temperature = [4,6,12,20,30,50,70,100:50:250];
        SusceptibilityData.Inverse_ChiAB = [10.1772   11.5966   15.7263   21.0341   27.3460   38.4183   47.9779   60.7196   80.8536  100.1618  119.1967];
        SusceptibilityData.Inverse_ChiC  = [17.4148   20.1512   25.1159   28.4515   30.8158   34.0760   37.5031   43.6232   57.1081   71.9755   87.9933];

        % To convert susceptibility from molar to mub per FU T, we need to
        % multiply by 10000 Oe/T and divide by (N_A * \muB in emu) where
        % N_A is avogadros number and muB [emu] = 9.27400968e-21 muBs per
        % emu.
        SusceptibilityData.ChiAB_meV_Per_T2 = (1./SusceptibilityData.Inverse_ChiAB).*...
            10000.*BohrMagInMevPerT./(AvogadrosNumber.*BohrMagInEMU);
        SusceptibilityData.ChiC_meV_Per_T2 = (1./SusceptibilityData.Inverse_ChiC).*...
            10000.*BohrMagInMevPerT./(AvogadrosNumber.*BohrMagInEMU);

        % extracted the data from a paper figure, leading to this ugly
        % section.
        PaperData = load('E:\IanComputer\Documents\Physics\Minhyea Lee Research\Chris_CEF\CsYbSe2_ChiData.mat');
        PaperData = PaperData.dc;
        SusceptibilityData.FullTemperature_AB =PaperData.Temperature_ab;
        SusceptibilityData.FullTemperature_C = PaperData.Temperature_c;
        SusceptibilityData.ChiAB_cm3_Per_mol = (PaperData.Chiab).*4.*pi;
        SusceptibilityData.ChiC_cm3_Per_mol = (PaperData.Chic).*4.*pi;
        SusceptibilityData.Inverse_ChiAB_Full_molOe_Per_emu= 1./PaperData.Chiab;
        SusceptibilityData.Inverse_ChiC_Full_molOe_Per_emu= 1./PaperData.Chic;

        LowerInterpLimit = max([min(SusceptibilityData.FullTemperature_AB),min(SusceptibilityData.FullTemperature_C)]);
        UpperInterpLimit = min([max(SusceptibilityData.FullTemperature_AB),max(SusceptibilityData.FullTemperature_C)]);
        SusceptibilityData.InterpTemp = linspace(LowerInterpLimit,UpperInterpLimit,2000);
        SusceptibilityData.Interp_Inverse_ChiAB_cm3_Per_mol = Interp1NonUnique(SusceptibilityData.FullTemperature_AB,...
            1./SusceptibilityData.ChiAB_cm3_Per_mol,SusceptibilityData.InterpTemp);
        SusceptibilityData.Interp_Inverse_ChiC_cm3_Per_mol = Interp1NonUnique(SusceptibilityData.FullTemperature_C,...
            1./SusceptibilityData.ChiC_cm3_Per_mol,SusceptibilityData.InterpTemp);
    case 'KYbSe2'
        [OutcellAB] = KYbSe2_LANL_Driver('AB T Dependence',...
            'Down Rough Delta rmoutliers then smooth no plot','Temperature','No');
        [OutcellC] = KYbSe2_LANL_Driver('C T Dependence',...
            'Down Rough Delta rmoutliers then smooth no plot','Temperature','No');

        % Amplitude for scaling the magnetotropic data to susceptibility
        % data? Updated on 2/4/2022. Scaled 4K magnetotropic data to Chi.
        ProcessParameters.Amp = 3.7509;

        % Pull susceptibility from relevant driver.
        SusceptibilityData.Temperature = [4,6,12,20,30,50,70,100:50:250];
        [Mab_T,~] = KYbSe2_Chi_Nov2021('HParAB',...
            'Skip','No','No');
        [Mc_T,~] = KYbSe2_Chi_Nov2021('HParC',...
            'Skip','No','No');


        ProcessParameters.FitRange = 1:5; % Which datasets are fit? This was hardcoded previously, but it should be variable.
        ProcessParameters.FitIterations = 3; % In new fitting routine, how many additional iterations of fminsearchbnd?

        % Need a weighting factor for the fitting of magnetotropic data.
        % It was previously set to 10, but I wonder if this changes if the
        % residuals are normalized -- check that!
        ProcessParameters.ChiFittingWeight = 0;
        ProcessParameters.RamanEnergyFittingWeight = 100.*[10 10 1 .2];
        ProcessParameters.Chi_PerturbationTerms = [-1.94e-3,-1.53e-4]; % [az0,ax0
        ProcessParameters.Chi_A0Vals            = [.9783,2.903]; % [A0z,A0x]
        ProcessParameters.ChiFitResultsWeight = 1000;

        % Eigenenergies determined from raman data. Bound eigenvalues to be
        % similar to these.
        ProcessParameters.Eigenenergies = [0 114.912 188.759 257.148]./8.06;

        % Interpolate at predetermined density in temperature.
        SusceptibilityData.Inverse_ChiAB = Interp1NonUnique(Mab_T{1}.TemperatureK,...
            Mab_T{1}.InvChi_molOe_Per_emu,SusceptibilityData.Temperature);
        SusceptibilityData.Inverse_ChiC = Interp1NonUnique(Mc_T{1}.TemperatureK,...
            Mc_T{1}.InvChi_molOe_Per_emu,SusceptibilityData.Temperature);

        % Define susceptibility in molar units for use in fitting HT B20
        % value.
        SusceptibilityData.FullTemperature_AB = Mab_T{1}.TemperatureK;
        SusceptibilityData.FullTemperature_C = Mc_T{1}.TemperatureK;
        SusceptibilityData.ChiAB_cm3_Per_mol = (1./Mab_T{1}.InvChi_molOe_Per_emu).*4.*pi;
        SusceptibilityData.ChiC_cm3_Per_mol = (1./Mc_T{1}.InvChi_molOe_Per_emu).*4.*pi;
        SusceptibilityData.Inverse_ChiAB_Full_molOe_Per_emu= 1./Mab_T{1}.InvChi_molOe_Per_emu;
        SusceptibilityData.Inverse_ChiC_Full_molOe_Per_emu= 1./Mc_T{1}.InvChi_molOe_Per_emu;

        SusceptibilityData.ChiAB_meV_Per_T2 = (1./SusceptibilityData.Inverse_ChiAB).*...
            10000.*BohrMagInMevPerT./(AvogadrosNumber.*BohrMagInEMU);
        SusceptibilityData.ChiC_meV_Per_T2 = (1./SusceptibilityData.Inverse_ChiC).*...
            10000.*BohrMagInMevPerT./(AvogadrosNumber.*BohrMagInEMU);

        LowerInterpLimit = max([min(SusceptibilityData.FullTemperature_AB),min(SusceptibilityData.FullTemperature_C)]);
        UpperInterpLimit = min([max(SusceptibilityData.FullTemperature_AB),max(SusceptibilityData.FullTemperature_C)]);
        SusceptibilityData.InterpTemp = linspace(LowerInterpLimit,UpperInterpLimit,2000);
        SusceptibilityData.Interp_Inverse_ChiAB_cm3_Per_mol = Interp1NonUnique(SusceptibilityData.FullTemperature_AB,...
            1./SusceptibilityData.ChiAB_cm3_Per_mol,SusceptibilityData.InterpTemp);
        SusceptibilityData.Interp_Inverse_ChiC_cm3_Per_mol = Interp1NonUnique(SusceptibilityData.FullTemperature_C,...
            1./SusceptibilityData.ChiC_cm3_Per_mol,SusceptibilityData.InterpTemp);

    case 'KYbSe2 18T'
        [OutcellAB,~, ~, ~] = RTMdriver('res vs B','fall KYbSeSmall', 'n', '90', 'n');
        [OutcellC,~, ~, ~] = RTMdriver('res vs B','fall KYbSeSmall', 'n', '0', 'n');

        % Amplitude for scaling the magnetotropic data to susceptibility
        % data? Updated on 2/4/2022. Scaled 4K magnetotropic data to Chi.
        ProcessParameters.Amp = 9.256760738;
%         ProcessParameters.Amp = 7.5489;
        % Pull susceptibility from relevant driver.
        SusceptibilityData.Temperature = [4,6,12,20,30,50,70,100:50:250];
        [Mab_T,~] = KYbSe2_Chi_Nov2021('HParAB',...
            'Skip','No','No');
        [Mc_T,~] = KYbSe2_Chi_Nov2021('HParC',...
            'Skip','No','No');


        ProcessParameters.FitRange = 2:10; % Which datasets are fit? This was hardcoded previously, but it should be variable.
        ProcessParameters.FitIterations = 0; % In new fitting routine, how many additional iterations of fminsearchbnd?

        % Need a weighting factor for the fitting of magnetotropic data.
        % It was previously set to 10, but I wonder if this changes if the
        % residuals are normalized -- check that!
        ProcessParameters.ChiFittingWeight = 3;
        ProcessParameters.RamanEnergyFittingWeight = 30.*[10 10 1 .2];

        ProcessParameters.Chi_PerturbationTerms = [-1.94e-3,-1.53e-4]; % [az0,ax0
%         ProcessParameters.Chi_A0Vals            = [.9783,2.903]; % [A0z,A0x]
        ProcessParameters.Chi_A0Vals            = [1.0124,2.808]; % [A0z,A0x]
        
        ProcessParameters.ChiFitResultsWeight = 50;

        % Eigenenergies determined from raman data. Bound eigenvalues to be
        % similar to these.
        ProcessParameters.Eigenenergies = [0 114.912 188.759 257.148]./8.06;

        % Interpolate at predetermined density in temperature.
        SusceptibilityData.Inverse_ChiAB = Interp1NonUnique(Mab_T{1}.TemperatureK,...
            Mab_T{1}.InvChi_molOe_Per_emu,SusceptibilityData.Temperature);
        SusceptibilityData.Inverse_ChiC = Interp1NonUnique(Mc_T{1}.TemperatureK,...
            Mc_T{1}.InvChi_molOe_Per_emu,SusceptibilityData.Temperature);

        % Define susceptibility in molar units for use in fitting HT B20
        % value.
        SusceptibilityData.FullTemperature_AB = Mab_T{1}.TemperatureK;
        SusceptibilityData.FullTemperature_C = Mc_T{1}.TemperatureK;
        SusceptibilityData.ChiAB_cm3_Per_mol = (1./Mab_T{1}.InvChi_molOe_Per_emu).*4.*pi;
        SusceptibilityData.ChiC_cm3_Per_mol = (1./Mc_T{1}.InvChi_molOe_Per_emu).*4.*pi;
        SusceptibilityData.Inverse_ChiAB_Full_molOe_Per_emu= 1./Mab_T{1}.InvChi_molOe_Per_emu;
        SusceptibilityData.Inverse_ChiC_Full_molOe_Per_emu= 1./Mc_T{1}.InvChi_molOe_Per_emu;

        SusceptibilityData.ChiAB_meV_Per_T2 = (1./SusceptibilityData.Inverse_ChiAB).*...
            10000.*BohrMagInMevPerT./(AvogadrosNumber.*BohrMagInEMU);
        SusceptibilityData.ChiC_meV_Per_T2 = (1./SusceptibilityData.Inverse_ChiC).*...
            10000.*BohrMagInMevPerT./(AvogadrosNumber.*BohrMagInEMU);

        LowerInterpLimit = max([min(SusceptibilityData.FullTemperature_AB),min(SusceptibilityData.FullTemperature_C)]);
        UpperInterpLimit = min([max(SusceptibilityData.FullTemperature_AB),max(SusceptibilityData.FullTemperature_C)]);
        SusceptibilityData.InterpTemp = linspace(LowerInterpLimit,UpperInterpLimit,2000);
        SusceptibilityData.Interp_Inverse_ChiAB_cm3_Per_mol = Interp1NonUnique(SusceptibilityData.FullTemperature_AB,...
            1./SusceptibilityData.ChiAB_cm3_Per_mol,SusceptibilityData.InterpTemp);
        SusceptibilityData.Interp_Inverse_ChiC_cm3_Per_mol = Interp1NonUnique(SusceptibilityData.FullTemperature_C,...
            1./SusceptibilityData.ChiC_cm3_Per_mol,SusceptibilityData.InterpTemp);

    otherwise
        disp('Invalid SampleOption in GenerateProcessParameters!');

end


end

function [outcell]=AddColumnsMagnetotropicData(incell,ProcessParameters)
% Internal function to add columns to magnetotropic data for the fitting.
% Adds a column to the loaded datasets for the magnetotropic coefficient
% calculated from the frequency. Also adds an interpolated magnetotropic
% coefficient and interpolated field values for the fitting.

outcell = incell;
for i=1:length(incell)
    outcell(i).k = outcell(i).C_frequency./(1e3.*ProcessParameters.Amp/(outcell(i).zf_frequency./1000));
    outcell(i).InterpBounds = [min(outcell(i).C_Field),max(outcell(i).C_Field)];
    outcell(i).InterpFields = linspace(outcell(i).InterpBounds(1),outcell(i).InterpBounds(2),20);
    outcell(i).k_interp = Interp1NonUnique(outcell(i).C_Field,outcell(i).k,...
        outcell(i).InterpFields);
end
end

function [Bmn,Jxx_Parameter,Jzz_Parameter,chi0_Parameter] = GenerateParameters(InputParameters)
% This internal function takes InputParameters and determines the validity
% of the input, and either gives zero for all the parameters if none are
% given, or appropriately assigns parameters given a valid InputParameter
% vector.

if isstring(InputParameters)
    Bmn           = zeros(1,6);
    Jxx_Parameter = 0;
    Jzz_Parameter = 0;
else
    Bmn = InputParameters(1:6);
    Jxx_Parameter = InputParameters(7);
    Jzz_Parameter = InputParameters(8);
    chi0_Parameter= InputParameters(9);
end


end

function [CEF_Hamiltonian] = GenerateCEFHamiltonian()
% This internal function returns an anonymous function  for the CEF
% Hamiltonian that takes as an input a set of six CEF parameters

% Pull relevant matrices for writing down the CEF Hamiltonian.
[ID,~,~,Jz,Jp,Jm] = GenerateSpinMatrices();

% System has total angular momentum 7/2, We'll need this factor and j(j+1)
% multiple times in the construction of the CEF Hamiltonian, so we define
% them here.
J = 7/2;
X = J*(J+1);

% Specific combination of raising and lowering operators that shows up in
% HCEF more than once. Written here for convenience.
A = Jp*Jp*Jp + Jm*Jm*Jm;


% The matrices that make up the CEF Hamiltonian.
O20 = 3*Jz*Jz - X*ID;
O40 = 35*power(Jz,4) - (30*X - 25)*Jz*Jz + (3*X*X - 6*X)*ID;
O60 = 231*power(Jz,6) - (315*X-735)*power(Jz,4) +...
    (105*X*X - 525*X +294)*power(Jz,2) -(5*X*X*X + 40*X*X -60*X)*ID;
O43 = (1/4)*( (A)*Jz + Jz*(A) );
O63 = (1/4)*( A*(11*power(Jz,3) - (3*X + 59)*Jz ) +...
    (11*power(Jz,3) -(3*X + 59)*Jz)*A );
O66 = (1/2)*(Jp*Jp*Jp*Jp*Jp*Jp + Jm*Jm*Jm*Jm*Jm*Jm);

% To scale the Stevens operators for effective fitting, we use the 'fac'
% vector to divide the relevant O matrices. When reporting the
% coefficients, this is eventually removed and is purely for computational
% purposes.
fac = [1,10,1,1e3,1e2,1e2];

% Define the anonymous function 'CEF_Hamiltonian' that depends on the
% matrices defined above and a set of 6 element vector Bmn of Steven's
% operator coefficients.
CEF_Hamiltonian = @(Bmn) Bmn(1)*O20/fac(1) + Bmn(2)*O40/fac(2) + Bmn(3)*O43/fac(3)...
    + Bmn(4)*O60/fac(4) + Bmn(5)*O63/fac(5) + Bmn(6)*O66/fac(6);


end

function [TotalError] = CumulativeError(Bmn_Parameter,Jxx_Parameter,Jzz_Parameter,chi0_Parameter,CEF_Hamiltonian,ID,Jx,Jz,MagnetotropicData,SusceptibilityData,muB,kB,gJ,A1,C0,ProcessParameters)
% This internal function calculates a residual error for the fitting of the
% magnetotropic coefficient. InputParameters should be a vector containing
% the six Steven's operator coefficients followed by Jxx and Jzz. The
% Steven's operator coefficients are in units of meV while Jxx and Jzz are
% in K. Lastly, it contains chi0 -- a chi shift.

% Declare initial versions of parameters/Hamiltonian before
CEF_Hamiltonian_Initial = CEF_Hamiltonian(Bmn_Parameter);
CEF_Hamiltonian_Initial = (CEF_Hamiltonian_Initial+...
    CEF_Hamiltonian_Initial.')/2; % The Hamiltonian picks up a small machine precision error. Round it away.
Jxx_Initial = Jxx_Parameter;
Jzz_Initial = Jzz_Parameter;
Chi0_Initial= chi0_Parameter;

% Find the eigenvalues and eigenvectors of the starting Hamiltonian. This
% was originally set to  eig(CEF_Hamiltonian_Initial + Jz*1e-10) but I have
% no idea why -- I suspect this was to lift some degeneracy or cure some
% rounding problem, but I'm not sure. NEEDS TO BE CHECKED 1/26/22.
[~,EigenValues_Initial] = eig(CEF_Hamiltonian_Initial);
[EnergyValues_Initial,~] = sort(diag(EigenValues_Initial));
ZeroEnergy = EnergyValues_Initial(1);

TotalError = 0; % Declare total error as zero to start.
FitRange = ProcessParameters.FitRange; % Which datasets are going to be fit? This is changed in 'LoadSampleData'
for q = FitRange
    % Magnetotropic Data [parameters].
    Temperature    = MagnetotropicData{1}(q).Temperature;
    kab_Data       = MagnetotropicData{1}(q).k_interp;
    FieldValues_ab = MagnetotropicData{1}(q).InterpFields;
    kc_Data        = MagnetotropicData{2}(q).k_interp;
    FieldValues_c  = MagnetotropicData{2}(q).InterpFields;

    % Calculated Magnetotropic coefficients at fixed temperature and
    % several field values (FieldValues_ab/c).
    kab_Calc = Magnetotropic_AB_vs_H(CEF_Hamiltonian_Initial,...
        ZeroEnergy,ID,Jx,Jz,Jxx_Initial,Jzz_Initial,Temperature,...
        FieldValues_ab,muB,kB,gJ,A1);
    kc_Calc = Magnetotropic_C_vs_H(CEF_Hamiltonian_Initial,...
        ZeroEnergy,ID,Jx,Jz,Jxx_Initial,Jzz_Initial,Temperature,...
        FieldValues_c,muB,kB,gJ,A1);

    % It seems like Chris was not previously doing a normalized residual
    % which would weight larger frequency shift data as more important to
    % fit. Maybe there was a reason for this, but perhaps consider changing
    % this to be mean residual normalized to get a least squares like sum.
    % In it's current form (1/27/22) it remains un-normalized.
    TotalError = TotalError + sum(power(kab_Data-kab_Calc,2)) +...
        sum(power(kc_Data-kc_Calc,2));
end

% Include an error associated with a fit of the susceptibility as well.
% First define Inverse chi data from input.
InverseChiAB_Data = SusceptibilityData.Inverse_ChiAB;
InverseChiC_Data  = SusceptibilityData.Inverse_ChiC;
ChiTemperature    = SusceptibilityData.Temperature;

% Calculate inverse chi. Include an optional vertical offset chi0.
InverseChiAB_Calc = InverseSusceptibility_AB_vs_T(CEF_Hamiltonian_Initial,...
    ZeroEnergy,ID,Jx,Jxx_Initial,ChiTemperature,muB,kB,gJ,C0);
InverseChiC_Calc  = InverseSusceptibility_C_vs_T(CEF_Hamiltonian_Initial,...
    ZeroEnergy,ID,Jz,Jzz_Initial,ChiTemperature,muB,kB,gJ,C0);
InverseChiAB_Calc =  1./(1./InverseChiAB_Calc + Chi0_Initial);
InverseChiC_Calc =  1./(1./InverseChiC_Calc + Chi0_Initial);

% Calculate the new total error including the chi residual.
TotalError = TotalError + ProcessParameters.ChiFittingWeight.*...
    (  sum(power(InverseChiAB_Data - InverseChiAB_Calc,2)) +...
    sum(power(InverseChiC_Data - InverseChiC_Calc,2))  );
end

function [TotalError] = CumulativeError_BoundEigenvalues(Bmn_Parameter,Jxx_Parameter,Jzz_Parameter,chi0_Parameter,CEF_Hamiltonian,ID,Jx,Jz,MagnetotropicData,SusceptibilityData,muB,kB,gJ,A1,C0,ProcessParameters)
% This internal function calculates a residual associated with the
% simultaneous fitting of susceptibility and magnetotropic data.
% Additionally, this function introduces a penalty for the eigenvalue
% energies deviating from the values obtained from Raman spectroscopy.


% Declare initial versions of parameters/Hamiltonian before
CEF_Hamiltonian_Initial = CEF_Hamiltonian(Bmn_Parameter);
CEF_Hamiltonian_Initial = (CEF_Hamiltonian_Initial+...
    CEF_Hamiltonian_Initial.')/2; % The Hamiltonian picks up a small machine precision error. Round it away.
Jxx_Initial = Jxx_Parameter;
Jzz_Initial = Jzz_Parameter;
Chi0_Initial= chi0_Parameter;


% Find the eigenvalues and eigenvectors of the starting Hamiltonian. This
% was originally set to  eig(CEF_Hamiltonian_Initial + Jz*1e-10) but I have
% no idea why -- I suspect this was to lift some degeneracy or cure some
% rounding problem, but I'm not sure. NEEDS TO BE CHECKED 1/26/22.
[~,EigenValues_Initial] = eig(CEF_Hamiltonian_Initial);
[EnergyValues_Initial,~] = sort(diag(EigenValues_Initial));
ZeroEnergy = EnergyValues_Initial(1);

TotalError = 0; % Declare total error as zero to start.
FitRange = ProcessParameters.FitRange; % Which datasets are going to be fit? This is changed in 'LoadSampleData'
for q = FitRange
    % Magnetotropic Data [parameters].
    Temperature    = MagnetotropicData{1}(q).Temperature;
    kab_Data       = MagnetotropicData{1}(q).k_interp;
    FieldValues_ab = MagnetotropicData{1}(q).InterpFields;
    kc_Data        = MagnetotropicData{2}(q).k_interp;
    FieldValues_c  = MagnetotropicData{2}(q).InterpFields;

    % Calculated Magnetotropic coefficients at fixed temperature and
    % several field values (FieldValues_ab/c).
    kab_Calc = Magnetotropic_AB_vs_H(CEF_Hamiltonian_Initial,...
        ZeroEnergy,ID,Jx,Jz,Jxx_Initial,Jzz_Initial,Temperature,...
        FieldValues_ab,muB,kB,gJ,A1);
    kc_Calc = Magnetotropic_C_vs_H(CEF_Hamiltonian_Initial,...
        ZeroEnergy,ID,Jx,Jz,Jxx_Initial,Jzz_Initial,Temperature,...
        FieldValues_c,muB,kB,gJ,A1);

    % Normalized as of 2/4/2022.
    TotalError = TotalError + (sum(power(kab_Data-kab_Calc,2))./...
        sum( power(kab_Data-mean(kab_Data),2) )+...
        sum(power(kc_Data-kc_Calc,2))./sum(power(kc_Data-mean(kc_Data),2)));
end

% Include an error associated with a fit of the susceptibility as well.
% First define Inverse chi data from input.
InverseChiAB_Data = SusceptibilityData.Inverse_ChiAB;
InverseChiC_Data  = SusceptibilityData.Inverse_ChiC;
ChiTemperature    = SusceptibilityData.Temperature;

% Calculate inverse chi. Include an optional vertical offset chi0.
InverseChiAB_Calc = InverseSusceptibility_AB_vs_T(CEF_Hamiltonian_Initial,...
    ZeroEnergy,ID,Jx,Jxx_Initial,ChiTemperature,muB,kB,gJ,C0);
InverseChiC_Calc  = InverseSusceptibility_C_vs_T(CEF_Hamiltonian_Initial,...
    ZeroEnergy,ID,Jz,Jzz_Initial,ChiTemperature,muB,kB,gJ,C0);
InverseChiAB_Calc =  1./(1./InverseChiAB_Calc + Chi0_Initial);
InverseChiC_Calc =  1./(1./InverseChiC_Calc + Chi0_Initial);

% Calculate the new total error including the chi residual, normalized as
% of 2/4/2022.
TotalError = TotalError + ProcessParameters.ChiFittingWeight.*...
    (  sum(power(InverseChiAB_Data - InverseChiAB_Calc,2))./...
    sum(power( InverseChiAB_Data - mean(InverseChiAB_Data),2)) +...
    sum(power(InverseChiC_Data - InverseChiC_Calc,2))./...
    sum(power( InverseChiC_Data - mean(InverseChiC_Data),2)) );

% Add a penalty for the eigenvalues being different than the values
% obtained from Raman.
UniqueEigenvalues = (EnergyValues_Initial(1:2:end)-ZeroEnergy)';
FixedEigenvalues  = ProcessParameters.Eigenenergies;
TotalError = TotalError + ...
    sum( ProcessParameters.RamanEnergyFittingWeight.*(FixedEigenvalues - UniqueEigenvalues).^2 )./...
    sum( (FixedEigenvalues-mean(FixedEigenvalues)).^2 );
end

function [TotalError] = CumulativeError_BoundEigenvalues_AndChiParams(Bmn_Parameter,Jxx_Parameter,Jzz_Parameter,chi0_Parameter,CEF_Hamiltonian,ID,Jx,Jz,MagnetotropicData,SusceptibilityData,muB,kB,gJ,A1,C0,ProcessParameters)
% This internal function calculates a residual associated with the
% simultaneous fitting of susceptibility and magnetotropic data.
% Additionally, this function introduces a penalty for the eigenvalue
% energies deviating from the values obtained from Raman spectroscopy.


% Declare initial versions of parameters/Hamiltonian before
CEF_Hamiltonian_Initial = CEF_Hamiltonian(Bmn_Parameter);
CEF_Hamiltonian_Initial = (CEF_Hamiltonian_Initial+...
    CEF_Hamiltonian_Initial.')/2; % The Hamiltonian picks up a small machine precision error. Round it away.
Jxx_Initial = Jxx_Parameter;
Jzz_Initial = Jzz_Parameter;
Chi0_Initial= chi0_Parameter;


% Find the eigenvalues and eigenvectors of the starting Hamiltonian. This
% was originally set to  eig(CEF_Hamiltonian_Initial + Jz*1e-10) but I have
% no idea why -- I suspect this was to lift some degeneracy or cure some
% rounding problem, but I'm not sure. NEEDS TO BE CHECKED 1/26/22.
% [~,EigenValues_Initial] = eig(CEF_Hamiltonian_Initial);
% [EnergyValues_Initial,~] = sort(diag(EigenValues_Initial));
% ZeroEnergy = EnergyValues_Initial(1);

[EigenVectors_Initial,EigenValues_Initial] = eig(CEF_Hamiltonian_Initial + Jz.*1e-10);
[EnergyValues_Initial,sinds] = sort(diag(EigenValues_Initial));
ZeroEnergy = EnergyValues_Initial(1);
EigenVectors_Initial = EigenVectors_Initial(:,sinds);

TotalError = 0; % Declare total error as zero to start.
FitRange = ProcessParameters.FitRange; % Which datasets are going to be fit? This is changed in 'LoadSampleData'
for q = FitRange
    % Magnetotropic Data [parameters].
    Temperature    = MagnetotropicData{1}(q).Temperature;
    kab_Data       = MagnetotropicData{1}(q).k_interp;
    FieldValues_ab = MagnetotropicData{1}(q).InterpFields;
    kc_Data        = MagnetotropicData{2}(q).k_interp;
    FieldValues_c  = MagnetotropicData{2}(q).InterpFields;

    % Calculated Magnetotropic coefficients at fixed temperature and
    % several field values (FieldValues_ab/c).
    kab_Calc = Magnetotropic_AB_vs_H(CEF_Hamiltonian_Initial,...
        ZeroEnergy,ID,Jx,Jz,Jxx_Initial,Jzz_Initial,Temperature,...
        FieldValues_ab,muB,kB,gJ,A1);
    kc_Calc = Magnetotropic_C_vs_H(CEF_Hamiltonian_Initial,...
        ZeroEnergy,ID,Jx,Jz,Jxx_Initial,Jzz_Initial,Temperature,...
        FieldValues_c,muB,kB,gJ,A1);

    % Normalized as of 2/4/2022.
    TotalError = TotalError + (sum(power(kab_Data-kab_Calc,2))./...
        sum( power(kab_Data-mean(kab_Data),2) )+...
        sum(power(kc_Data-kc_Calc,2))./sum(power(kc_Data-mean(kc_Data),2)));
end

% Include an error associated with a fit of the susceptibility as well.
% First define Inverse chi data from input.
InverseChiAB_Data = SusceptibilityData.Inverse_ChiAB;
InverseChiC_Data  = SusceptibilityData.Inverse_ChiC;
ChiTemperature    = SusceptibilityData.Temperature;

% Calculate inverse chi. Include an optional vertical offset chi0.
InverseChiAB_Calc = InverseSusceptibility_AB_vs_T(CEF_Hamiltonian_Initial,...
    ZeroEnergy,ID,Jx,Jxx_Initial,ChiTemperature,muB,kB,gJ,C0);
InverseChiC_Calc  = InverseSusceptibility_C_vs_T(CEF_Hamiltonian_Initial,...
    ZeroEnergy,ID,Jz,Jzz_Initial,ChiTemperature,muB,kB,gJ,C0);
InverseChiAB_Calc =  1./(1./InverseChiAB_Calc + Chi0_Initial);
InverseChiC_Calc =  1./(1./InverseChiC_Calc + Chi0_Initial);

% Calculate the new total error including the chi residual, normalized as
% of 2/4/2022.
TotalError = TotalError + ProcessParameters.ChiFittingWeight.*...
    (  sum(power(InverseChiAB_Data - InverseChiAB_Calc,2))./...
    sum(power( InverseChiAB_Data - mean(InverseChiAB_Data),2)) +...
    sum(power(InverseChiC_Data - InverseChiC_Calc,2))./...
    sum(power( InverseChiC_Data - mean(InverseChiC_Data),2)) );

% Add a penalty for the eigenvalues being different than the values
% obtained from Raman.
UniqueEigenvalues = (EnergyValues_Initial(1:2:end)-ZeroEnergy)';
FixedEigenvalues  = ProcessParameters.Eigenenergies;
TotalError = TotalError + ...
    sum( ProcessParameters.RamanEnergyFittingWeight.*(FixedEigenvalues - UniqueEigenvalues).^2 )./...
    sum( (FixedEigenvalues-mean(FixedEigenvalues)).^2 );

% Add a penalty for the chi parameters (A0 and perturb coeff) being
% different than what is obtained from susceptibiltiy fitting.
EV1N = EigenVectors_Initial(:,1:2:8);
EV2N = EigenVectors_Initial(:,2:2:8);

AZ0_Calc = -power(muB*gJ,2).*((power(EV1N(:,1)'*Jz*EV1N(:,2),2) + power(EV1N(:,1)'*Jz*EV2N(:,2),2) )./(EnergyValues_Initial(3) - ZeroEnergy)...
    +(power(EV1N(:,1)'*Jz*EV1N(:,3),2) + power(EV1N(:,1)'*Jz*EV2N(:,3),2) )./(EnergyValues_Initial(5) - ZeroEnergy)...
    +(power(EV1N(:,1)'*Jz*EV1N(:,4),2) + power(EV1N(:,1)'*Jz*EV2N(:,4),2) )./(EnergyValues_Initial(7) - ZeroEnergy));
%
AX0_Calc = -power(muB*gJ,2).*((power(EV1N(:,1)'*Jx*EV1N(:,2),2) + power(EV1N(:,1)'*Jx*EV2N(:,2),2) )./(EnergyValues_Initial(3) - ZeroEnergy)...
    +(power(EV1N(:,1)'*Jx*EV1N(:,3),2) + power(EV1N(:,1)'*Jx*EV2N(:,3),2) )./(EnergyValues_Initial(5) - ZeroEnergy)...
    +(power(EV1N(:,1)'*Jx*EV1N(:,4),2) + power(EV1N(:,1)'*Jx*EV2N(:,4),2) )./(EnergyValues_Initial(7) - ZeroEnergy));

A0_AB_Calc = power(EV1N(:,1)'*Jx*EV2N(:,1),2);
A0_C_Calc = power(EV1N(:,1)'*Jz*EV1N(:,1),2);

Chi_PerturbationTerms = ProcessParameters.Chi_PerturbationTerms;
Chi_A0Vals = ProcessParameters.Chi_A0Vals;

TotalError = TotalError + ...
    sum( ProcessParameters.ChiFitResultsWeight.*...
    (Chi_PerturbationTerms - [AZ0_Calc,AX0_Calc]).^2 )./...
    sum( (Chi_PerturbationTerms-mean(Chi_PerturbationTerms)).^2 );

TotalError = TotalError + ...
    sum( ProcessParameters.ChiFitResultsWeight.*...
    (Chi_A0Vals - [A0_C_Calc,A0_AB_Calc]).^2 )./...
    sum( (Chi_A0Vals-mean(Chi_A0Vals)).^2 );



end



function [AX0,AZ0,A0_AB,A0_C] = Calculate_AZ(Bmn_Parameter,Jxx_Parameter,Jzz_Parameter,chi0_Parameter,CEF_Hamiltonian,ID,Jx,Jz,muB,gJ)
% This internal function calculates a residual associated with the
% simultaneous fitting of susceptibility and magnetotropic data.
% Additionally, this function introduces a penalty for the eigenvalue
% energies deviating from the values obtained from Raman spectroscopy.


% Declare initial versions of parameters/Hamiltonian before
CEF_Hamiltonian_Initial = CEF_Hamiltonian(Bmn_Parameter);
CEF_Hamiltonian_Initial = (CEF_Hamiltonian_Initial+...
    CEF_Hamiltonian_Initial.')/2; % The Hamiltonian picks up a small machine precision error. Round it away.
Jxx_Initial = Jxx_Parameter;
Jzz_Initial = Jzz_Parameter;
Chi0_Initial= chi0_Parameter;


% Find the eigenvalues and eigenvectors of the starting Hamiltonian. This
% was originally set to  eig(CEF_Hamiltonian_Initial + Jz*1e-10) but I have
% no idea why -- I suspect this was to lift some degeneracy or cure some
% rounding problem, but I'm not sure. NEEDS TO BE CHECKED 1/26/22.
[EigenVectors_Initial,EigenValues_Initial] = eig(CEF_Hamiltonian_Initial + Jz.*1e-10);
[EnergyValues_Initial,sinds] = sort(diag(EigenValues_Initial));
ZeroEnergy = EnergyValues_Initial(1);
EigenVectors_Initial = EigenVectors_Initial(:,sinds);

% % % Code from Chris Old Version, checked that it produced matching results
% % % with new code...3/8/2022
% % Ev1 = EnergyValues_Initial;
% % En = diag(EigenValues_Initial);
% % P1 = EigenVectors_Initial;
% % E1 = ZeroEnergy;
% % for n=1:4
% %     ind1 = find(En==Ev1(2*(n-1) + 1));
% %     ind2 = find(En==Ev1(2*n));
% %
% %     ev1{n} = P1(:,ind1(1));
% %     ev2{n} = P1(:,ind2(end));
% % end
% % AZ0_old = -power(muB*gJ,2).*((power(ev1{1}'*Jz*ev1{2},2) + power(ev1{1}'*Jz*ev2{2},2) )./(Ev1(3) - E1)...
% %     +(power(ev1{1}'*Jz*ev1{3},2) + power(ev1{1}'*Jz*ev2{3},2) )./(Ev1(5) - E1)...
% %     +(power(ev1{1}'*Jz*ev1{4},2) + power(ev1{1}'*Jz*ev2{4},2) )./(Ev1(7) - E1));
% % AX0_old = -power(muB*gJ,2).*((power(ev1{1}'*Jx*ev1{2},2) + power(ev1{1}'*Jx*ev2{2},2) )./(Ev1(3) - E1)...
% %     +(power(ev1{1}'*Jx*ev1{3},2) + power(ev1{1}'*Jx*ev2{3},2) )./(Ev1(5) - E1)...
% %     +(power(ev1{1}'*Jx*ev1{4},2) + power(ev1{1}'*Jx*ev2{4},2) )./(Ev1(7) - E1));
% %

EV1N = EigenVectors_Initial(:,1:2:8);
EV2N = EigenVectors_Initial(:,2:2:8);

AZ0 = -power(muB*gJ,2).*((power(EV1N(:,1)'*Jz*EV1N(:,2),2) + power(EV1N(:,1)'*Jz*EV2N(:,2),2) )./(EnergyValues_Initial(3) - ZeroEnergy)...
    +(power(EV1N(:,1)'*Jz*EV1N(:,3),2) + power(EV1N(:,1)'*Jz*EV2N(:,3),2) )./(EnergyValues_Initial(5) - ZeroEnergy)...
    +(power(EV1N(:,1)'*Jz*EV1N(:,4),2) + power(EV1N(:,1)'*Jz*EV2N(:,4),2) )./(EnergyValues_Initial(7) - ZeroEnergy));
%
AX0 = -power(muB*gJ,2).*((power(EV1N(:,1)'*Jx*EV1N(:,2),2) + power(EV1N(:,1)'*Jx*EV2N(:,2),2) )./(EnergyValues_Initial(3) - ZeroEnergy)...
    +(power(EV1N(:,1)'*Jx*EV1N(:,3),2) + power(EV1N(:,1)'*Jx*EV2N(:,3),2) )./(EnergyValues_Initial(5) - ZeroEnergy)...
    +(power(EV1N(:,1)'*Jx*EV1N(:,4),2) + power(EV1N(:,1)'*Jx*EV2N(:,4),2) )./(EnergyValues_Initial(7) - ZeroEnergy));

A0_AB = power(EV1N(:,1)'*Jx*EV2N(:,1),2);
A0_C = power(EV1N(:,1)'*Jz*EV1N(:,1),2);
end


function Magnetotropic_Out = Magnetotropic_AB_vs_H(CEF_Hamiltonian_Matrix,ZeroEnergy,ID,Jx,Jz,Jxx_Parameter,Jzz_Parameter,Temperature,Field,muB,kB,gJ,A1)
% Internal function for calculating the magnetotropic coefficient with H
% along ab. Takes as an input the CEF Hamiltonian, the zeroeth eigenenergy,
% the Jx and Jz angular momentum matrices, the temperature, the applied
% field, and some constants.

%
Beta = 1./(kB.*Temperature); % Stat Mech Beta 1/k_bT.
hs = 0.01; %Small field in c direction for calculation?
LengthFields = length(Field); % Loop counter for different applied fields.

% Declare these forr storage purposes.
mx   = zeros(1,LengthFields);
mzT0 = zeros(1,LengthFields);
mzT1 = zeros(1,LengthFields);

% Loop calculations over all fields.
for q=1:LengthFields
    % Calculate Hamiltonian with an included Zeeman term in the x and z
    % directions.
    CurrentFieldValue = Field(q);
    Hamiltonian_WithZeeman = CEF_Hamiltonian_Matrix - ...
        gJ*muB*CurrentFieldValue*Jx;

    % Calculate x-magnetization.
    mx(q) = trace( Jx*expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)) )./...
        trace( expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)));

    % Calculate z-magnetization.
    mzT0(q) = trace( Jz*expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)) )./...
        trace( expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)));

    % Recalculate Hamiltonian with a small applied field in z-direction,
    % then recalculate mz.
    Hamiltonian_WithZeeman = Hamiltonian_WithZeeman - gJ*muB*hs*Jz;
    mzT1(q) = trace( Jz*expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)) )./...
        trace( expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)));
end


X1 = Field + Jxx_Parameter*A1*mx; X2 = power(X1,2);
chiT = (mzT1-mzT0)./(hs + Jzz_Parameter*A1*(mzT1 - mzT0));
Magnetotropic_Out = Interp1NonUnique([0 X1],...
    [0 gJ*muB*(X1.*mx - X2.*chiT)],Field);
end

function Magnetotropic_Out = Magnetotropic_C_vs_H(CEF_Hamiltonian_Matrix,ZeroEnergy,ID,Jx,Jz,Jxx_Parameter,Jzz_Parameter,Temperature,Field,muB,kB,gJ,A1)
% Internal function for calculating the magnetotropic coefficient with H
% along ab. Takes as an input the CEF Hamiltonian, the zeroeth eigenenergy,
% the Jxx and Jzz angular momentum matrices, the temperature, the applied
% field, and some constants.

%
Beta = 1./(kB.*Temperature); % Stat Mech Beta 1/k_bT.
hs = 0.01; %Small field in c direction for calculation?
LengthFields = length(Field); % Loop counter for different applied fields.

% Declare these forr storage purposes.
mz   = zeros(1,LengthFields);
mxT0 = zeros(1,LengthFields);
mxT1 = zeros(1,LengthFields);

% Loop calculations over all fields.
for q=1:LengthFields
    % Calculate Hamiltonian with an included Zeeman term in the x and z
    % directions.
    CurrentFieldValue = Field(q);
    Hamiltonian_WithZeeman = CEF_Hamiltonian_Matrix - ...
        gJ*muB*CurrentFieldValue*Jz;

    % Calculate z-magnetization.
    mz(q) = trace( Jz*expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)) )./...
        trace( expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)));

    % Calculate x-magnetization
    mxT0(q) = trace( Jx*expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)) )./...
        trace( expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)));

    % Calculate x magnetization for small field applied in x direction.
    Hamiltonian_WithZeeman = Hamiltonian_WithZeeman - gJ*muB*hs*Jx;
    mxT1(q) = trace( Jx*expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)) )./...
        trace( expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)));

end

% Calculate magnetotropic coefficient.
X1 = Field + Jzz_Parameter*A1*mz; X2 = power(X1,2);
chiT = (mxT1-mxT0)./(hs + Jxx_Parameter*A1*(mxT1 - mxT0));
Magnetotropic_Out = Interp1NonUnique([0 X1],...
    [0 gJ*muB*(X1.*mz - X2.*chiT)], Field);
end

function Susceptibility_Out = InverseSusceptibility_AB_vs_T(CEF_Hamiltonian_Matrix,ZeroEnergy,ID,Jx,Jxx_Parameter,Temperature,muB,kB,gJ,C0)
% Internal function to calculate the susceptibility in the AB direction
% at several temperatures.

% The susceptibility is defined at zero applied field. Apply a 'small'
% field for the sake of calculation.Calculate Hamiltonian.
hs = 0.01;
Hamiltonian_WithZeeman = CEF_Hamiltonian_Matrix - ...
    gJ*muB*hs*Jx;
LengthTemperature = length(Temperature); % Number of temperatures for loop.

Magnetization = zeros(1,LengthTemperature); % Declare magnetization.

% Calculate magnetization at desired temperatures.
for q=1:LengthTemperature
    CurrentTemperature = Temperature(q); % Set the current temperature.
    Beta = 1./(kB.*CurrentTemperature); % Thermo beta=1/kB*T

    % Calculate the magnetization.
    Magnetization(q) = trace( Jx*expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)) )./...
        trace( expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)));
end

Susceptibility_Out = C0*((gJ*muB/kB)*hs./Magnetization + 6*Jxx_Parameter);
end

function Susceptibility_Out = InverseSusceptibility_C_vs_T(CEF_Hamiltonian_Matrix,ZeroEnergy,ID,Jz,Jzz_Parameter,Temperature,muB,kB,gJ,C0)
% Internal function to calculate the susceptibility in the AB direction
% at several temperatures.

% The susceptibility is defined at zero applied field. Apply a 'small'
% field for the sake of calculation.Calculate Hamiltonian.
hs = 0.01;
Hamiltonian_WithZeeman = CEF_Hamiltonian_Matrix - ...
    gJ*muB*hs*Jz;
LengthTemperature = length(Temperature); % Number of temperatures for loop.

Magnetization = zeros(1,LengthTemperature); % Declare magnetization.

% Calculate magnetization at desired temperatures.
for q=1:LengthTemperature
    CurrentTemperature = Temperature(q); % Set the current temperature.
    Beta = 1./(kB.*CurrentTemperature); % Thermo beta=1/kB*T

    % Calculate the magnetization.
    Magnetization(q) = trace( Jz*expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)) )./...
        trace( expm(-Beta*(Hamiltonian_WithZeeman - ZeroEnergy*ID)));
end

Susceptibility_Out = C0*((gJ*muB/kB)*hs./Magnetization + 6*Jzz_Parameter);
end




