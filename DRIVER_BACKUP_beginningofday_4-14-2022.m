% % BACKUP - beginning of day 4/14/2022
% % Sarah C. Jones
% % 3/29/2022
% % Function to analyze data and output different plots
% 
% function [MS] = Sarah_KVO_Driver(Calibration, PlotOpt)
%     MS = LoadData();
%     MS = ProcessData(.2,MS);  % Process data with a tolerance of 20%
%     
%     switch Calibration
%         case 'Recalculate'
%             % New Hot & Cold Cernox Calibration
%             for i = 1:length(MS)
%                 MS(i).BathTempCalibrate = MS(i).BathTemp(:,1:MS(i).SampPerStep) % New field: first step Bath Temp
%                 MS(i).HotResistanceCalibrate = MS(i).HotRes(:,1:MS(i).SampPerStep);
%                 MS(i).ColdResistanceCalibrate = MS(i).ColdRes(:,1:MS(i).SampPerStep);
%             end
%             
%             SarahsPlot(MS,'BathTempCalibrate','HotResistanceCalibrate','r','Hot Cernox')
%             xlabel('T [K]'); ylabel('Resistance [\Omega]');
%             SarahsPlot(MS,'BathTempCalibrate','ColdResistanceCalibrate','b','Cold Cernox')
%             xlabel('T [K]'); ylabel('Resistance [\Omega]');
%             
%             deg = 7;
%             [~,sinds] = sort([MS.BathTempCalibrate]);
%             y = [MS.BathTempCalibrate];
%             x = [MS.HotResistanceCalibrate];
%             x = x(sinds); y = y(sinds);
%             HotFit = RtoT_cal_inputRT(x,y,'fitRvT',deg,'Yes','Polynomial',[]);
%             x = [MS.ColdResistanceCalibrate];
%             x = x(sinds);
%             ColdFit = RtoT_cal_inputRT(x,y,'fitRvT',deg,'Yes','Polynomial',[]);
%             
%         case 'Calculated'
%             % Applying new calibration
%             temp = load('C:\Users\LeeLabLaptop\Documents\Sarah\KVO_Nitrogen\KVO_InSitu_MN_CalibrationCurves_4_11_22')
%             CXCell_MN = temp.CXCell_MN
%             for i = 1:length(MS)
%                 MS(i).HotTemp = CXCell_MN{1}(MS(i).HotRes);
%                 MS(i).ColdTemp = CXCell_MN{2}(MS(i).ColdRes);
%             end
%             
%         case 'None'
%             % do nothing
%             
%         otherwise
%             disp('Invalid Calibration')
%     end
%     
%     MS = Kappa(MS);
%     
%     switch PlotOpt
%         case 'Plot Raw and Calculated Data'
%             temp = load('C:\Users\LeeLabLaptop\Documents\Sarah\KVO_Nitrogen\KVO_InSitu_MN_CalibrationCurves_4_11_22');
%             CXCell_MN = temp.CXCell_MN;
%             SarahPlot_RawAndCal_Kappa(CXCell_MN);
%             
%         case 'Kappa v T' 
%             xarg = 'TotalAverageTemp';
%             yarg = 'Kappa';
%             MyColor = [0 0 1]; 
%             PlotTitle = '\kappa vs Temperature' 
%             SarahsPlot(MS,xarg,yarg,MyColor,PlotTitle);
%             xlabel('T [K]'); ylabel('\kappa [W/m\cdotK]');
%             
%         case 'Delta T vs Power'
%            figure('name','Delta T vs Power');
%            colorz = jet(length(MS));
%            displayname = num2str(MS(i).TotalAverageTemp);
%            for i = 1:25:length(MS)
%                plot(MS(i).HeaterPowerStepsCut,MS(i).DeltaTStepsCut,'o','color',colorz(i,:),...  %% Fix this to use calibrated CX temp vals!!!
%                    'DisplayName',displayname); hold on;
%            end
%            Label_Plot('Power','W','Delta T','K');
%            MyColorbar = colorbar;                   % Add color bar
%            colormap(jet(256));                      % Change colormap to match plot color scheme.
%            MyColobar.Limits = [80,290];             % Set color bar values to min and max temp
%            MyAxes = gca;
%            MyAxes.CLim = [80,290];
%            %MyColorbar.Ticks = [80:20:290];
%            
%         case 'None'
%            % No output
%            
%         otherwise 
%             disp('Invalid PlotOpt'); 
%     end
%     
% end
% 
% function [MS] = LoadData()
%     fileroot = 'C:\Users\LeeLabLaptop\Documents\Sarah\KVO_Nitrogen\TswpPractice';
%     cd(fileroot);
%     filedata = dir('*.mat');
%     for i = 1:length(filedata)
%        MS(i) = load(filedata(i).name);  % Load all data first before adding fields so they are all the same length
%     end
%     for i = 1:length(filedata)
%        MS(i).name = filedata(i).name;
%        MS(i).Fileroot = fileroot;
%        MS(i).DateStr = filedata(i).date;
%        MS(i).Datenum = filedata(i).datenum;
%     end
% end
% 
% function [MS] = ProcessData(Tolerance,MS)
%     for i = 1:length(MS)
%         HeaterCurrentSteps = 4;             % This should be modified to be a calculation
%         LogicalVector2 = logical([ones(1,HeaterCurrentSteps*MS(i).SampPerStep),zeros(1,MS(i).SampPerStep)]);                                  % Logic vector with 0's for I = 0
%         MS(i).DeltaTStepsOnly = MS(i).DeltaT(LogicalVector2);                    % Eliminate all Delta T where I = 0
%         MS(i).HeaterPowerStepsOnly = MS(i).HeaterPower(LogicalVector2);          % Eliminate all Power where I = 0
%         SingleStepLogical = logical([zeros(1,MS(i).SampPerStep*(1-Tolerance)),...
%             ones(1,MS(i).SampPerStep*(Tolerance))]);                            % Logical to cut each individual step
%         CutSteps = logical(repmat(SingleStepLogical,1,...
%             HeaterCurrentSteps));                                               % Repeating Logical to cut all steps
%         MS(i).DeltaTStepsCut = MS(i).DeltaT(CutSteps);                          % Cut Delta T to tolerance
%         MS(i).HeaterPowerStepsCut = MS(i).HeaterPower(CutSteps);                % Cut Power to tolerance
%     end
%     
%     % Getting sample avg temp and sorting it
%     for i = 1:length(MS)
%         MS(i).TotalAverageTemp = mean(MS(i).SampleAverageTemp);     % Create new field for single Total Average Temp per .mat data file
%     end
%     TotalAverageTemp = vertcat(MS.TotalAverageTemp);                % Create Total Average Temp vector from structure field
%     [Temp,Index] = sort(TotalAverageTemp);                           % Sort Total Average Temp vector and index values
%     MS = MS(Index);
%    
% end
%     
% function [MS] = Kappa(MS)
%     CrossSectArea = 1.2e-3 * 0.3e-3;
%     Length = 3e-3;
%     SampleGeometry = Length/CrossSectArea;
%     for i = 1:length(MS)
%         CoeffValues = polyfit(MS(i).DeltaTStepsCut,MS(i).HeaterPowerStepsCut,1);
%         MS(i).Kappa = CoeffValues(1)*SampleGeometry;
%         MS(i).Intercept = CoeffValues(2);
%     end
% end
% 
% function SarahsPlot(MS,xarg,yarg,PlotColor,PlotTitle)
%     figure; hold on;
%     plot([MS.(xarg)],[MS.(yarg)],'.','Color',PlotColor);
%     title(PlotTitle);
% end