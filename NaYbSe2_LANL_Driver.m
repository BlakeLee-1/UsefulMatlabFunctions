function [MC]=NaYbSe2_LANL_Driver(loadopt,plotopt,sortopt,newfigopt)
h2 = NewFigureCheck(newfigopt);

MC = LoadPFDat(loadopt);
MC = SortPFDat(MC,sortopt);
switch plotopt
    case 'All'
        figure(h2); hold on;
        for i=1:length(MC)
            plot(MC(1,i).Field,MC(1,i).frequency,'-','Color',...
                MC(1,i).UColor,'DisplayName',MC(1,i).LegendAngle)
        end
    case 'All Rough Delta'
        figure(h2); hold on;
        for i=1:length(MC)
            ainds = and(MC(1,i).Field>2,MC(1,i).Field<4);
            f0 = mean(MC(1,i).frequency(ainds)); 
            plot(MC(1,i).Field,MC(1,i).frequency-f0,'-','Color',...
                MC(1,i).UColor,'DisplayName',MC(1,i).LegendAngle)
        end       
        xlabel('\mu_0H [T]'); ylabel('\Deltaf [Hz]'); 
    case 'Down Rough Delta'
        figure(h2); hold on;
        lim = -.0001;
        for i=1:length(MC)
            kinds= diff(MC(1,i).Field)<lim;
            x = MC(1,i).Field(kinds);
            y = MC(1,i).frequency(kinds);
            [x,sinds] = sort(x);
            y = y(sinds);

            ainds = and(x>1,x<3);
            f0 = mean(y(ainds)); 
            NewColor = ColormapInterpolate_InBounds(MC(1,i).Temperature,[3.8,160]);
            plot(x,y-f0,'-','Color',...
                NewColor,'DisplayName',MC(1,i).LegendAngle);
        end       
        xlabel('\mu_0H [T]'); ylabel('\Deltaf [Hz]'); 
    case 'Down'
        lim = -.0001;
        figure(h2); hold on;
        for i=1:length(MC)
            kinds= diff(MC(1,i).Field)<lim;
            x = MC(1,i).Field(kinds);
            y = MC(1,i).frequency(kinds);
            [x,sinds] = sort(x);
            y = y(sinds);
            plot(x,y,'-','Color',MC(1,i).UColor,'DisplayName',MC(1,i).LegendAngle);
        end
    case 'Angle Dependence'
        cfields = 1:1:9;
        colorz = jet(length(cfields)); 
        for j=1:length(cfields)
            cfield = cfields(j);
            clear anglevals; clear freqvals; clear sinds; 
            for i=1:length(MC)
                anglevals(i) = MC(1,i).Angle;
                freqvals(i) = Interp1NonUnique(MC(1,i).Field,MC(1,i).frequency,cfield);
                [anglevals,sinds] = sort(anglevals);          
            end
                   hold on; plot(anglevals,freqvals(sinds),'-o','Color',...
                    colorz(j,:),'MarkerFaceColor',brc(colorz(j,:),.5))
          
        end
    case 'Angle Dependence HSquared'
        figure(h2); 
        for i=1:length(MC)
           kinds = and(MC(1,i).Field>1,MC(1,i).Field<10); 
           x = MC(1,i).Field(kinds).^2;
           y = MC(1,i).frequency(kinds);
           pfit = polyfit(x,y,1);
           y0(i) = pfit(2); 
           mval(i) = pfit(1); 
           anglevals(i) = MC(1,i).Angle;
        end
        [anglevals,svinds] = sort(anglevals); 
        mval = mval(svinds); 
        
    otherwise
        disp('Invalid plotopt!');
end

end


function h2 = NewFigureCheck(newfigopt)
if strcmp('New',newfigopt)
    h2 = figure;
else
    h2 = gcf;
end
end

function MC = LoadPFDat(loadopt)
switch loadopt
    case 'Phi1'
        fr = 'E:\IanComputer\Documents\Physics\Minhyea Lee Research\LosAlamos_Nov2021\NaYbSe2\TDOFiles\';
        [DateID,ID,Angle,Temperature] = LU_RTM_Nov2021(loadopt);
    case 'AB T Dependence'
        fr = 'E:\IanComputer\Documents\Physics\Minhyea Lee Research\LosAlamos_Nov2021\NaYbSe2\Phi1_AB\';
        [DateID,ID,Angle,Temperature] = LU_RTM_Nov2021(loadopt);

    otherwise
        disp('Invalid loadopt in LoadPFDat!');
end
cd(fr);
DStore = dir('*.dat');
colorz = jet(length(DStore));
for i=1:length(DStore)
    temp = readmatrix(DStore(i).name);
    MC(1,i).Field = temp(:,1);
    MC(1,i).frequency = temp(:,2);
    MC(1,i).ID = str2num(DStore(i).name(2:4));
    MC(1,i).Date_ID = str2num(DStore(i).name(8:9));
    MC(1,i).Name = DStore(i).name;
    MC(1,i).Angle = Angle(and(ID==MC(1,i).ID,DateID==MC(1,i).Date_ID));
    MC(1,i).Temperature = Temperature(and(ID==MC(1,i).ID,DateID==MC(1,i).Date_ID));
    MC(1,i).LegendAngle = ['ID ',num2str(MC(1,i).ID),': ',num2str(MC(1,i).Temperature),' K; ',num2str(MC(1,i).Angle),' deg'];
    %     MC(1,i).UColor = ColormapInterpolate_InBounds(MC(1,i).Angle,[70,200]);
    MC(1,i).UColor = colorz(i,:);
end
end

function [DateID,ID,Angle,Temperature] = LU_RTM_Nov2021(loadopt)
switch loadopt
    case 'Phi1'
        DateID = [repmat(1,1,28),2,2,2,2,2,2];
        ID =    [5,7,9,11,13,14,15,16,17,18,19,20,21,22,23,24,26,27,28,29,30,31,32,34,35,36,37,50,...
            15,17,32,34,97,98];
        Temperature = [repmat(4,1,25),10,10,20,40,40,80,80,160,160];
        Angle = [90,90,90,100.13,79.86,69.49,59.46,49.23,38.869,28.478,15.96,.83,...
            -9.29,-19.5498,-29.75,-.965,-39.62,-49.003,-59.203,-68.973,...
            -78.98,-89.28,-99.45,-85.3,-85.3,-85.8,-85.8,-85.8,...
            -86.57,-86.45,-86.25,-86.15,-85.06,-84.96];
    case 'AB T Dependence'
        DateID = [1,1,1,2,2,2];
        ID = [35,37,50,17,34,98];
        Temperature = [4,10,20,40,80,160]; 
        Angle = [-85.3,-85.7,-85.7,-86.5,-86.15,-84.97]; 
    otherwise
        disp('Invalid loadopt in LU_RTM_Nov2021!');
end
end

function outcell = SortPFDat(incell,sortopt)
for i=1:length(incell)
   sortvals(i) = mean(incell(1,i).(sortopt));
end
[~,sinds] = sort(sortvals); 
outcell = incell(1,sinds); 
end


function OutColormap = ColormapInterpolate_InBounds(InFields,InBounds)

CMAP=pmkmp(256,'Swtth');
% CMAP = varycolor(256); 
% CMAP = jet(256); 
InterpDensity = 256;
% Trange = logspace(-3,.85,256);
% Trange = linspace(-.305,1.26,256);
Trange = linspace(InBounds(1),InBounds(2),InterpDensity);
OutSize = length(InFields);
for i=1:OutSize
    try
        [~,ind] = min(abs(Trange-(InFields(i))));
%         OutColormap(i,:) = drc(CMAP(ind,:),.1); 
        OutColormap(i,:)= CMAP(ind,:);
    catch
        OutColormap(i,:) = [0 0 0];
    end
end

end