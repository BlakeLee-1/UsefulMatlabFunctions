function plotPaper
M = LoadData()
plotboi(M)
end


function [outcell] = LoadData()
loadopt = 'Phi1';
switch loadopt
    case 'Phi1'
        fr = ['/home/blake/Documents/MATLAB/MPMS Data/T'];
        Phi_Angle = 0;
        headerlines = 1;
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
    outcell(i).x = temp(headerlines:end,1);
    outcell(i).y = temp(headerlines:end,2);
    
end
end

function plotboi(incell)
plot(incell(1).x,incell(1).y.*10^4*52.6/38.8, 'ko')
end