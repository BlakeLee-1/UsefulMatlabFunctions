% Ian Leahy
% 4/15/2021
% Colormap colors Uniform in temperature between 0 and 310K.

function OutColormap = ColormapInterpolate_InBounds(InFields,InBounds)

% CMAP=pmkmp(256,'Swtth');
% CMAP = varycolor(256);
CMAP = jet(256);
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