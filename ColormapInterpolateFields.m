% Ian Leahy
% 4/15/2021
% Colormap colors Uniform in temperature between 0 and 310K.

function OutColormap = ColormapInterpolateFields(InFields)

CMAP=pmkmp(256,'Swtth');
InterpDensity = 256;
Trange = logspace(-3,.85,256);
OutSize = length(InFields);
for i=1:OutSize
    try
        [~,ind] = min(abs(Trange-(InFields(i))));
        OutColormap(i,:) = drc(CMAP(ind,:),.1);
    catch
        OutColormap(i,:) = [0 0 0];
    end
end

end