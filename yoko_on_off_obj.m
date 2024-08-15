function yoko_on_off_obj(obj1,SourceState)
%Turns on output from Yokogawa GS200 source meter
% SourceState can be 'ON' or 'OFF'.
fprintf(obj1,['OUTPut ',SourceState]); %Sets voltage mode
end