function yokoampset_obj(current,obj1)
%DC current output from Yokogawa GS200 source meter
%Sets output current level to current (amps)
fprintf(obj1,'SOURce:FUNCtion CURRent'); 
fprintf(obj1,['SOURce:RANGe ' sprintf('%0.9f',current)]); 
fprintf(obj1,['SOURce:LEVel ' sprintf('%0.9f',current)]); %Sets current to i
end