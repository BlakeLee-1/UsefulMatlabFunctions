%Ian Leahy
%Jan 10, 2016
%New MPMS load code -- much simpler than loadMPMSBatch.

%fileroot - Directory where all datafiles are located
%MeasurementType - RvsT
%                - RvsH
%                - MvsT
%                - MvsH
function outdata=SnatchMPMSData(filename,MeasurementType)

switch MeasurementType
    case 'RvsT'
        headerlines=24;
    case 'RvsH'
        headerlines=24;
    case 'MvsT'
        headerlines=31;
    case 'MvsH'
        headerlines=31;
    case 'MvsH Raw'
        headerlines=31;
    case 'IV'
        headerlines=24;
    otherwise
        error('Pick valid measurement type');
end

datahold=importdata(filename,',',headerlines);
try
    names=datahold.textdata(end,:);
catch
    if strcmp(datahold{end}(1:4),'Time')
        display(['It seems like the filename is empty: ',filename])
    else
        for i=headerlines:headerlines+20
            datahold=importdata(filename,',',i);
            try
                names=datahold.textdata(end,:);
                display(['Headerlines should be: ',num2str(headerlines)]);
                display(filename);
            catch
            end
        end
    end
end
if length(fieldnames(datahold))<3
    headerstring=datahold.textdata{end};
    names=strsplit(headerstring,',');
end


try
    for i=1:length(names)
        temphold=char(names(i));
        temphold(ismember(temphold,' ,.:;!')) = [];
        temphold(ismember(temphold,'(,.:;!')) = [];
        temphold(ismember(temphold,'),.:;!')) = [];
        temphold(ismember(temphold,'%,.:;!')) = [];
        temphold(ismember(temphold,'/,.:;!')) = [];
        temphold(ismember(temphold,'[,.:;!')) = [];
        temphold(ismember(temphold,'],.:;!')) = [];
        temphold(ismember(temphold,'-,.:;!')) = [];
        namescell{i}=temphold;
    end
    for i=1:length(names)
        try
        outdata.(namescell{i})=datahold.data(:,i);
        catch
        end
    end
catch
    outdata=[filename];
end

end