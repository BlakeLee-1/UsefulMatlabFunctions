%Load a directory of MPMS data
%Ian Leahy
%Jan 10,2016
% Input fileroot and MeasurementType and return outcell, a cell structure
% containing the magnetization data. 
%fileroot - Directory where all datafiles are located
%MeasurementType - RvsT
%                - RvsH
%                - MvsT
%                - MvsH
function outcell=DirectorySnatchMPMSData(fileroot,MeasurementType)

cd(fileroot);
DirectoryStore=dir;
for i=1:length(DirectoryStore)
    isdirlogical(i)=DirectoryStore(i).isdir;
    if ~isdirlogical(i)
        if strcmp(DirectoryStore(i).name(end-3:end),'.dat')
            datfilelogical(i)=true;
        else
            datfilelogical(i)=false;
        end
    else
        datfilelogical(i)=false;
    end
end
DirectoryStore=DirectoryStore(datfilelogical);
[~,sortind]=sort([DirectoryStore.datenum]);
DirectoryStore=DirectoryStore(sortind);
for i=1:length(DirectoryStore)
    outcell{i}=SnatchMPMSData(char([fileroot,char(DirectoryStore(i).name)]),MeasurementType);
    outcell{i}.Fileroot=char(fileroot);
    outcell{i}.Filename=char(DirectoryStore(i).name);
end


end
