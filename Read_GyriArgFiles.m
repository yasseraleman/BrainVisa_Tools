function [Lines, StNames] = Read_GyriArgFiles(ArgFile);

fio = fopen(ArgFile,'rt');lines = '';cont = 0;
contmax = 0;
Lines = '';
StNames = '';

while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    if strfind(line,'name ')
        contmax = contmax + 1;
        temp = strread(line,'%s');
        StNames = strvcat(StNames,char(temp(2)));
        
    end
    Lines = strvcat(Lines,line);
end
fclose(fio);
return