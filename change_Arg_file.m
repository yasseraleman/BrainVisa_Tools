function Multi_Change_Arg_file
labels = strvcat('','S.F.sup._left')
orig = '/media/BORRARYA/LASPER_00001__101-20060510_default_session_auto_original.arg'
touc = '/media/BORRARYA/LASPER_00001__101-20060510_default_session_auto_tocado.arg'
for i = 1:size(labels,2)
    label = deblank(labels(i,:));
    Nname = change_Arg_file(orig,touc,label)
    touc = Nname;
end
return


function Nname = change_Arg_file(orig,touc,label)

[pth,name,ext] = fileparts(touc);
Nname = [pth filesep name '_mod.arg'];
fio = fopen(orig);
cont = 0;
fiotxt = '';
conti =0;stro = '';
contil = 0; strol = '';
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    %fiotxt = strvcat(fiotxt,line);
    if ~isempty(strfind(lower(line),'*begin node'));
       conti = conti+1;
       nidfio(conti) = cont;
       stro = strvcat(stro,line);
    end
    if cont == 74
        a =1
    end
    if strfind(lower(line),'label');
       if strcmp(lower(line(1:5)),'label')
           contil = contil+1;
           lidfio(contil) = cont;
       strol = strvcat(strol,line);
       end
    end
    
end
fclose(fio);
clear tline cont;




fit = fopen(touc);
fitwn = fopen(Nname,'wt');
cont = 0;
fiotxt = '';
conti =0;strt = '';
contil = 0; strtl = '';
while 1
    cont = cont + 1;
    line = fgetl(fit);
    if ~ischar(line),   break,   end
    
    %%
    if ~isempty(strfind(lower(line),'*begin node'));
       %conti = conti+1;
       %nidfit(conti) = cont;
       %strt = strvcat(strt,line);
       strt = line;
    end
%     if strfind(lower(line),'label');
%        if strcmp(lower(line(1:5)),'label')
%            contil = contil+1;
%            lidfit(contil) = cont;
%        strtl = strvcat(strtl,line);
%        end
%     end
    %%
    
    if ~isempty(strfind(lower(line),lower(label)))
        if strcmp(lower(line(1:5)),'label')
            indtemp = find(sum(ismember(stro(:,1:length(strt)),strt)') == length(strt));
            for i = 1:length(indtemp)
                if strcmp(lower(deblank(stro(indtemp(i),:))),lower(strt))
                    ind = indtemp(i);
                end
            end
            pos = nidfio(ind);
            indt = find((lidfio - pos)==min((lidfio - pos)));
            nline  = deblank(strol(indt,:));
        end
    else
        nline = line;
    end
    if isempty(nline)
        fprintf(fitwn,'\r');
    else
        fprintf(fitwn,'%s\r',nline);
    end
end
fclose all;
clear tline cont;stro