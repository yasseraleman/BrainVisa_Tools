function Multi_Change_Arg_file
labels = strvcat('S.F.inter._right','S.F.sup._right','S.F.inf._right');

toucs = sel_files('/media/LA-PUBLIC/temp','R*edited.arg');
origs = sel_files('/media/LA-PUBLIC/temp','R*orig.arg');
%orig = '/media/BORRARYA/LASPER_00001__101-20060510_default_session_auto_original.arg'
%touc = '/media/BORRARYA/LASPER_00001__101-20060510_default_session_auto_tocado.arg'

for j = 1:size(toucs,1)
    disp(['==== Runing Subject ===== ' num2str(j) ' of ' num2str(size(toucs,1))]);
    touc = deblank(toucs(j,:));
    orig = deblank(origs(j,:));
    for i = 1:size(labels,1)
        disp(['Changing Label ======>  ' labels(i,:) ' >> : ' num2str(i) ' of ' num2str(size(labels,1))]);
        label = deblank(labels(i,:));
        Nname = change_Arg_file(orig,touc,label);
        touc = Nname;
    end
    
    delete([Nname(1:end-12) '.arg']);
    delete([Nname(1:end-8) '.arg']);
    %delete([Nname(1:end-16) '.arg']);
    movefile(Nname,[Nname(1:end-12) '.arg'])
    disp(['!!! Finished:  ' Nname]);
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
    if strfind(lower(line),'label');
       if strcmp(lower(line(1:5)),'label')
           contil = contil+1;
           lidfio(contil) = cont;
       stro = strvcat(stro,line);
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
            
            pos = find(ismember(stro,strt,'rows'));
            
            %indt = find((lidfio - nidfio(pos))==min((lidfio - nidfio(pos))));
            nline  = stro(pos+1,:);
        end
    else
        nline = line;
    end
    if isempty(nline)
        fprintf(fitwn,'\n');
    else
        fprintf(fitwn,'%s\n',nline);
    end
end
fclose all;
clear tline cont;

function OutFiles = sel_files(InputDir,filtro);
TempDir = genpath(InputDir);

if ~isunix
    Pthst = strread(TempDir,'%s','delimiter',';'); Pthst = char(Pthst);
else
    Pthst = strread(TempDir,'%s','delimiter',':'); Pthst = char(Pthst);
end
pthold = pwd;
OutFiles = '';
cont = 0;
for i = 1:size(Pthst,1)
    cd(deblank(Pthst(i,:)));
    a = dir(filtro);
    if ~isempty(a)
    [names{1:size(a,1),1}]=deal(a.name);
    files = [repmat([deblank(Pthst(i,:)) filesep],[size(names,1) 1]) char(names)];
    OutFiles = strvcat(OutFiles,files);
    end
end
cd(pthold);
return