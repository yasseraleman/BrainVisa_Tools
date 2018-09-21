function Multi_Replace_Gyri_name_Argfile(ChangeTxtFile,InputDir, hemi ,Outdir);
%
% Syntax :
% Multi_Replace_Gyri_name_Argfile(ChangeTxtFile,InputDir, hemi ,Outdir);
%
% Script file to replace the gyri names in brainvisa Arg files.
%
% Input Parameters:
%   ChangeTxtFile     : Modified Nomenclature text file
%   InputDir          : Input directory where Arg files are stored (It can be
%                       the brainvisa root directory).
%   hemi              : Hemisphere
%   Outdir            : Output Directory.
%
% Output Parameters:
%
%
%
% Related references:
%
%
% See also: Replace_Gyri_name_Argfile
% 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 20th 2012
% Version $1.0

% ============= Input parameters ================================%
if strcmp(hemi,'left')
    cad = 'L';
elseif strcmp(hemi,'right')
    cad = 'R';
end
Args = sel_files(InputDir,[cad '*_FS2BV_Lgyri_default_session_auto.arg']);
Cname = strvcat('surface_area'); 
% If you want to put new characteristics you have to add it this way
% Cname = strvcat('surface_area';'new_charact1';'new_charact2';...etc); 
Cth = [5000000000];
% Cth = [50 22 9 ...etc];
NLabel =strvcat('unknown');
% NLabel = strvcat('unknown';'new_name1';'new_name2';...etc); 
% ==============================================================%
% ============= Reading Nomenclature text file  ==========================%
fio = fopen(ChangeTxtFile,'rt');lines = '';cont = 0;
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    lines = strvcat(lines,line);
end
fclose(fio);
% ========================================================================%
for j= 1:size(lines,1)
    sts = strread(deblank(lines(j,:)),'%s','delimiter',' ');
    Gyrilist = '';
    for k=1:size(sts,1)
        if k==1
            Ngyriname = char(sts{k});
            Ngyriname = [Ngyriname(1:end-1) hemi];
        else
            Gyrilist = strvcat(Gyrilist,[char(sts{k}) hemi]);
        end
    end
    stsl{j} = strvcat(Ngyriname,Gyrilist);
end
Ns = size(Args,1);
for i = 1:Ns
    
    InArg=deblank(Args(i,:));
    disp(['Renaming:  ==> ' InArg]); 
    disp(['Case ' num2str(i) ' of ' num2str(Ns)]);
    disp(['                    ']);
    
    [pth, nm, ext] = fileparts(InArg);
    if nargin<4
        Outdir = pth;
    end
    OutArg = [Outdir filesep nm '_mod.arg'];
    TempArg = [Outdir filesep nm '_temp.arg'];
    for j = 1:size(stsl,2)
        Temp = char(stsl{j});
        Gyrilist = Temp(2:end,:);
        Ngyriname = deblank(Temp(1,:));
        OutArg = Replace_Gyri_name_Argfile(InArg,Gyrilist,Ngyriname,OutArg);
        copyfile(OutArg,TempArg);
        InArg = TempArg;
    end
%     for k = 1:size(Cname,1)
%         if size(NLabel,1)==1
%             NLabelname = NLabel;
%         else 
%             NLabelname = deblank(NLabel(k,:));
%         end
%         OutArg = Rename_gyri(TempArg,deblank(Cname(k,:)),Cth(k),OutArg,NLabelname);
%         copyfile(OutArg,TempArg);
%     end
    delete(TempArg);
end
return

function OutArg = Replace_Gyri_name_Argfile(InArg,Gyrilist,Ngyriname,OutArg);
%
% Syntax :
% OutArg = Replace_Gyri_name_Argfile(InArg,Gyrilist,Ngyriname,Outdir);
%
% This function replace the gyri names, contained in Gyrilist, with a new
% gyri name in the Arg File (Brainvisa gyri organization format)
%
% Input Parameters:
%   InArg             : Input Arg File
%   Gyrilist          : List names for gyriis that will be renamed.
%   Ngyriname         : New Gyriis name
%   Outdir            : Output Directory.
%
% Output Parameters:
%   OutArg            : New/Output Arg file with the gyriis names changed
%
%
% Related references:
%
%
% See also: Multi_Replace_Gyri_name_Argfile
% 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 20th 2012
% Version $1.0


%%
[pth, nm, ext] = fileparts(InArg);
Ns = size(Gyrilist,1);
Sult = Gyrilist';Sult = Sult(:)'; inds = isspace(Sult);Sult(inds) = [];
%=====================Checking Input Parameters===========================%
% if nargin<4
%     Outdir = pth;
% end
%=========================================================================%
%==========================  Reading Arg File  ===========================%
fio = fopen(InArg);
cont = 0;
fiotxt = '';
conti =0;stro = '';
contil = 0; strol = '';
fir = fopen(OutArg,'wt');
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    if strfind(lower(line),'name');
        if strcmp(lower(line(1:4)),'name')
            ind = find(isspace(line)==1);
            sname = deblank(line(ind(end)+1:end));
            if ~isempty(strfind(Sult,sname));
                nline = [line(1:ind(end)) deblank(Ngyriname)];
            else
                nline = line;
            end
            
        else
            nline= line;
        end
    else
        nline = line;
    end
    
    fprintf(fir,'%s\n',nline);
end
fclose(fio);fclose(fir);
return

function OutArg = Rename_gyri(InArg,Cname,Cth,OutArg,NLabel);
%
% Syntax :
% OutArg = Replace_Gyri_name_Argfile(InArg,Gyrilist,Ngyriname,Outdir);
%
% This function replace the gyri names with surface area lower than a
% specified threshold. The label is replaced by unknown
%
% Input Parameters:
%   InArg             : Input Arg File.
%   Cname             : Characteristics Name.
%   Cth               : Characteristics Threshold.
%   OutArg            : Output Arg filename.
%
% Output Parameters:
%   OutArg            : New/Output Arg file with the gyriis labels changed
%
%
% Related references:
%
%
% See also: Multi_Replace_Gyri_name_Argfile
% 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 20th 2012
% Version $1.0

%%
[pth, nm, ext] = fileparts(InArg);
%=====================Checking Input Parameters===========================%

if nargin<5
    NLabel = 'unknown';
end
if nargin<4
    OutArg = [pth filesep nm '_mod.arg'];
    NLabel = 'unknown';
end
%=========================================================================%
%==========================  Reading Arg File  ===========================%
fio = fopen(InArg);
fir = fopen(OutArg,'wt');
l = length(Cname);
cont = 0;
conti =0;conto = 0;stlines = '';
inds = 0; contind = 0; st = 0; conts =0;
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    if ~isempty(strfind(lower(line),'*begin node'));
        conti = conti+1;
    end
    if ~isempty(strfind(lower(line),'*begin uedge'));
        st = 1;
    end
    if st ==0
    if conti-conto ==1
        inds(1) = [];
        % =======  Changing Characteristic ===============================
        ind = find(ismember(stlines(:,1:l+1),[Cname ' '],'rows') == 1);
        indl = find(ismember(stlines(:,1:6),['label '],'rows') == 1);
        if str2num(deblank(stlines(ind,l+1:end)))<Cth;
            maxind = max(find(isspace(deblank(stlines(indl,:)))));
            Newline = [stlines(indl,1:maxind) NLabel];
            if size(Newline,2)>size(stlines,2)
                stlines = [stlines repmat(' ', [size(stlines,1) (size(Newline,2)-size(stlines,2))])];
            end
            stlines(indl,1:size(Newline,2)) = Newline;
            stlines(indl,size(Newline,2)+1:end) = repmat(' ',[1  length(size(Newline,2)+1:size(stlines,2))]);
        end
        
        % ================================================================
        for i = 1:size(stlines,1)
            if sum(ismember(inds,i))==0
                fprintf(fir,'%s\n',stlines(i,:));
            else
                fprintf(fir,'\n');
            end
        end
        conto = conti;
        contind = 0;
        inds = 0;        
        stlines = '';
        fprintf(fir,'%s\n',line);
    else
        contind = contind+1;
        if isempty(line)
            inds = [inds contind];
            stlines = strvcat(stlines,' ');
        else
            stlines = strvcat(stlines,line);
        end
        
    end
    else
        conts = conts+1;
        if conts == 1
        for i = 1:size(stlines,1)
            if sum(ismember(inds,i))==0
                fprintf(fir,'%s\n',stlines(i,:));
            else
                fprintf(fir,'\n');
            end
        end
        fprintf(fir,'%s\n',line);
        else
            fprintf(fir,'%s\n',line);
        end
    end
end
fclose(fio);fclose(fir);
clear tline cont;
return;

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