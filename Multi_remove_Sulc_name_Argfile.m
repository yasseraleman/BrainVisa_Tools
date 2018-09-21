function Multi_remove_Sulc_name_Argfile(InputDir, hemi ,Outdir);
%
% Syntax :
% Multi_remove_Sulc_name_Argfile(InputDir, hemi ,Outdir);
%
% Script file to replace the sulcis names in brainvisa Arg files.
%
% Input Parameters:
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
% See also: Remove_Sulc_name_Argfile
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
Args = sel_files(InputDir,[cad '*_FS2BV_default_session_auto.arg']);
% ==============================================================%
% ============= Reading Nomenclature text file  ==========================%

Ns = size(Args,1);
for i = 1:Ns
    
    InArg=deblank(Args(i,:));
    disp(['Processing:  ==> ' InArg]); 
    disp(['Case ' num2str(i) ' of ' num2str(Ns)]);
    disp(['                    ']);
    [pth, nm, ext] = fileparts(InArg);
    if nargin<4
        Outdir = pth;
    end
    OutArg = [Outdir filesep nm '_rem.arg'];
    OutArg = Remove_sulci(InArg,OutArg);
end
return

function OutArg = Remove_sulci(InArg,OutArg);
%
% Syntax :
% OutArg = Replace_Sulc_name_Argfile(InArg,Sulclist,Nsulcname,Outdir);
%
% This function deletes the sulcis with labels unknown.
%
% Input Parameters:
%   InArg             : Input Arg File.
%   Cname             : Characteristics Name.
%   Cth               : Characteristics Threshold.
%   OutArg            : Output Arg filename.
%
% Output Parameters:
%   OutArg            : New/Output Arg file with the sulcis labels changed
%
%
% Related references:
%
%
% See also: Multi_Replace_Sulc_name_Argfile
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

if nargin<2
    OutArg = [pth filesep nm '_rem.arg'];
end
%=========================================================================%
%==========================  Reading Arg File  ===========================%
fio = fopen(InArg);
fir = fopen(OutArg,'wt');
cont = 0;
conti =0;conto = 0;stlines = '';
inds = 0; contind = 0; st = 0; conts =0; cadant = '';
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    if ~isempty(strfind(lower(line),'*begin node'));
        conti = conti+1;
        cadant = strvcat(cadant,line);
    end
    if ~isempty(strfind(lower(line),'*begin uedge'));
        st = 1;
    end
    if st ==0
    if conti-conto ==1
        inds(1) = [];
        % =======  Changing Characteristic ===============================
        indl = find(ismember(stlines(:,1:6),['label '],'rows') == 1);
        maxind = max(find(isspace(deblank(stlines(indl,:)))));
        sname = deblank(stlines(indl,maxind+1:end));
        if ~strcmp(sname,'unknown')
            if conti>1
                fprintf(fir,'%s\n',deblank(cadant(conti-1,:)));
            end
            for i = 1:size(stlines,1)
                if sum(ismember(inds,i))==0
                    fprintf(fir,'%s\n',stlines(i,:));
                else
                    fprintf(fir,'\n');
                end
            end
        end
        
        % ================================================================
       
        conto = conti;
        contind = 0;
        inds = 0;        
        stlines = '';
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