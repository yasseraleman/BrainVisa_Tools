function OutArg = Replace_Sulc_name_Argfile(InArg, Sulclist, Nsulcname, OutArg);
%
% Syntax :
% OutArg = Replace_Sulc_name_Argfile(InArg,Sulclist,Nsulcname,Outdir);
%
% This function replace the sulcis names, contained in Sulclist, with a new
% sulc name in the Arg File (Brainvisa Sulci organization format)
%
% Input Parameters:
%   InArg             : Input Arg File
%   Sulclist          : List names for sulcis that will be renamed.
%   Nsulcname         : New Sulcis name
%   Outdir            : Output Directory.
%
% Output Parameters:
%   OutArg            : New/Output Arg file with the sulcis names changed
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
Ns = size(Sulclist,1);
Sult = Sulclist';Sult = Sult(:)'; inds = isspace(Sult);Sult(inds) = [];
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
    if strfind(lower(line),'label');
        if strcmp(lower(line(1:5)),'label')
            ind = find(isspace(line)==1);
            sname = deblank(line(ind(end)+1:end));
            if ~isempty(strfind(Sult,sname));
                nline = [line(1:ind(end)) deblank(Nsulcname)];
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
