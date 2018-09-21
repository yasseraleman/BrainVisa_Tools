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