function OutArg = Rename_sulci(InArg,Cname,Cth,OutArg,NLabel);
%
% Syntax :
% OutArg = Replace_Sulc_name_Argfile(InArg,Sulclist,Nsulcname,Outdir);
%
% This function replace the sulcis names with surface area lower than an
% specified threshold. The label is replaced by unknown
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
