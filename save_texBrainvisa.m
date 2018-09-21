function [Txtfile] = save_texBrainvisa(Text,Txtfile, Filetype);
%
% Syntax :
% [Txtfile] = save_texBrainvisa(Text,Txtfile,Filetype);
%
% Saves texture values to Brainvisa Format. They can be saved in
% Binary or ASCII format
%
% Input Parameters:
%   Text: Struct Variable containing 2 fields.
%           Datatype: Original brainvisa datatype. 
%           Values:   Texture values. Each surface vertex contains a
%           texture value saved in the same order.
%   Txtfile     : Texture File (Binary or ASCII) Name
%   Filetype    : Save in binary or ascii text format
%
% Output Parameters:
%   Txtfile: Texture Filename
%
% Related references:
%
% See also:  read_texBrainvisa
%__________________________________________________
% Authors:  Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% September 16th 2011
% Version $1.0
%=====================Checking Input Parameters===========================%

if nargin<3
    Filetype = 'binar';
    %[Txtfile,sts] = spm_select([1],'image','Selecting Texture File','',cd);
elseif nargin<2
    errormdlg('Please specify a valid Texture filename','Texture Save Error...');
end
%=====================Main Program =======================================%
dat = Text.Datatype;
if strcmp(lower(dat),'uint16')
    Datatype = 'S16';
elseif strcmp(lower(dat)','uint32')
    Datatype = 'U32';
elseif strcmp(lower(dat),'float32')
    Datatype = 'FLOAT';
end
if strcmp(Filetype,'binar');
    fid = fopen(Txtfile,'wb');
    fwrite(fid, 'binar', 'uchar');
    fwrite(fid, hex2dec('41424344'), 'uint32');
    fwrite(fid,length(Datatype),'uint32');
    fwrite(fid,Datatype,'uchar');
    fwrite(fid,size(Text.Values,2),'uint32');
    for t=0:size(Text.Values,2)-1
        fwrite(fid,t,'uint32');
        fwrite(fid,size(Text.Values,1),'uint32');
        fwrite(fid,Text.Values(:,t+1),dat);
    end
elseif strcmp(Filetype,'ascii')
    fid = fopen(Txtfile,'wt');
    fprintf(fid,'%s\n','ascii');
    fprintf(fid,'%s\n',Datatype);
    fprintf(fid,'%u\n',size(Text.Values,2));
    temp=[[0:size(Text.Values,2)-1];repmat(size(Text.Values,1),[1 size(Text.Values,2)]);Text.Values];temp = temp(:)';
    for i = 1:size(temp,2);
        if i~=size(temp,2)
            fprintf(fid,'%u ',temp(i));
        else
            fprintf(fid,'%u',temp(i));
        end
    end;
end
fclose(fid);
return
%=========================================================================%