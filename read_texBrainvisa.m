function [Text,Format] = read_texBrainvisa(Txtfile)
%
% Syntax :
% Text = read_texBrainvisa(Txtfile)
%
% Reads texture files saved in Brainvisa Format. They can be saved in
% Binary or ASCII format
%
% Input Parameters:
%   Txtfile     : Texture File (Binary or ASCII)
%
% Output Parameters:
%   Text: Struct Variable containing 2 fields.
%           Datatype: Original brainvisa datatype. 
%           Values:   Texture values. Each surface vertex contains a
%           texture value saved in the same order.
%
% Related references:
%
% See also:  save_texBrainvisa
%__________________________________________________
% Authors:  Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% September 16th 2011
% Version $1.0
%=====================Checking Input Parameters===========================%

if nargin==0
    [Txtfile,sts] = spm_select([1],'image','Selecting Texture File','',cd);
end
%=====================Main Program =======================================%
fid = fopen(Txtfile,'r');
Type = char(fread(fid, 5, 'uchar'));  %- 'ascii' or 'binar'
if strcmp(Type','binar');
    Format = 'binar';
    [byte_swapping, COUNT]     = fread(fid, 1, 'uint32'); %- 'ABCD' or 'DCBA'
    ff = strcmp(dec2hex(byte_swapping),'41424344');
    if ~ff
        [fn, pm, mf] = fopen(1); %- machine format
        fclose(fid);
        if strmatch(mf,'ieee-le');
            fid = fopen(deblank(Sfile),'r','ieee-be');
        else
            fid = fopen(deblank(Sfile),'r','ieee-le');
        end
        [file_format, COUNT]   = fread(fid, 5, 'uchar');
        [byte_swapping, COUNT] = fread(fid, 1, 'uint32');
    end
    Val = fread(fid,1,'uint32');
    dat = char(fread(fid,Val,'char'));
    if strcmp(lower(dat)','s16')
        Datatype = 'uint16';
    elseif strcmp(lower(dat)','u32')
        Datatype = 'uint32';
    elseif strcmp(lower(dat)','float')
        Datatype = 'float32';
    end
    Pol = fread(fid,1,'uint32');
    Text.Datatype = Datatype;
    Values = 0;
    for t=1:Pol
        Tsteps = fread(fid,1,'uint32');
        Npoints = fread(fid,1,'uint32');
        vert = fread(fid,Npoints,Datatype)';
        Values = [Values; vert'];
    end
    Text.Values = reshape(Values(2:end),[Npoints Pol] );
elseif strcmp(Type','ascii')
   Format = 'ascii';
   temp = fgetl(fid);
   dat = fgetl(fid);
   if strcmp(lower(dat),'s16')
        Datatype = 'uint16';
        dt = '%u';
    elseif strcmp(lower(dat),'u32')
        Datatype = 'uint32';
        dt = '%u';
    elseif strcmp(lower(dat),'float')
        Datatype = 'float32';
        dt = '%f';
   end
   Text.Datatype = Datatype;
   Pol = str2num(fgetl(fid));
   Values = textread(Txtfile,dt,'delimiter',' ','headerlines',3);
   Npoints = Values(2);
   Text.Values = reshape(Values,[Npoints+2 Pol] );Text.Values(1:2,:) = [];
else
    fclose(fid);
    Values = textread(Txtfile);
    Nstep = (size(Values,1)-2)/Values(2);
    Text.Datatype = 'float32';
    Text.Values = reshape(Values,[Values(2)+2 Nstep]);Text.Values = Text.Values(3:end,:);
end
fclose(fid);
return
%=========================================================================%