function Outmeshes = Multi_Surf_plus_tal(Inputdir);
[Meshfiles] = sel_files(Inputdir,'*lh.white+tal.aims.mesh');
[Txtfiles] = sel_files(Inputdir,'*lh.aparc.annot.tex');
[ArgFiles] = sel_files(Inputdir,'*_Lgyri_default_session_auto.arg');
Ns = size(Meshfiles,1);
Outmeshes = '';
for i = 1:Ns
    disp(['Processing mesh ' num2str(i) ' of ' num2str(Ns) ]);
    Meshfile = deblank(Meshfiles(i,:));
    Txtfile = deblank(Txtfiles(i,:));
    ArgFile = deblank(ArgFiles(i,:));
    
    [pth, nm, ext] = fileparts(Meshfile);
    Outmesh = [pth filesep nm '_Tmtkmtri.mesh'];
    NewArgFile = [pth filesep nm '_Tmtkmtri.arg'];
    Outmesh = mesh2tmtkmtri(Meshfile,Txtfile,ArgFile,Outmesh,NewArgFile);
    Outmeshes = strvcat(Outmeshes,Outmesh);
end
return

function Outmesh = mesh2tmtkmtri(Meshfile,Txtfile,ArgFile,Outmesh,NewArgFile,Taltransf,Imfile);
%
% Syntax :
% Outmesh = mesh2tmtkmtri(Meshfile,Txtfile,Outmesh,Taltransf);
%
% Reads mesh and Texture files saved in Brainvisa Format and creates an output
% tmtkmtri file. 
%
% Input Parameters:
%   Meshfile    : Individual Mesh file
%   Txtfile     : Texture File (Binary or ASCII)
%   ArgFile     : Reference Arg file to create a new argfile asociated to
%                 tmtktri mesh file
%  NewArgFile   : New argfile asociated to tmtktri mesh file
%  Outmesh      : Output tmtktri filename
%  Taltransf    : Talairach center of coordinates in Freesurfer format
%  Imfile       : Image file in freesurfer space.
%  
% Output Parameters:
%   Outmesh: : Output tmtktri filename
%
% Related references:
%
% See also:  save_texBrainvisa
%__________________________________________________
% Authors:  Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% April 16th 2012
% Version $1.0

% Meshfile = '/media/COSAS/Test/Joost/aquivan/lh.white+tal.aims.mesh';
% %Taltransf = '/media/COSAS/Test/Joost/aquivan/talairach.lta';
% Txtfile = '/media/COSAS/Test/Joost/aquivan/lh.aparc.annot.tex';
% ArgFile = '/media/COSAS/Test/Joost/aquivan/ASPER_00001__101-20060510_FS2BV_Lgyri_default_session_auto.arg';
% NewArgFile = '/media/COSAS/Test/Joost/aquivan/lh.white_Tmtkmtri.arg';
% Outmesh = '/media/COSAS/Test/Joost/aquivan/lh.white_Tmtkmtri.mesh';
% Imfile =  '/media/COSAS/Test/Joost/aquivan/ASPER_00001__101-20060510_T1.nii'
%================== Checking input parameters ============================%
if nargin ==4
    [pth, nm, txt] = fileparts(deblank(Outmesh));
    NewArgFile = [pth filesep nm '.arg'];
end
if nargin <4
    [pth, nm, txt] = fileparts(deblank(Meshfile));
    Outmesh = [pth filesep nm '_Tmtkmtri.mesh'];
    NewArgFile = [pth filesep nm '_Tmtkmtri.arg'];
end
%=========================================================================%
%==================== Reading Brainvisa surface ==========================%
Surf = load_mesh(Meshfile);

%==================== Reading Talairach center ===========================%
if nargin >6
    cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
    Surf.SurfData.vertices =Surf.SurfData.vertices+repmat(cras,[size(Surf.SurfData.vertices,1) 1]);
    Surf = freesCS2brainvisaCS(Surf,Imfile,'f2b');
end
%=========================================================================%
% ============= Reading Arg text file  ==========================%
[ptho,nmo,exto] = fileparts(Outmesh);
fio = fopen(ArgFile,'rt');lines = '';cont = 0;
fiw = fopen(NewArgFile,'wt');
conts = 0
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    linet= line;
    if isempty(line), conts = conts+1;indsp(conts) = cont;line = ' '; end
    lines = strvcat(lines,line);
    if conts<2
        if strfind(linet,'roi.global.tri  roi ');
            linet = ['roi.global.tri  roi ' nmo '.mesh Tmtktri_label'];
        end
        fprintf(fiw,'%s\n',linet);
    end
end
fclose(fio);
fprintf(fiw,'\n');
% ========================================================================%

[Text,Format] = read_texBrainvisa(Txtfile);
sts = unique(Text.Values);
sts(sts ==-1) = [];
Ns = length(sts);
for i = 1:Ns
    disp(['Processing ROI ' num2str(i) '  of ' num2str(Ns)]);
    ind = find(Text.Values == sts(i));
    indf = ismember(Surf.SurfData.faces,ind);
    indfaces = find(sum(indf')>0);
    Surftemp = Surf;
    Surftemp.SurfData.faces = Surf.SurfData.faces(indfaces,:);
    [Surft] = Reorg_Surf(Surftemp);
    Surftot(i) = Surft;
    [At] = Area_Comp(Surft.SurfData);
    Surftot(i).Area = sum(At);
    fprintf(fiw,'%s\n',['*BEGIN NODE roi ' num2str(i)]);
    if i<=3
        fprintf(fiw,'%s\n',['name          ' 'gyrus_' num2str(i)]);
    else
        fprintf(fiw,'%s\n',['name          ' 'gyrus_' num2str(i+1)]);
    end
    fprintf(fiw,'%s\n',['Tmtktri_label ' num2str(i)]);
    fprintf(fiw,'%s\n',['roi_label     ' num2str(i)]);
    fprintf(fiw,'%s\n',['surface_area  ' num2str(Surftot(i).Area)]);
    fprintf(fiw,'%s\n','*END');
    fprintf(fiw,'\n');
end
fprintf(fiw,'%s','*END');
Outmesh = save_mesh(Surftot,Outmesh);
return;

function [Surft] = Reorg_Surf(Surf);
%This scripts reorganize Surfaces in case of deleted faces. 

indf = unique(Surf.SurfData.faces(:));
%%%%%%%%%%%%%%%%%%%%%%%Corriegiendo Superficies para quedarme solo
%%%%%%%%%%%%%%%%%%%%%%%con la parte que me interesa
Surft = Surf;
Surft.SurfData.vertices = zeros(size(indf,1),3);
Surft.SurfData.VertexNormals = zeros(size(indf,1),3);
Surft.SurfData.faces = 0*Surft.SurfData.faces;
for i =1:size(indf,1)
    Surft.SurfData.vertices(i,:) = Surf.SurfData.vertices(indf(i),:);
    if isfield(Surf.SurfData,'VertexNormals')
        Surft.SurfData.VertexNormals(i,:) = Surf.SurfData.VertexNormals(indf(i),:);
    end
    indn = find(Surf.SurfData.faces ==indf(i));Surft.SurfData.faces(indn) = i;
end
[Tri] = Vert_Neib(double(Surft.SurfData.faces),size(Surft.SurfData.vertices,1),size(Surft.SurfData.faces,1));
Temp = sum(Tri);
Tri(:,Temp==0) = [];
Surft.Tri = Tri;Surft.Tri = Tri;
return;

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

function Surf = load_mesh(Filename)
[pth,nm,ext] = fileparts(Filename);
Sfile = [pth filesep nm ext(1:5)];
fid = fopen(deblank(Sfile),'rb');
Type = char(fread(fid, 5, 'uchar'));  %- 'ascii' or 'binar'
if strcmp(Type','binar');
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
    Text = char(fread(fid,4,'char'));
    Pol = fread(fid,1,'uint32');
    Tsteps = fread(fid,1,'uint32');
    for t=1:Tsteps
        Tinst = fread(fid,1,'uint32');
        
        % Reading vertices
        Npoints = fread(fid,1,'uint32');
        vert = fread(fid,Npoints*3,'float32')';
        vert = reshape(vert,[3,Npoints])';xyzn=vert;
        Surf(t).SurfData.vertices =vert;
        % Reading normals
        T = fread(fid,1,'uint32');
        normals = fread(fid,Npoints*3,'float32');
        normals = reshape(normals,[3,Npoints])';
        
        % Reading Faces
        T = fread(fid,1,'uint32');
        Nfaces = fread(fid,1,'uint32');
        faces = fread(fid,Nfaces*Pol,'uint32');
        faces = reshape(faces,[Pol,Nfaces])'+1;
        
        Surf(t).SurfData.VertexNormals = normals;clear normals
        Surf(t).SurfData.faces = faces;clear faces vert;
        if size(Surf(t).SurfData.faces,2) ==3
            [Tri] = Vert_Neib(double(Surf(t).SurfData.faces),Npoints,Nfaces);
            Temp = sum(Tri);
            Tri(:,Temp==0) = [];
            Surf(t).Tri = Tri;
        end
        if Tsteps~=1
            Surf(t).Name = [nm '_' sprintf('%.3d',t)];
        else
            Surf(t).Name = nm;
        end
        Surf(t).Area = 0;
        Surf(t).Type = 'Mask';
    end
elseif strcmp(Type(1:5,1)','ascii');
    Text = fscanf(fid,'%s',1);
    Pol = fscanf(fid,'%d',1);
    Tsteps = fscanf(fid,'%d',1);
    for t=1:Tsteps
        ms = fscanf(fid,'\n%d',1);
        
        % Reading vertices
        Npoints = fscanf(fid,'\n%d\n',1);
        vert = fscanf(fid,'(%f ,%f ,%f) ',3*Npoints);
        vert = reshape(vert,[3,Npoints])';
        
        % Reading Normals
        if Pol==3
            T = fscanf(fid,'\n%d\n',1);
            normals = fscanf(fid,'(%f ,%f ,%f) ',3*Npoints);
            normals = reshape(normals,[3,Npoints])';
            T = fscanf(fid,'\n%d\n',1);
        end
        
        % Reading Faces
        Nfaces = fscanf(fid,'\n%d\n',1);
        faces = fscanf(fid,'(%d ,%d ,%d) ',Pol*Nfaces);
        faces = reshape(faces,[Pol,Nfaces])'+1;
        Surf(t).Name = nm;
        Surf(t).SurfData.vertices = vert;
        Surf(t).SurfData.VertexNormals = normals;clear normals
        Surf(t).SurfData.faces = faces;clear faces vert;
        if size(Surf(t).SurfData.faces,2) ==3
            [Tri] = Vert_Neib(double(Surf(t).SurfData.faces),Npoints,Nfaces);
            Temp = sum(Tri);
            Tri(:,Temp==0) = [];
            Surf(t).Tri = Tri;
        end
        Surf(t).Area = 0;
        Surf(t).Type = 'Mask';
    end
end
fclose(fid);
return

function Out_mesh = save_mesh(Surf,Out_mesh);
%
% Syntax :
% Out_mesh = save_mesh(Surf,Out_mesh);
%
% Script file to save brainvisa mesh files
%
% Input Parameters:
%   Surf              :  Surface variable
%   Out_mesh          :  Output mesh filename
%
% Output Parameters:
%   Out_mesh          :  Output mesh filename
%
%
% Related references:
%
%
% See also: load_mesh
% 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% April 9th 2012
% Version $1.0

%=================== Main Program ========================================%

fid = fopen(Out_mesh,'wb');
fwrite(fid, 'binar', 'uchar');
fwrite(fid, hex2dec('41424344'), 'uint32');
fwrite(fid, 4, 'uint32');
fwrite(fid, 'VOID', 'uchar');
fwrite(fid, size(Surf(1).SurfData.faces,2), 'uint32');
fwrite(fid, size(Surf,2), 'uint32');
for k=1:size(Surf,2)
    vert = Surf(k).SurfData.vertices;
    face = Surf(k).SurfData.faces;
    if isfield(Surf(k).SurfData,'VertexNormals')
        normals = Surf(k).SurfData.VertexNormals;
    elseif (~isfield(Surf(k).SurfData,'VertexNormals'))
        fv.vertices = vert;
        fv.faces = face;
        h = figure;h1 = patch(fv);  normals = get(h1,'VertexNormals');close(h);clear fv;
        norma = sqrt(sum((normals').^2));
        normals = normals./repmat(norma',[1 3]);
    end
    Nfaces = size(Surf(k).SurfData.faces,1);
    Npoints = size(Surf(k).SurfData.vertices,1);
    face = face';face = face(:)-1;
    normals = normals';normals = normals(:);
    vert = vert';vert = vert(:);
    if size(Surf,2) == 1
        fwrite(fid, k-1, 'uint32');
    else
        fwrite(fid, k, 'uint32');
    end
    fwrite(fid, size(vert,1)/3, 'uint32');
    fwrite(fid, vert', 'float32');
    fwrite(fid, size(normals,1)/3, 'uint32');
    fwrite(fid, normals', 'float32');
    fwrite(fid, 0, 'uint32');
    fwrite(fid, size(face,1)/3, 'uint32');
    fwrite(fid, face', 'uint32');
end
fclose(fid);
%========================================================================%
return