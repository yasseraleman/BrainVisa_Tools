function Outmeshes = Multi_Surf_plus_tal(Inputdir);
[Meshfiles] = sel_files(Inputdir,'*lh.white.mesh');
[Taltransfiles] = sel_files(Inputdir,'*talairach.lta');
Ns = size(Meshfiles,1);
Outmeshes = '';
for i = 1:Ns
    disp(['Processing mesh ' num2str(i) ' of ' num2str(Ns) ]);
    Meshfile = deblank(Meshfiles(i,:));
    Taltransf = deblank(Taltransfiles(i,:));
    [pth, nm, ext] = fileparts(Meshfile);
    Outmesh = [pth nm '+tal.mesh'];
    Outmesh = Surf_plus_tal(Meshfile,Taltransf,Outmesh);
    Outmeshes = strvcat(Outmeshes,Outmesh);
end
return

function Outmesh = Surf_plus_tal(Meshfile,Taltransf,Outmesh);
%
% Syntax :
% Outmesh = Surf+tal(Meshfile,Taltransf,Outmesh);
%
% Reads mesh and Texture files saved in Brainvisa Format and creates an output
% tmtkmtri file. 
%
% Input Parameters:
%   Meshfile    : Individual Mesh file
%  Taltransf    : Talairach center of coordinates in Freesurfer format
%  Outmesh      : Output tmtktri filename
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

%================== Checking input parameters ============================%

if nargin <3
    [pth, nm, txt] = fileparts(deblank(Meshfile));
    Outmesh = [pth filesep nm '+tal.mesh'];
end
%=========================================================================%
%==================== Reading Brainvisa surface ==========================%
Surf = load_mesh(Meshfile);

%==================== Reading Talairach center ===========================%
%if nargin ==6
cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
Surf.SurfData.vertices =Surf.SurfData.vertices+repmat(cras,[size(Surf.SurfData.vertices,1) 1]);
%end
%=========================================================================%
% ============= Reading Arg text file  ==========================%

Outmesh = save_mesh(Surf,Outmesh);
return;


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