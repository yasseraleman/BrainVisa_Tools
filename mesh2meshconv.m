function [OutMeshs,Gi] = mesh2meshconv(InputDir);
%
% Syntax :
%  [OutMeshs] = mesh2meshconv(MeshFile,OutMesh);
%
% This script file converts mesh files to point files to obtain convex hull 
% representations of the cortex
%
% Input Parameters:
%   MeshFile          :  Brainvisa Mesh file
%   OutMesh           :  Output ascii filename
%
% Output Parameters:
%
%  OutMesh           :  Output Ascii filename
%
% Related references:
%
%
% See also: 
% 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% April 9th 2012
% Version $1.0

% ============= Input parameters ================================%

% ========================================================================%
MeshFiles = sel_files(InputDir,'*.pial.mesh');
%MeshFiles = '/media/Data/ASPER_00001__101-20060510.lh.r.aims.pial.mesh';
opts.method = 'matlab';
% ======================= Loading Mesh file  =============================%
Ns = size(MeshFiles,1);
OutMeshs = '';
for i = 1:Ns   
    MeshFile = deblank(MeshFiles(i,:));
    [pths,nms,exts] = fileparts(MeshFile);
    Surf = load_mesh_lines(MeshFile);
    
    % ======= Using Brainvisa Method =====================================%
    switch opts.method
        case 'brainvisa'
            OutMesh = [pths filesep nms '_ascii.txt'];
            fid = fopen(OutMesh,'wt');
            fprintf(fid,'%s\n',num2str(size(Surf.SurfData.vertices,1)));
            fprintf(fid,'%f %f %f\n',Surf.SurfData.vertices');
            fclose all;
            OutMeshs = strvcat(OutMeshs,Out_mesh);
            Gi = 0;
        case 'matlab'
            a = convhulln(Surf.SurfData.vertices);
            [A] = Area_Comp(Surf.SurfData);
            Surft = Surf;
            Surft.SurfData.faces = a;
            [Surft] = Reorg_Surf(Surft);
            [At] = Area_Comp(Surft.SurfData);
            Gi = sum(A(:))/sum(At(:));
            Out_mesh = save_mesh(Surft,[pths filesep nms '_convhull.mesh']);
            OutMeshs = strvcat(OutMeshs,Out_mesh);
    end
    %=====================================================================%
end


%=========================================================================%
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

function Surf = load_mesh_lines(Filename)
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
    fwrite(fid, k, 'uint32');
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
