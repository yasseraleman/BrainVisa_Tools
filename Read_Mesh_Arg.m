function [Surf] = Read_Mesh_Arg(TmtktriFile,ArgFile);
%
% Syntax :
% [Surf] = Read_Mesh+Arg(TmtktriFile,ArgFile);
%
% Script file to read sulcal mesh files in Brainvisa Tmtktri file. The
% ArgFile is employed to label the meshes
%
% Input Parameters:
%   TmtktriFile       :  Brainvisa Tmtktri file
%   ArgFile           :  ArgFile containing the surface relationships in
%                        Tmtktri
%   NewTmtktriFile    : Unified Brainvisa Tmtktri file.
%   NewArgFile        : ArgFile corresponding to the new Tmtktri.
%
% Output Parameters:
%
%   NewTmtktriFile    : Unified Brainvisa Tmtktri file.
%   NewArgFile        : ArgFile corresponding to the new Tmtktri.
%
% Related references:
%
%
% See also: Replace_Sulc_name_Argfile
% 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% April 9th 2012
% Version $1.0

% ============= Input parameters ================================%
% ==============================================================%
% ============= Reading Arg text file  ==========================%
fio = fopen(ArgFile,'rt');lines = '';cont = 0;
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    lines = strvcat(lines,line);
    
end
fclose(fio);
% ========================================================================%
% ======================= Loading Mesh file  =============================%
Surf = load_mesh_lines(TmtktriFile);
%bindex = find(sum(isspace(lines'))==size(lines,2)); % Index of blank lines
% ========================================================================%
names = lines(ismember(lines(:,1:6),'label ','rows'),:);
[ordnames,c,order] = unique(names(:,6:end),'rows');
tkmlabel = lines(ismember(lines(:,1:13),'Tmtktri_label','rows'),:);
tkmlabel = tkmlabel(2:end,:);
labels = tkmlabel(:,14:end);
labels = str2num(labels);
%============== Creating new structure and joining Surfaces ==============%
cont = 0;
for i = 1:size(ordnames,1)
    st(i).name = ordnames(i,:);
    ind = find(isspace(st(i).name));
    st(i).name(ind) = [];
    st(i).stcodes = find(order==i);
    vert = [0 0 0];faces = [0 0 0]; normals = [0 0 0];
    for j = 1:length(st(i).stcodes)
        faces = [faces;Surf(st(i).stcodes(j)).SurfData.faces+size(vert,1)];
        vert = [vert;Surf(st(i).stcodes(j)).SurfData.vertices];
        normals = [normals;Surf(st(i).stcodes(j)).SurfData.VertexNormals];
    end
    faces = faces-1;faces(1,:) = [];vert(1,:) = [];normals(1,:) = [];
    Surft(i).Name = st(i).name; 
    Surft(i).SurfData.vertices = vert;
    Surft(i).SurfData.faces = faces;
    Surft(i).SurfData.VertexNormals = normals;
end
%=========================================================================%
Surf = Surft;
return;


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
return;