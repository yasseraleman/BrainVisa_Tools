function Multi_Unify_Tmtktri(TmtkmtriDir,ArgsDir);
%
% Syntax :
% Multi_Unify_Tmtktri(TmtkmtriDir,ArgsDir);
%
% Script file to unify giry surfaces in Brainvisa Tmtktri file according to 
% their name in the arg file.
%
% Input Parameters:
%   TmtkmtriDir       :  Directory to find any Brainvisa Tmtktri file
%   ArgsDir           :  Directory to find any Brainvisa arg file
%
% Output Parameters:
%
%
% Related references:
%
%
% See also: Unify_Tmtktri
% 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% April 9th 2012
% Version $1.0

[Tmtkmtrifiles] = sel_files(TmtkmtriDir,'aims_Tmtktri.mesh');
[ArgFilesfiles] = sel_files(ArgsDir,'*Rgyri_default_session_auto.arg');
Outargs = '';
OutTmtkmtris  = '';
Ns = size(ArgFilesfiles,1);
for i = 1:Ns
    disp(['Computing Subject ==> ' num2str(i) ' of ' num2str(Ns)]);
    TmtktriFile = deblank(Tmtkmtrifiles(i,:));
    ArgFile = deblank(ArgFilesfiles(i,:));
    [NewTmtktriFile,NewArgFile] = Unify_Tmtktri(TmtktriFile,ArgFile);
    OutTmtkmtris = strvcat(OutTmtkmtris,NewTmtktriFile);
    Outargs = strvcat(Outargs,ArgFile);
end
return


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

function [NewTmtktriFile,NewArgFile] = Unify_Tmtktri(TmtktriFile,ArgFile, NewTmtktriFile ,NewArgFile);
%
% Syntax :
% [NewTmtktriFile,NewArgFile] = Unify_Tmtktri(TmtktriFile,ArgFile, NewTmtktriFile ,NewArgFile);
%
% Script file to unify giry surfaces in Brainvisa Tmtktri file.
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
[pths,nms,exts] = fileparts(TmtktriFile);
[ptha,nma,exta] = fileparts(ArgFile);
if nargin <4
    NewArgFile = [ptha filesep nma '_mod.arg'];
end
if nargin <3
    NewTmtktriFile = [pths filesep nms '_mod.mesh'];
    NewArgFile = [ptha filesep nma '_mod.arg'];
end
[ptho,nmo,exto] = fileparts(NewTmtktriFile);
% ==============================================================%
% ============= Reading Arg text file  ==========================%
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
% ======================= Loading Mesh file  =============================%
Surf = load_mesh_lines(TmtktriFile);
%bindex = find(sum(isspace(lines'))==size(lines,2)); % Index of blank lines
% ========================================================================%
names = lines(ismember(lines(:,1:4),'name','rows'),:);
[ordnames,c,order] = unique(names(:,5:end),'rows');
tkmlabel = lines(ismember(lines(:,1:13),'Tmtktri_label','rows'),:);
ids = str2num(tkmlabel(:,14:end));
roilabel = lines(ismember(lines(:,1:9),'roi_label','rows'),:);

surfarea = lines(ismember(lines(:,1:12),'surface_area','rows'),:);
area = str2num(surfarea(:,15:end));
%============== Creating new structure and joining Surfaces ==============%
cont = 0;
for i = 1:size(ordnames,1)
    st(i).name = ordnames(i,:);
    ind = find(isspace(st(i).name));
    st(i).name(ind) = [];
    st(i).stcodes = find(order==i);
    if ~isempty(strfind(ordnames(i,:),'backgroundleft'))|~isempty(strfind(ordnames(i,:),'backgroundright'))
        st(i).roilabel = 0;
    else
        cont = cont+1;
        st(i).roilabel = cont;
    end
    areat = 0;vert = [0 0 0];faces = [0 0 0]; normals = [0 0 0];
    for j = 1:length(st(i).stcodes)
        faces = [faces;Surf(st(i).stcodes(j)).SurfData.faces+size(vert,1)];
        vert = [vert;Surf(st(i).stcodes(j)).SurfData.vertices];
        normals = [normals;Surf(st(i).stcodes(j)).SurfData.VertexNormals];
        areat = areat+area(st(i).stcodes(j));
    end
    faces = faces-1;faces(1,:) = [];vert(1,:) = [];normals(1,:) = [];
    Surft(i).Name = st(i).name; 
    Surft(i).SurfData.vertices = vert;
    Surft(i).SurfData.faces = faces;
    Surft(i).SurfData.VertexNormals = normals;
    Surft(i).Area = areat;
    fprintf(fiw,'%s\n',['*BEGIN NODE roi ' num2str(i)]);
    fprintf(fiw,'%s\n',['name          ' st(i).name]);
    fprintf(fiw,'%s\n',['Tmtktri_label ' num2str(i)]);
    fprintf(fiw,'%s\n',['roi_label     ' num2str(st(i).roilabel)]);
    fprintf(fiw,'%s\n',['surface_area  ' num2str(Surft(i).Area)]);
    fprintf(fiw,'%s\n','*END');
    fprintf(fiw,'\n');
end
fprintf(fiw,'%s','*END');
%=========================================================================%

Out_mesh = save_mesh(Surft,NewTmtktriFile);
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
