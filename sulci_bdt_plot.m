function sulci_bdt_plot
close all
clc

Hemi = '/media/COSAS/Test/TOM/Data/brainvisa/ASPER_00001__101-20060510_Lhemi.mesh';
Surfh = load_mesh_lines(Hemi);
Surfh = rotatemat(Surfh,180,0,0);
Surfh.SurfData.vertices(:,1) = -1*Surfh.SurfData.vertices(:,1);
% Plot, Hemisphere
h = figure('numbertitle','off','Color', 'black','name',['IBASPM Surface Ploting...:  Plot ']);
strsurf=patch(Surfh.SurfData,'edgecolor','none','facecolor',[1 1 1], 'tag','patch','facelighting','gouraud');
set(strsurf,'FaceAlpha',0.8);
view([-114 10]); axis off; axis equal
camlight;

Sulci = '/media/COSAS/Test/TOM/Data/brainvisa/S.C._left.mesh';
Surfs = load_mesh_lines(Sulci); 
Surfs = rotatemat(Surfs,180,0,0);
Surfs.SurfData.vertices(:,1) = -1*Surfs.SurfData.vertices(:,1);
hold on;
strsurf=patch(Surfs.SurfData,'edgecolor','none','facecolor',[0 1 0], 'tag','patch','facelighting','gouraud');
set(strsurf,'FaceAlpha',0.5);




% Botline = '/media/COSAS/Test/TOM/Data/brainvisa/S.C._left_bottom.mesh';
% Surfb =  load_mesh_lines(Botline);
% Surfb = rotatemat(Surfb,180,0,0);
% Surfb.SurfData.vertices(:,1) = -1*Surfb.SurfData.vertices(:,1);
% hold on
% for k = 1:length(Surfb.SurfData.lines)
%     poss = Surfb.SurfData.lines{k};
%     line(Surfb.SurfData.vertices(poss,1),Surfb.SurfData.vertices(poss,2),Surfb.SurfData.vertices(poss,3),'Color',[1 0 0],'Linewidth',20);
% end

Depthline = '/media/COSAS/Test/TOM/Data/brainvisa/S.C._left_Depth.mesh';
Surfd = load_mesh_lines(Depthline)
Surfd = rotatemat(Surfd,180,0,0);
Surfd.SurfData.vertices(:,1) = -1*Surfd.SurfData.vertices(:,1);
hold on
for k = 1:length(Surfd.SurfData.lines)
    poss = Surfd.SurfData.lines{k};
    line(Surfd.SurfData.vertices(poss,1),Surfd.SurfData.vertices(poss,2),Surfd.SurfData.vertices(poss,3),'Color',[0 0 1],'Linewidth',5);
end


Topline = '/media/COSAS/Test/TOM/Data/brainvisa/S.C._left_top.mesh';
Surft = load_mesh_lines(Topline)
Surft = rotatemat(Surft,180,0,0);
Surft.SurfData.vertices(:,1) = -1*Surft.SurfData.vertices(:,1);
hold on
for k = 1:length(Surft.SurfData.lines)
    poss = Surft.SurfData.lines{k};
    line(Surft.SurfData.vertices(poss,1),Surft.SurfData.vertices(poss,2),Surft.SurfData.vertices(poss,3),'Color',[1 0 0],'Linewidth',20);
end


return

function Surf = rotatemat(Surf,x,y,z)
x = x*pi/180;
y = y*pi/180;
z = z*pi/180;
Rx = [1 0 0 0; 0 cos(x) -sin(x) 0; 0 sin(x) cos(x) 0; 0 0 0 1];
Ry = [cos(y) 0 sin(y) 0; 0 1 0 0; -sin(y) 0 cos(y) 0; 0 0 0 1];
Rz = [cos(z) -sin(z) 0 0; sin(z) cos(z) 0 0; 0 0 1 0; 0 0 0 1];
Mat = Rx*Ry*Rz;
t = [Surf.SurfData.vertices ones(size(Surf.SurfData.vertices,1),1)]*Mat';
Surf.SurfData.vertices = t(:,1:3);


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
        vert = reshape(vert,[3,Npoints])';i
        
        % Reading Normals
        Polt = fscanf(fid,'%d\n',1);Polt = fscanf(fid,'%d\n',1);Nfaces = fscanf(fid,'%d\n',1);
        
        % Reading Faces
        if Pol == 3
            faces = fscanf(fid,'(%d ,%d ,%d) ',Pol*Nfaces);
            
            Surf(t).SurfData.faces = faces
        elseif Pol ==2
            faces = fscanf(fid,'(%d ,%d) ',Pol*Nfaces); 
            faces = reshape(faces,[Pol Nfaces])'+1; 
            if faces(2,1)-faces(1,1)==1
                Surf(t).SurfData.lines{t} = unique(faces);
            else
                for k =1:size(faces,1)
                    Surf(t).SurfData.lines{k} = faces(k,:)';
                end
            end
        end
        Surf(t).Name = nm;
        Surf(t).SurfData.vertices = vert;clear faces vert;

      end
end
fclose(fid);