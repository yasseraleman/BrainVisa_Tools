function Surf = load_vtk(Filename)
%Filename = '/media/Data/R_Thal1_Fthresh.vtk';
% Filename = '/media/BORRARYA/joost_DTI/JOOST/BUNDLES_LIN/ASPER_00001/ASPER_00001_arcurate_L.vtk';
[pth,nm,ext] = fileparts(deblank(Filename));
Sfile = [pth filesep nm ext];
% Surf.Name = nm;
% fid = fopen(Sfile, 'rt');
% line = fgetl(fid);line = fgetl(fid);ftype = deblank(fgetl(fid));
% line = fgetl(fid);
% [Info,typ] = strread(line,'%s%s','delimiter',' ');
% if ~strcmp(lower(typ),'polydata')
%     errordlg('Please select a correct surface format');
%     return;
% end
% if strcmp(lower(ftype),'ascii')
%     line = fgetl(fid);
%     % Reading Vertices
%     [txt,Npoints,typ] = strread(line,'%s%n%s','delimiter',' ');
%     Surf.SurfData.vertices = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',5),3,Npoints)';
%
%     % Reading Faces
%     [line,Nfaces] = textread(Sfile,'%s %u',1,'headerlines',5+Npoints);
%     pol = textread(Sfile,'%u',1,'headerlines',6+Npoints);pol = pol+1;
%     Surf.SurfData.faces = reshape(textread(Sfile,'%f',Nfaces*pol,'headerlines',6+Npoints),pol,Nfaces)';
%     Surf.SurfData.faces(:,1)=[]; Surf.SurfData.faces= Surf.SurfData.faces+1;
%     Surf.Is = textread(Sfile,'%f',Npoints,'headerlines',9+Npoints+Nfaces);
%     if ~isempty(Surf.Is)
%         Surf.SurfData.VertexNormals = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',10+2*Npoints+Nfaces),3,Npoints)';
%     else
%         Surf =rmfield(Surf,'Is');
%     end
%     [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
%     Temp = sum(Tri);
%     Tri(:,Temp==0) = [];
%     Surf.Tri = Tri; clear Tri;
%
% else
%     V = 0;
%     fclose(fid);
%     % open file (OBS! big endian format)
%     fid = fopen(Sfile,'r','b');
%
%     if( fid == -1 )
%         return
%     end
%
%     fgetl(fid); % # vtk DataFile Version x.x
%     fgetl(fid); % comments
%     fgetl(fid); % BINARY
%     fgetl(fid); % DATASET STRUCTURED_POINTS
%
%     s = fgetl(fid); % DIMENSIONS NX NY NZ
%     sz = sscanf(s, '%*s%d%d%d').'
%
%     fgetl(fid); % ORIGIN OX OY OZ
%     fgetl(fid); % SPACING SX SY SZ
%     fgetl(fid); % POINT_DATA NXNYNZ
%
%     s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
%     svstr = sscanf(s, '%s', 1)
%     dtstr = sscanf(s, '%*s%*s%s')
%
%     if( strcmp(svstr,'SCALARS') > 0 )
%         fgetl(fid); % the lookup table
%         if( strcmp(dtstr,'unsigned_char') > 0 )
%             % read data
%             V = fread(fid,prod(sz),'*uint8');
%             V = reshape(V,sz);/media/COSAS/scripts/BrainVisa_myTools/load_vtk.m
%         elseif( strcmp(dtstr,'char') > 0 )
%             % read data
%             V = fread(fid,prod(sz),'*int8');
%             V = reshape(V,sz);
%         elseif( strcmp(dtstr,'unsigned_short') > 0 )
%             % read data
%             V = fread(fid,prod(sz),'*uint16');
%             V = reshape(V,sz);
%         elseif( strcmp(dtstr,'short') > 0 )
%             % read data
%             V = fread(fid,prod(sz),'*int16');
%             V = reshape(V,sz);
%         elseif( strcmp(dtstr,'unsigned_int') > 0 )
%             % read data
%             V = fread(fid,prod(sz),'*uint32');
%             V = reshape(V,sz);
%         elseif( strcmp(dtstr,'int') > 0 )
%             % read data
%             V = fread(fid,prod(sz),'*int32');
%             V = reshape(V,sz);
%         elseif( strcmp(dtstr,'float') > 0 )
%             % read data
%             V = fread(fid,prod(sz),'*single');
%             V = reshape(V,sz);
%         elseif( strcmp(dtstr,'double') > 0 )
%             % read data
%             V = fread(fid,prod(sz),'*double');
%             V = reshape(V,sz);
%         end
%
%     elseif( strcmp(svstr,'VECTORS') > 0 )
%         if( strcmp(dtstr,'float') > 0 )
%             % read data
%             V = fread(fid,3*prod(sz),'*single');
%             V = reshape(V,[3 sz]);
%             V = permute(V,[2 3 4 1]);
%         end
%     end
%
%     fclose(fid);
%
% end

fid = fopen(Sfile);
Surf.Name = nm;
cont = 0;
cadt = '';pointpos = 0;
while 1
    cont = cont + 1;
    line = fgetl(fid);
    
    if ~ischar(line),   break,   end
    %cadt = strvcat(cadt,line);
    ind = strfind(lower(line),'dataset');
    if ~isempty(ind)
        types = strread(line,'%s','delimiter',' ');
        typ = deblank(char(types{2}));
    end
    ind = strfind(lower(line),'points');
    if ~isempty(ind)
        pointpos = cont;
        types = strread(line,'%s','delimiter',' ');
        Npoints = str2num(deblank(char(types{2})));
        typ = deblank(char(types{3}));
        if strcmp(typ,'float')
            charac = '%f';
        end
    end
    ind = strfind(lower(line),'lines');
    if ~isempty(ind)
        linepos = cont;
        types = strread(line,'%s','delimiter',' ');
        Nlines = str2num(deblank(char(types{2})));
    end
    if exist('linepos','var')
        if (cont >=linepos+1)&(cont <=linepos+Nlines+1)
            cadt = strvcat(cadt,line);
        end
    end
    ind = strfind(lower(line),'polygons');
    if ~isempty(ind)
        polypos = cont;
        types = strread(line,'%s','delimiter',' ');
        Npoly = str2num(deblank(char(types{2})));
    end
    ind = strfind(lower(line),'point_data');
    if ~isempty(ind)
        charcpos = cont+2;
    end
    ind = strfind(lower(line),'vectors');
    if ~isempty(ind)
        vectorspos = cont;
    end
end
fclose(fid);
clear tline cont;
%% Reading Points
if exist('pointpos','var')
    Surf.SurfData.vertices = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',pointpos),3,Npoints)';
end
%% Reading Faces
if exist('polypos','var')
    pol = textread(Sfile,'%u',1,'headerlines',polypos);pol = pol+1;
    Surf.SurfData.faces = reshape(textread(Sfile,'%u',Npoly*pol,'headerlines',polypos),pol,Npoly)';
    Surf.SurfData.faces(:,1)=[]; Surf.SurfData.faces= Surf.SurfData.faces+1;
end
%% Reading Lines
if exist('linepos','var')
    cont = 0;
    for i = 1:Nlines
        pol = strread(cadt(i,:),'%u','delimiter',' ');
        Surf.SurfData.lines{i} = pol(2:end)+1;
    end
end
%% Reading Characteristics
if exist('charcpos','var')
    Surf.Is = textread(Sfile,'%f',Npoints,'headerlines',charcpos);
end
%% Reading Normals
if exist('vectorspos','var')
    Surf.SurfData.VertexNormals = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',vectorspos),3,Npoints)';
end
return





