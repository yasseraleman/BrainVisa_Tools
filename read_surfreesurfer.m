function Surf= read_surfreesurfer(Sfile);
[pth,nm,ext] = fileparts(deblank(Sfile));
form =ext(2:4);

Tmn =  16777214 ;
Qmn =  16777215 ;
fid = fopen(deblank(Sfile),'rb','b');
if fid >0
    MNumber = bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar');
    if(MNumber == Qmn)
        Surf.Imp = 'frb';
        Npoints = bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar');
        Nfaces = bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar');
        Surf.SurfData.vertices = fread(fid, Npoints*3, 'int16') ./ 100 ;
        for t=1:Npoints
            for k=1:4
                Surf.SurfData.faces(t,k) =  bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar'); ;
            end
        end
        Surf.SurfData.faces = Surf.SurfData.faces+1;
    elseif (MNumber == Tmn)
        Surf.Imp = 'frb';
        tline = fgets(fid);
        tline = fgets(fid);
        Np = fread(fid,1,'int32');
        Nf = fread(fid,1,'int32');
        vertices = fread(fid,Np*3,'float32');
        Surf.Name = [nm ext];
        Surf.SurfData.vertices = reshape(vertices,3,Np)';
        faces = fread(fid,Nf*3,'int32');
        Surf.SurfData.faces = reshape(faces,3,Nf)'+1;
        fclose(fid);
    else
        errordlg('This file is not a known Surface File. Please try again');
        return
    end
    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Np,Nf);
    Temp = sum(Tri);
    Tri(:,Temp==0) = [];
    Surf.Tri = Tri;
    Surf.Name = nm;
end

