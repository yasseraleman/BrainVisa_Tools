close all;
V = spm_vol('/media/COSAS/Test/freesurfer/ch2/tmp/T1.nii');
I = spm_read_vols(V);
Surf= read_surfreesurfer('/media/COSAS/Test/freesurfer/ch2/surf/lh.sphere');
txt = read_cfiles('/media/COSAS/Test/freesurfer/ch2/surf/lh.curv');
Surf.Is = txt;
Plot_Surf(Surf);
norma = sqrt(sum(Surf.SurfData.vertices'.^2))';
ind = find(txt > 0);
indd = find(sum(ismember(Surf.SurfData.faces,ind)') >0);
Surfd = Surf;
Surf.SurfData.vertices = Surf.SurfData.vertices./repmat(norma,[1 3]) - 0.5*[txt txt txt].*Surf.SurfData.vertices./repmat(norma,[1 3]);

Surfd.SurfData.faces(indd,:) = [];
Plot_Surf(Surf);
Plot_Surf(Surfd);

orig = repmat([0.5000 -16.5000 19.5000],[size(Surf.SurfData.vertices,1) 1]);
Surf.SurfData.vertices = Surf.SurfData.vertices+orig;
% vertvox = (inv(V.mat)*[Surf.SurfData.vertices ones(size(Surf.SurfData.vertices),1)]')';
% Surf.SurfData.vertices = vertvox(:,1:3);
Surfa(1) = Surf;

Surf= read_surfreesurfer('/media/COSAS/Test/freesurfer/ch2/surf/lh.white');
orig = repmat([0.5000 -16.5000 19.5000],[size(Surf.SurfData.vertices,1) 1]);
Surf.SurfData.vertices = Surf.SurfData.vertices+orig;
% vertvox = (inv(V.mat)*[Surf.SurfData.vertices ones(size(Surf.SurfData.vertices),1)]')';
% Surf.SurfData.vertices = vertvox(:,1:3);
Surfa(2) = Surf;


hf = Image_plus_surface_nosave('/media/COSAS/Test/freesurfer/ch2/tmp/fT1.nii',Surfa, [0 0 120]);

% 
showcs3(I);
strsurf=patch(Surf.SurfData,'edgecolor','none', 'tag','patch','facelighting','gouraud');
camlight
view([90 0])