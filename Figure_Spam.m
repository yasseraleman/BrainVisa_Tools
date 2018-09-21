%% For Aparc
% 
% close all; clear all;
% Inputdir = '/media/Data/yasser/ch_brainvisa';
% Surf = load_mesh('/media/Data/yasser/ch_brainvisa/ch2_L_gyralWM_skeleton.mesh');
% Surf = freesCS2brainvisaCS(Surf,'/media/COSAS/Test/freesurfer/ch2/tmp/T1.nii','b2f');
% [txt,ctab] = read_cfiles('/media/COSAS/Test/freesurfer/ch2/label/lh.aparc.annot');
% % figure;
% % strsurf=patch(Surf.SurfData,'facecolor',[1 1 1],'edgecolor','none','tag', 'model0','facelighting','gouraud');
% % set(strsurf,'FaceAlpha',1);
% % axis image;view(3);camlight;
% Plot_Surf(Surf);axis off;box off;
% a = dir([Inputdir filesep '*_wmspan.mesh']);
% if ~isempty(a)
%     [names{1:size(a,1),1}]=deal(a.name);
%     SpamFiles = [repmat([deblank(Inputdir) filesep],[size(names,1) 1]) char(names)];
% end
% for i = 1:size(SpamFiles,1)
%     Snames = deblank(SpamFiles(i,:));
%     [pth,name,ext] = fileparts(Snames);
%     names = char(ctab.struct_names);
%     ind = strfind(name,'_');
%     sname = name(ind(2)+1:ind(3)-1);
%     ind = find(ismember(names(:,1:length(sname)),sname,'rows'));
%     col = ctab.table(ind,1:3)/255;
%     hold on;
%     Surf1 = load_mesh(deblank(SpamFiles(i,:)));
%     Surf1 = freesCS2brainvisaCS(Surf1,'/media/COSAS/Test/freesurfer/ch2/tmp/T1.nii','b2f');
%     V1 = Surf1.SurfData.vertices(Surf1.SurfData.faces(:,1),:);
%     V2 = Surf1.SurfData.vertices(Surf1.SurfData.faces(:,2),:);
%     line([V1(:,1) V2(:,1)]',[V1(:,2) V2(:,2)]',[V1(:,3) V2(:,3)]','Color',col);
% end

%% For Lobes

close all; clear all;
Inputdir = '/media/Data/yasser/ch_brainvisa/Lobes';
Surf = load_mesh('/media/Data/yasser/ch_brainvisa/ch2_L_gyralWM_skeleton.mesh');
Surf = freesCS2brainvisaCS(Surf,'/media/COSAS/Test/freesurfer/ch2/tmp/T1.nii','b2f');
cols = [74 95 148 ;108 177 210;219 171 145;84 176 113]; 
% figure;
% strsurf=patch(Surf.SurfData,'facecolor',[1 1 1],'edgecolor','none','tag', 'model0','facelighting','gouraud');
% set(strsurf,'FaceAlpha',1);
% axis image;view(3);camlight;
Plot_Surf(Surf);axis off;box off;
a = dir([Inputdir filesep '*_wmspan.mesh']);
if ~isempty(a)
    [names{1:size(a,1),1}]=deal(a.name);
    SpamFiles = [repmat([deblank(Inputdir) filesep],[size(names,1) 1]) char(names)];
end
for i = 1:size(SpamFiles,1)
    Snames = deblank(SpamFiles(i,:));
    [pth,name,ext] = fileparts(Snames);
    col = cols(i,:)/255;
    hold on;
    Surf1 = load_mesh(deblank(SpamFiles(i,:)));
    Surf1 = freesCS2brainvisaCS(Surf1,'/media/COSAS/Test/freesurfer/ch2/tmp/T1.nii','b2f');
    V1 = Surf1.SurfData.vertices(Surf1.SurfData.faces(:,1),:);
    V2 = Surf1.SurfData.vertices(Surf1.SurfData.faces(:,2),:);
    line([V1(:,1) V2(:,1)]',[V1(:,2) V2(:,2)]',[V1(:,3) V2(:,3)]','Color',[ 1 0 0]);
end