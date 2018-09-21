function [cads] = MeanFA_SpamVecs(InputFolder, Id,FreeSDir,ConnectDir);
%
% Syntax :
% OutAtlasFile = GSpam2DiffRoi(GvecFolders,ArgFile,Imfile,Famap,Transf);
%
% Script file to unify giry surfaces in Brainvisa Tmtktri file.
%
% Input Parameters:
%   Lspamvecs         :  Left Gyral spam files
%   Rspamvecs         :  Right Gyral spam files
%   Imfile            : Image file in FreeSurfer space
%   Famap             : FA Map filename
%   Transf            : Transformation from T1 Space to Fa Space
%   ttype             : Transformation type(ie Fa2Gy: Fa Map to T1 space, 
%                                    Gy2Fa: Gyral Spam to FA Space)
%  
%
% Output Parameters:
%
%   OutAtlasFile      : Output Atlas File of Gyral Spam structures.
%      Oufiles        : 1st row FA, 2nd Gyral Spam Atlas.
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
% April 27th 2012
% Version $1.0

warning off;
% if nargin <6
%     [pth, nm, ext] = fileparts(Imfile); 
%     OutAtlasFile = [pth nm(1:end-3) '_GyralSpam_Atlas.nii'];
% end
%  InputFolder = '/media/Data/Joost/gspan/gyvecs/gyri';
%  Id = '0923x-20101105';
%  FreeSDir = '/media/Data/PROCESSING_RESULTS/PEPS/5-freesurfer_processing';
%  ConnectDir = '/media/Data/PROCESSING_RESULTS/PEPS/7-connectome';

ttype = 'Fa2Gy';
%Id = '0054x-20100604'
OutAtlasFile = [FreeSDir filesep Id filesep 'tmp' filesep 'tempatlas.nii' ];
Famap = [ConnectDir filesep Id filesep 'dtreconst' filesep 'fsl' filesep Id '_fa.nii'];
Transf =  [ConnectDir filesep Id filesep 'transforms' filesep Id '_T2Bet_2_diffAffine.txt'];

ctab = [FreeSDir filesep Id filesep 'label' filesep 'aparc.annot.ctab' ];
%ctab = ['/home/yaleman/lobes.ctab' ];
%if ~exist( [FreeSDir filesep Id filesep 'tmp' filesep 'aparc+aseg.nii'],'file')
    cad = ['mri_convert -i ' FreeSDir filesep Id filesep 'mri' filesep 'aparc+aseg.mgz  -o ' FreeSDir filesep Id filesep 'tmp' filesep 'aparc+aseg.nii'];
    system(cad);
%end
talfile = [FreeSDir filesep Id filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
cras1 = textread(talfile,'%s',5,'headerlines',20);
cras = char(cras1);
cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];

Imfile = [FreeSDir filesep Id filesep 'tmp' filesep 'aparc+aseg.nii'];
[sid,snames,r,g,b,tr] = textread(ctab,'%u%s%u%u%u%u','delimiter',' ');
 snames([1 5 36]) = [];
sid([1 5 36]) = [];
%% ====================== Left Hemisphere ================================%
Lspamvecs = Select_Files(InputFolder,{[Id '_L_gyrus*']},{[Id '_L_gyrus']});
% ----------------- Reading White Matter Spams ---------------------------%
for i = 1:size(Lspamvecs,1)
    name = deblank(Lspamvecs(i,:));
    ind = strfind(name,'_');
    ord(i) = str2num(name(ind(6)+1:ind(7)-1));
end
temp(ord,:) = Lspamvecs;
temp(4,:) = [];
Lspamvecs = temp;

cont = 0; inds = 0;
for i = 1:size(Lspamvecs,1)
    try
        Surfgs(i) = load_mesh_lines(deblank(Lspamvecs(i,:)));
        Surfgs(i).Name = ['LH_' char(snames(i))];
    catch
        inds = [inds;i];
    end
end

% ----------------- End of Reading White Matter Spams --------------------%
% ----------- Changing from Brainvisa Space to freesurfer Space ----------%
Surfgs = freesCS2brainvisaCS(Surfgs,Imfile,'b2f');

V = spm_vol(Imfile);
Ns = length(Surfgs);
N = 10;
Ilab = zeros(V.dim(1),V.dim(2),V.dim(3));
for i = 1:Ns
    if ~isempty(Surfgs(i).SurfData.vertices)
        %Surfgs(i).SurfData.vertices = Surfgs(i).SurfData.vertices + repmat(cras,[size(Surfgs(i).SurfData.vertices,1) 1]);
        vertvox = (inv(V.mat)*[Surfgs(i).SurfData.vertices ones(size(Surfgs(i).SurfData.vertices),1)]')';
        Surfgs(i).SurfData.vertices = vertvox(:,1:3);
        Vect = Surfgs(i).SurfData.vertices(Surfgs(i).SurfData.faces(:,2),:)-Surfgs(i).SurfData.vertices(Surfgs(i).SurfData.faces(:,1),:);
        Cp = Surfgs(i).SurfData.vertices(Surfgs(i).SurfData.faces(:,2),:);
        Cw = Surfgs(i).SurfData.vertices(Surfgs(i).SurfData.faces(:,1),:);
        Vect = Cp-Cw; dista = (sqrt((Cp(:,1)-Cw(:,1)).^2+(Cp(:,2)-Cw(:,2)).^2+(Cp(:,3)-Cw(:,3)).^2));Normals = Vect./repmat(dista,[1 3]);
        steps = [0:1/N:1];dista = repmat(steps,[size(dista,1) 1]).*repmat(dista,[1 N+1]);
        CintX = [repmat(Cw(:,1),[1 size(dista,2)])+repmat(Normals(:,1),[1 size(dista,2)]).*dista ]';siz = size(CintX);CintX = CintX(:);
        CintY = [repmat(Cw(:,2),[1 size(dista,2)])+repmat(Normals(:,2),[1 size(dista,2)]).*dista ]';CintY = CintY(:);
        CintZ = [repmat(Cw(:,3),[1 size(dista,2)])+repmat(Normals(:,3),[1 size(dista,2)]).*dista ]';CintZ = CintZ(:);
        ind = sub2ind(size(Ilab),round(CintX),round(CintY),round(CintZ));
        %inds = strfind(names(i,:),'_'); lab = str2num(deblank(names(i,inds+1:end)));
        Ilab(unique(ind)) = sid(i)+3000;
    end
    
end
%% =================== End of Left Hemisphere ============================%
%% ====================== Right Hemisphere ================================%
Rspamvecs = Select_Files(InputFolder,{[Id '_R_gyrus*']},{[Id '_R_gyrus']});
% ----------------- Reading White Matter Spams ---------------------------%
for i = 1:size(Rspamvecs,1)
    name = deblank(Rspamvecs(i,:));
    ind = strfind(name,'_');
    ord(i) = str2num(name(ind(6)+1:ind(7)-1));
end
temp(ord,:) = Rspamvecs;
temp(4,:) = [];
Rspamvecs = temp;

cont = 0; inds = 0;
for i = 1:size(Rspamvecs,1)
    try
        Surfgs(i) = load_mesh_lines(deblank(Rspamvecs(i,:)));
        Surfgs(i).Name = ['RH_' char(snames(i))];
    catch
        inds = [inds;i];
    end
end

% ----------------- End of Reading White Matter Spams --------------------%
% ----------- Changing from Brainvisa Space to freesurfer Space ----------%
Surfgs = freesCS2brainvisaCS(Surfgs,Imfile,'b2f');

V = spm_vol(Imfile);
Ns = length(Surfgs);
N = 10;
%Ilab = zeros(V.dim(1),V.dim(2),V.dim(3));
for i = 1:Ns
    if ~isempty(Surfgs(i).SurfData.vertices)
        %Surfgs(i).SurfData.vertices = Surfgs(i).SurfData.vertices + repmat(cras,[size(Surfgs(i).SurfData.vertices,1) 1]);
        vertvox = (inv(V.mat)*[Surfgs(i).SurfData.vertices ones(size(Surfgs(i).SurfData.vertices),1)]')';
        Surfgs(i).SurfData.vertices = vertvox(:,1:3);
        Vect = Surfgs(i).SurfData.vertices(Surfgs(i).SurfData.faces(:,2),:)-Surfgs(i).SurfData.vertices(Surfgs(i).SurfData.faces(:,1),:);
        Cp = Surfgs(i).SurfData.vertices(Surfgs(i).SurfData.faces(:,2),:);
        Cw = Surfgs(i).SurfData.vertices(Surfgs(i).SurfData.faces(:,1),:);
        Vect = Cp-Cw; dista = (sqrt((Cp(:,1)-Cw(:,1)).^2+(Cp(:,2)-Cw(:,2)).^2+(Cp(:,3)-Cw(:,3)).^2));Normals = Vect./repmat(dista,[1 3]);
        steps = [0:1/N:1];dista = repmat(steps,[size(dista,1) 1]).*repmat(dista,[1 N+1]);
        CintX = [repmat(Cw(:,1),[1 size(dista,2)])+repmat(Normals(:,1),[1 size(dista,2)]).*dista ]';siz = size(CintX);CintX = CintX(:);
        CintY = [repmat(Cw(:,2),[1 size(dista,2)])+repmat(Normals(:,2),[1 size(dista,2)]).*dista ]';CintY = CintY(:);
        CintZ = [repmat(Cw(:,3),[1 size(dista,2)])+repmat(Normals(:,3),[1 size(dista,2)]).*dista ]';CintZ = CintZ(:);
        ind = sub2ind(size(Ilab),round(CintX),round(CintY),round(CintZ));
        %inds = strfind(names(i,:),'_'); lab = str2num(deblank(names(i,inds+1:end)));
        Ilab(unique(ind)) = sid(i)+4000;
    end
    
end
%% =================== End of Left Hemisphere ============================%
%% =================== Saving Gyral Spam Atlas ===========================%
Nc = size(Ilab,3);
for i = 1:Nc
    T = squeeze(Ilab(:,:,i));
    Tf = imfill(T,4,'holes');
    Ilab(:,:,i) = Tf;
end
Nc = size(Ilab,1);
for i = 1:Nc
    T = squeeze(Ilab(i,:,:));
    %Tf = BWMORPH(T,'fill',Inf);
    Tf = imfill(T,4,'holes');
    Ilab(i,:,:) = Tf;
end
Nc = size(Ilab,2);
for i = 1:Nc
    T = squeeze(Ilab(:,i,:));
    %Tf = BWMORPH(T,'fill',Inf);
    Tf = imfill(T,4,'holes');
    Ilab(:,i,:) = Tf;
end

V = spm_vol(Imfile);
Vol = V;

Vol.fname = OutAtlasFile;
spm_write_vol(Vol,Ilab);
%% ================= End of Saving Gyral Spam Atlas ====================%
switch ttype
    case 'Fa2Gy'
        %% =========== Saving Results (FA Maps to T1 Space) ======================%
        Outfile =  [ConnectDir filesep Id filesep 'tmp' filesep 'Temp_Image.nii'];
        cad = ['WarpImageMultiTransform 3 ' Famap ' ' Outfile ' -R ' OutAtlasFile ' -i ' Transf ' --use-BSpline'];
        system(cad);
        Oufiles = strvcat(Outfile,OutAtlasFile);
    case 'Gy2Fa'
        %% =========== Saving Results (Gyral Spams to FA Maps Space) =============%
        Outfile =  [ConnectDir filesep Id filesep 'tmp' filesep 'Temp_Image.nii'];
        cad = ['WarpImageMultiTransform 3 ' OutAtlasFile ' ' Outfile ' -R ' Famap ' ' Transf ' --use-NN'];
        system(cad);
        Oufiles = strvcat(Famap,Outfile);
end

Vfa = spm_vol(Outfile);
Ifa = spm_read_vols(Vfa);

LHs = strvcat('LH-MeanFA-','RH-MeanFA-');
temp = repmat(snames,[1 2])';
aaaa = [repmat(LHs,[size(temp(:),1)/2 1]) char(temp(:)) repmat(';',[size(temp(:),1) 1])];
aba = aaaa';aba = aba(:)';aba(isspace(aba))=[];%varnames=char2cell(aba,{':'});
valuescad = '';
for i = 1:length(sid)
    indl = find(Ilab == sid(i)+3000);
    indr = find(Ilab == sid(i)+4000);
    if ~isempty(indl)
        a = mean(Ifa(indl));
%         b = std(Ifa(indl));
    else
        a = 0;
        b = 0;
    end
        if ~isempty(indr)
        c = mean(Ifa(indr));
%         d = std(Ifa(indr));
    else
        c = 0;
        d = 0;
    end
%     valuescad = [valuescad ';' num2str(a) ';' num2str(b) ';' num2str(c) ';' num2str(d)];
 valuescad = [valuescad ';' num2str(a) ';' num2str(c)];
end
% colnames = 
cads = ['Subject_ID;' aba];
cads = strvcat(cads,[Id  valuescad]);
delete(OutAtlasFile);
delete(Outfile);
return;

