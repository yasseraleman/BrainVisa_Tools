function [OutAtlasFile,Oufiles] = GSpam2DiffRoi(Lspamvecs,Rspamvecs,ArgFiles,Imfile,Famap,Transf,ttype,OutAtlasFile);
%
% Syntax :
% OutAtlasFile = GSpam2DiffRoi(GvecFolders,ArgFile,Imfile,Famap,Transf);
%
% Script file to unify giry surfaces in Brainvisa Tmtktri file.
%
% Input Parameters:
%   Lspamvecs         :  Left Gyral spam files
%   Rspamvecs         :  Right Gyral spam files
%   ArgFiles          :  ArgFile containing the surface relationships in
%                        Tmtktri(1st row. Left hemisphere, 2nd row Right Hemisphere)
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
% ArgFiles =  strvcat('/media/COSAS/Test/Joost/aquivan/aquivienen/ASPER_00001__101-20060510_FS2BV.lh.white+tal.aims_Tmtkmtri.arg','/media/COSAS/Test/Joost/aquivan/aquivienen/ASPER_00001__101-20060510_FS2BV.rh.white+tal.aims_Tmtkmtri.arg'); %Mandatory 
% Imfile = '/media/COSAS/Test/Joost/aquivan/aquivienen/ASPER_00001__101-20060510_T1.nii'; %Mandatory 
% Famap = '/media/Data/PROCESSING_RESULTS/ASPERGER/7-connectome/ASPER_00001__101-20060510/dtreconst/fsl/ASPER_00001__101-20060510_fa.nii'; %Mandatory
% Transf = '/media/Data/PROCESSING_RESULTS/ASPERGER/7-connectome/ASPER_00001__101-20060510/transforms/ASPER_00001__101-20060510_T2Bet_2_diffAffine.txt';
if nargin <7
    [pth, nm, ext] = fileparts(Imfile); 
    OutAtlasFile = [pth nm(1:end-3) '_GyralSpam_Atlas.nii'];
end
%% ====================== Left Hemisphere ================================%
ArgFile = deblank(ArgFiles(1,:));
% ----------------------- Reading Arg File -------------------------------%

fio = fopen(ArgFile,'rt');lines = '';cont = 0;
conts = 0;
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    lines = strvcat(lines,line);
end
fclose(fio);
names = lines(ismember(lines(:,1:4),'name','rows'),:);names = names(:,15:end);names(:,find(sum(isspace(names))==size(names,1))) = [];
labels = lines(ismember(lines(:,1:14),'Tmtktri_label','rows'),:);labels = labels(:,15:end);labels(:,find(sum(isspace(labels))==size(labels,1))) = [];

% ------------------- End of reading Arg File ----------------------------%
% ----------------- Reading White Matter Spams ---------------------------%

cont = 0; inds = 0;
for i = 1:size(names,1)
    spamvecst = Lspamvecs';spamvecst = spamvecst(:)';
    ind = strfind(spamvecst,['_' deblank(names(i,:)) '_']);[row, col] = ind2sub([size(Lspamvecs,2) size(Lspamvecs,1)],ind);
    try
        cont = cont+1;
        Surfgs(cont) = load_mesh_lines(deblank(Lspamvecs(col,:)));
        Surfgs(cont).Name = deblank(names(i,:));
    catch
        inds = [inds;i];
    end
end
inds(1) = [];
names(inds,:) = [];
labels(inds,:) = [];

% ----------------- End of Reading White Matter Spams --------------------%
% ----------- Changing from Brainvisa Space to freesurfer Space ----------%
Surfgs = freesCS2brainvisaCS(Surfgs,Imfile,'b2f');
V = spm_vol(Imfile);
Ns = length(Surfgs);
N = 10;
Ilab = zeros(V.dim(1),V.dim(2),V.dim(3));
for i = 1:Ns
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
     inds = strfind(names(i,:),'_'); lab = str2num(deblank(names(i,inds+1:end)));
     Ilab(unique(ind)) = lab;
end
%% =================== End of Left Hemisphere ============================%
%% ===================== Right Hemisphere ================================%
ArgFile = deblank(ArgFiles(2,:));
% ----------------------- Reading Arg File -------------------------------%

fio = fopen(ArgFile,'rt');lines = '';cont = 0;
conts = 0;
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    lines = strvcat(lines,line);
end
fclose(fio);
names = lines(ismember(lines(:,1:4),'name','rows'),:);names = names(:,15:end);names(:,find(sum(isspace(names))==size(names,1))) = [];
labels = lines(ismember(lines(:,1:14),'Tmtktri_label','rows'),:);labels = labels(:,15:end);labels(:,find(sum(isspace(labels))==size(labels,1))) = [];

% ------------------- End of reading Arg File ----------------------------%
% ----------------- Reading White Matter Spams ---------------------------%

cont = 0; inds = 0;
for i = 1:size(names,1)
    spamvecst = Rspamvecs';spamvecst = spamvecst(:)';
    ind = strfind(spamvecst,['_' deblank(names(i,:)) '_']);[row, col] = ind2sub([size(Rspamvecs,2) size(Rspamvecs,1)],ind);
    try
        cont = cont+1;
        Surfrs(cont) = load_mesh_lines(deblank(Rspamvecs(col,:)));
        Surfrs(cont).Name = deblank(names(i,:));
    catch
        inds = [inds;i];
    end
end
inds(1) = [];
names(inds,:) = [];
labels(inds,:) = [];

% ----------------- End of Reading White Matter Spams --------------------%
% ----------- Changing from Brainvisa Space to freesurfer Space ----------%

Surfrs = freesCS2brainvisaCS(Surfrs,Imfile,'b2f');
Ns = length(Surfrs);
for i = 1:Ns
    vertvox = (inv(V.mat)*[Surfrs(i).SurfData.vertices ones(size(Surfrs(i).SurfData.vertices),1)]')';
    Surfrs(i).SurfData.vertices = vertvox(:,1:3);
     Vect = Surfrs(i).SurfData.vertices(Surfrs(i).SurfData.faces(:,2),:)-Surfrs(i).SurfData.vertices(Surfrs(i).SurfData.faces(:,1),:);
     Cp = Surfrs(i).SurfData.vertices(Surfrs(i).SurfData.faces(:,2),:);
     Cw = Surfrs(i).SurfData.vertices(Surfrs(i).SurfData.faces(:,1),:);
     Vect = Cp-Cw; dista = (sqrt((Cp(:,1)-Cw(:,1)).^2+(Cp(:,2)-Cw(:,2)).^2+(Cp(:,3)-Cw(:,3)).^2));Normals = Vect./repmat(dista,[1 3]);
     steps = [0:1/N:1];dista = repmat(steps,[size(dista,1) 1]).*repmat(dista,[1 N+1]);
     CintX = [repmat(Cw(:,1),[1 size(dista,2)])+repmat(Normals(:,1),[1 size(dista,2)]).*dista ]';siz = size(CintX);CintX = CintX(:);
     CintY = [repmat(Cw(:,2),[1 size(dista,2)])+repmat(Normals(:,2),[1 size(dista,2)]).*dista ]';CintY = CintY(:);
     CintZ = [repmat(Cw(:,3),[1 size(dista,2)])+repmat(Normals(:,3),[1 size(dista,2)]).*dista ]';CintZ = CintZ(:);
%      plot3(CintX,CintY,CintZ,'.','Markersize',10,'Color',Color(i,:)/255);
     ind = sub2ind(size(Ilab),round(CintX),round(CintY),round(CintZ));
     inds = strfind(names(i,:),'_'); lab = str2num(deblank(names(i,inds+1:end)));
     Ilab(unique(ind)) = lab+100;
     
     %Plot_Normals(Surfrs(i).SurfData.vertices(Surfrs(i).SurfData.faces(:,1),:), Vect,0);
end
%% =================== End of Right Hemisphere ============================%   
%% =================== Saving Gyral Spam Atlas ===========================%
V = spm_vol(Imfile);
Vol = V;
Vol.fname = OutAtlasFile;
spm_write_vol(Vol,Ilab);
%% ================= End of Saving Gyral Spam Atlas ====================%
switch ttype
    case 'Fa2Gy'
        %% =========== Saving Results (FA Maps to T1 Space) ======================%
        [pth, nm, ext] = fileparts(Imfile);
        mkdir(pth,'tmp');
        Outfile =  [pth filesep 'tmp' filesep 'Temp_Image.nii'];
        [pth, nm, ext] = fileparts(Famap);
        cad = ['WarpImageMultiTransform 3 ' Famap ' ' Outfile ' -R ' OutAtlasFile ' -i ' Transf ' --use-BSpline'];
        system(cad);
        Oufiles = strvcat(Outfile,OutAtlasFile);
    case 'Gy2Fa'
        %% =========== Saving Results (Gyral Spams to FA Maps Space) =============%
        [pth, nm, ext] = fileparts(Imfile);
        mkdir(pth,'tmp');
        Outfile =  [pth filesep 'tmp' filesep 'Temp_Image.nii'];
        cad = ['WarpImageMultiTransform 3 ' OutAtlasFile ' ' Outfile ' -R ' Famap ' ' Transf ' --use-NN'];
        system(cad);
        Oufiles = strvcat(Famap,Outfile);
end
return;

