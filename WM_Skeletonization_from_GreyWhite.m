function [OutFile] = Gyral_WM_Skeletonization(Id);
%
% Syntax :
% [OutFile, cads] = Gyral_WM_Skeletonization(FreeSDir, Id, OutAtlasFile);
%
% This script skeletonize the fa image and label the FA inside each gyri
% using the FreeSurfer cortical parcellation.
%
% Input Parameters:
%   OutputDir         : Pipeline Output directory 
%   Id                : Subject Id
%  
%
% Output Parameters:
%
%     OutFile         : Labeled skeleton
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
% July 27th 2013
% Version $1.0

% OutputDir = '/media/Data/PROCESSING_RESULTS/PEPS';
% OutAtlasFile = '/media/Data/Joost/gspan/volumes/gyri/0019x-20091120.wmparc.gyri.nii';


AtlasFile = '/media/MyDisk/Test/subjects/1Test_HCP_899885-20140807-T1wMPR1/t1mri/default_acquisition/RII/Lgrey_white_1Test_HCP_899885-20140807-T1wMPR1.nii';

cad = ['mri_convert /media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing/1Test_HCP_899885-20140807-T1wMPR1/mri/ribbon.mgz ' TempDir filesep 'ribbon.nii'];
system(cad);


% try
%     cad = ['gunzip -df ' ImFile];
%     system(cad);
% end

%AtlasFile = [FreeSDir filesep Id filesep 'tmp' filesep AtlasId '.nii'];
TempDir = '/media/MyDisk/Test/subjects/1Test_HCP_899885-20140807-T1wMPR1/t1mri/default_acquisition/RII';
mkdir(TempDir);


V = spm_vol([TempDir filesep 'ribbon.nii']);
V.mat(1:3,4) = V.mat(1:3,4) - [-4.024276733398438e+00 4.005416870117188e+01 -6.548309326171875e-02]';
%      V.mat = [1     0     0  -128;
%      0     1     0  -129;
%      0     0     1  -129;
%      0     0     0     1];

V = spm_create_vol(V);




I = spm_read_vols(V);
ind = find(I <200);
I(ind) = 0;
V.fname =  [TempDir filesep 'wm.nii'];
spm_write_vol(V,logical(I));

 cad = ['fslmaths ' V.fname ' -s 2 ' TempDir filesep 'swm.nii'];
 system(cad);
try
    cad = ['gunzip -df ' TempDir filesep 'swm.nii.gz'];
    system(cad);
end

cad = ['tbss_skeleton -i ' TempDir filesep 'swm.nii' ' -o ' TempDir filesep 'swm_skel.nii'];
system(cad);
try
    cad = ['gunzip -df ' TempDir filesep 'swm_skel.nii.gz'];
    system(cad);
end


V = spm_vol([ FreeSDir filesep Id filesep 'tmp' filesep Id '_swm_skel.nii']);
I = spm_read_vols(V);
VA = spm_vol(AtlasFile);
IA = spm_read_vols(VA);
IA = logical(I).*IA;
ind = find(IA<3000);
IA(ind) = 0;
ind = find(IA>5000);
IA(ind) = 0;
VA.fname =  [FreeSDir filesep Id filesep 'mri' filesep Id '_skelwmatlas.nii'];
spm_write_vol(VA,IA);
delete([FreeSDir filesep Id filesep 'tmp' filesep Id '_wm.nii']);
delete([FreeSDir filesep Id filesep 'tmp' filesep Id '_swm.nii.gz']);
delete([ FreeSDir filesep Id filesep 'tmp' filesep Id '_swm_skel.nii']);
delete(AtlasFile);

OutFile = [FreeSDir filesep Id filesep 'mri' filesep Id '_skelwmatlas.nii'];
cad = ['mri_convert ' OutFile ' ' FreeSDir filesep Id filesep 'mri' filesep Id '_skelwmatlas.mgz'];
system(cad);
delete(OutFile);
OutFile = [FreeSDir filesep Id filesep 'mri' filesep Id '_skelwmatlas.mgz'];

return;
