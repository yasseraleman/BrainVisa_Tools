function OutTextFile = Import_FreeSurfer_into_BrainVisa_Inputs_and_Outputs(FreesDir, FreesDatabaseDir, IdFile);
%
% Syntax :
%    Import_FreeSurfer_into_BrainVisa(FreesDir, FreesDatabaseDir, IdFile);
%
% This script generates the brainvisa inputs and outputs to facilitate the importation process.
%
% Input Parameters:
%     FreeSurferDatabaseDir          : FreeSurfer Database directory.
%     BrainVisaDatabaseDir           : Brainvisa Database Directory.
%     IdFile                         :Ids File.
%
% Output Parameters:
%    OutTextFile          : Output Text File.
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
% August 8th 2012
% Version $1.0

%% ======================== Main Program =================================%

FreeSurferDatabaseDir = '/home/yaleman/FreesurferDatabase';
BrainVisaDatabaseDir = ' /home/yaleman/BrainVisaDataBase';
IdFile = '/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';
OutTextFile = [BrainVisaDatabaseDir filesep 'Config_Import_File.txt'];


Ids = char(textread(IdFile,'%s'));

Ids = strvcat('URNC_IVI_F_17_POST1.long.IVI_F_17_template','URNC_IVI_F_17_PRE.long.IVI_F_17_template','URNC_IVI_M_17_POST1.long.IVI_M_17_template','URNC_IVI_M_17_PRE.long.IVI_M_17_template');
% Ids = '1SUBJECT_100206 _TestSubject';

% Ids = Ids([401:end],:);
Ns = size(Ids,1);
cadmri = '';
cadribbon = '';
cadtala = '';

cadt1out= '';
cadrefe= '';
cadsca= '';
cadbias= '';
cadnorm= '';
cadacpct= '';
cadacpc= '';
cadhan= '';
cadbrain = '';
cadsegm= '';
cadrgrey= '';
cadlgrey= '';
cadbf = '';
cadhf = '';
cadwr = '';
cadmc = '';
cadvar = '';
cadedg = '';

for i = 1:Ns
    Id = deblank(Ids(i,:));
    disp(['Preparing Subject: ' Id '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    
    
    %% Inputs
    
    % MRI orig
    ncad = ['''' FreeSurferDatabaseDir filesep Id filesep 'mri' filesep 'orig.mgz' '''' ' '];
    cadmri= [cadmri ncad];
    
    %     '/home/yaleman/FreeSurferDataBase/HCP_100408-20140807-T1wMPR1/mri/orig.mgz' '/home/yaleman/FreeSurferDataBase/HCP_101107-20140807-T1wMPR1/mri/orig.mgz';
    %
    %     % MRI ribbon
    ncad = ['''' FreeSurferDatabaseDir filesep Id filesep 'mri' filesep 'ribbon.mgz' '''' ' '];
    cadribbon= [cadribbon ncad];
    %     '/home/yaleman/FreeSurferDataBase/HCP_100408-20140807-T1wMPR1/mri/orig.mgz' '/home/yaleman/FreeSurferDataBase/HCP_101107-20140807-T1wMPR1/mri/ribbon.mgz';
    %
    %     % Talairach
    ncad = ['''' FreeSurferDatabaseDir filesep Id filesep 'mri' filesep 'transforms' filesep 'talairach.auto.xfm' '''' ' '];
    cadtala= [cadtala ncad];
    %     '/home/yaleman/FreeSurferDataBase/HCP_101107-20140807-T1wMPR1/mri/transforms/talairach.auto.xfm' '/home/yaleman/FreeSurferDataBase/HCP_101107-20140807-T1wMPR2/mri/transforms/talairach.auto.xfm'
    %
%     
%     
%     %% Outputs
%     % T1_outputç
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep Id '.nii.gz' '''' ' '];
    cadt1out= [cadt1out ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/HCP_100307-20140807-T1wMPR1.nii.gz
%     
%     % T1_referential    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'registration' filesep 'RawT1-' Id '_default_acquisition.referential' '''' ' '];
    cadrefe= [cadrefe ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/registration/RawT1-HCP_100307-20140807-T1wMPR1_default_acquisition.referential
%     
%     % Transform_to_scanner_based
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'registration' filesep 'RawT1-' Id '_default_acquisition_TO_Scanner_Based.trm' '''' ' '];
    cadsca= [cadsca ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/registration/RawT1-HCP_100307-20140807-T1wMPR1_default_acquisition_TO_Scanner_Based.trm
%     
%     % Bias_corrected_output
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'nobias_' Id '.nii.gz' '''' ' '];
    cadbias= [cadbias ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/default_analysis/nobias_HCP_100307-20140807-T1wMPR1.nii.gz
%     
%     % Normalization_transformation
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'registration' filesep 'RawT1-' Id '_default_acquisition_TO_Talairach-MNI.trm' '''' ' '];
    cadnorm= [cadnorm ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/registration/RawT1-HCP_100307-20140807-T1wMPR1_default_acquisition_TO_Talairach-MNI.trm
%     
%     % Talairach_transform
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'registration' filesep 'RawT1-' Id '_default_acquisition_TO_Talairach-ACPC.trm' '''' ' '];
    cadacpct= [cadacpct ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/registration/RawT1-HCP_100307-20140807-T1wMPR1_default_acquisition_TO_Talairach-ACPC.trm
%     
%     % Commisure_coordinates
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep  Id '.APC' '''' ' '];
    cadacpc= [cadacpc ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/HCP_100307-20140807-T1wMPR1.APC
%     
%     % Histo_analysis
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/default_analysis/nobias_HCP_100307-20140807-T1wMPR1.han
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'nobias_' Id '.han' '''' ' '];
    cadhan= [cadhan ncad];
%     
%     % Split_brain_output
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'brain_' Id '.nii.gz' '''' ' '];
    cadbrain= [cadbrain ncad];

%     % Split_brain_output
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'voronoi_' Id '.nii.gz' '''' ' '];
    cadsegm= [cadsegm ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/default_analysis/segmentation/voronoi_HCP_100307-20140807-T1wMPR1.nii.gz
%     
%     % Rgrey_white_output
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Rgrey_white_' Id '.nii.gz' '''' ' '];
    cadrgrey= [cadrgrey ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/default_analysis/segmentation/Rgrey_white_HCP_100307-20140807-T1wMPR1.nii.gz
%     
%     % Lgrey_white_output
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Lgrey_white_' Id '.nii.gz' '''' ' '];
    cadlgrey= [cadlgrey ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/default_analysis/segmentation/Lgrey_white_HCP_100307-20140807-T1wMPR1.nii.gz

%     % Bias Field
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'biasfield_' Id '.nii.gz' '''' ' '];
    cadbf= [cadbf ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/100206/t1mri/default_acquisition/default_analysis/biasfield_100206.nii.gz

    
    %     % H Filtered 
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'hfiltered_' Id '.nii.gz' '''' ' '];
    cadhf= [cadhf ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/100206/t1mri/default_acquisition/default_analysis/hfiltered_100206.nii.gz

    
    %     % White Ridge
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'whiteridge_' Id '.nii.gz' '''' ' '];
    cadwr= [cadwr ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/100206/t1mri/default_acquisition/default_analysis/whiteridge_100206.nii.gz

    
    %     % Mean Curvature
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'mean_curvature_' Id '.nii.gz' '''' ' '];
    cadmc= [cadmc ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/100206/t1mri/default_acquisition/default_analysis/mean_curvature_100206.nii.gz

    
    %     % Variance
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'variance_' Id '.nii.gz' '''' ' '];
    cadvar= [cadvar ncad];
% /home/yaleman/BrainVisaDataBase/subjects/100206/t1mri/default_acquisition/default_analysis/variance_100206.nii.gz

    %     % Edges
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'edges_' Id '.nii.gz' '''' ' '];
    cadedg= [cadedg ncad];


end
cadmri(end) = [];
cadribbon(end) = [];
cadtala(end) = [];

cadt1out(end) = [];
cadrefe(end) = [];
cadsca(end) = [];
cadbias(end) = [];
cadnorm(end) = [];
cadacpct(end) = [];
cadacpc(end) = [];
cadhan(end) = [];
cadbrain(end) = [];
cadsegm(end) = [];
cadrgrey(end) = [];
cadlgrey(end) = [];

cadbf(end) = [];
cadhf(end) = [];
cadwr(end) = [];
cadmc(end) = [];
cadvar(end) = [];
cadedg(end) = [];

cads = strvcat('## =========== Inputs ===============','# MRI orig', cadmri,' ' ,'# MRI ribbon' ,cadribbon, ' ' , '# Talairach Transform', cadtala, ' ' ,' ' , '## =========== Outputs ===============', '# T1 output' ,cadt1out, ' ' , '# T1 referential', cadrefe, ' ' , '# Transform to scanner_based' , cadsca, ...
    ' ' , '# Bias corrected output', cadbias, ' ' , '# Normalization transformation' ,cadnorm, ' ' , '# Talairach transform', cadacpct, ' ' , '# Commisure coordinates', cadacpc, ' ' , '# Histo analysis', cadhan, ' ' , '# Brain output', cadbrain, ' ' , '# Split brain output', cadsegm,...
    ' ' , '# Rgrey white output', cadrgrey, ' ' , '# Lgrey white output', cadlgrey, ' ', '# Bias field output', cadbf, ' ' , '# H Filtered', cadhf,' ', '# White Ridges', cadwr,' ', '# Mean Curvature', cadmc,' ', '# Variance', cadvar,' ', '# Edges', cadedg);
fid = fopen('/home/yaleman/test.txt','wt');
for i = 1:size(cads,1)
    fprintf(fid,'%s\n',deblank(cads(i,:)));
end
fclose(fid);

return;