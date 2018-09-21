function My_WhiteMatterSPAM_From_Brainvisa(FreeSurferDatabaseDir,BrainVisaDatabaseDir,IdFile);

FreeSurferDatabaseDir = '/media/Data/PROCESSING_RESULTS/URNC/Longitudinal_analysis/Structural/5-freesurfer_processing/New_annot/subjects_long_temp/';
BrainVisaDatabaseDir = '/media/COSAS/URNC-BrainVISADataBase/';
IdFile = '/media/Data/PROCESSING_RESULTS/URNC/Longitudinal_analysis/Structural/5-freesurfer_processing/New_annot/subjects_long_temp/URNC_Long_IDs.txt';
if exist(IdFile,'file')
    Ids = char(textread(IdFile,'%s'));
else
    Ids = IdFile;
end

clust = parcluster('local');
numWorkers = clust.NumWorkers;

try
    parpool(numWorkers-2);
catch
    parpool(numWorkers-2);
end


%Ids = '1Test_HCP_899885-20140807-T1wMPR1';
Ns = size(Ids, 1);
for i = 23:28
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    %subjId = '1Test_HCP_899885-20140807-T1wMPR1';
    % ---- Creating Output Directory
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII']);
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'spams']);
    Outdir = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'spams'];
    
    %% %%%%%%%%%%%%%%%%%%%%%%%  Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Reading ArgFile
    LArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep subjId '_Lgyri_default_session_auto.arg' ];
    [Lines, StNames] = Read_GyriArgFiles(LArgFile);
    StNames(ismember(StNames(:,1:11),'medial_wall','rows'),:) = []; % Removing Meadial Wall from the structures list
    
    % ---- Skeleton Mesh
    LSkelMesh = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh' filesep 'L' subjId '_skel.mesh'];
    
    % ---- Computing Gyral White Matter Spams
    Nstruct = size(StNames,1);
  
    parfor j = 1:Nstruct
        strname = deblank(StNames(j,:));
        cad = ['RicGyralSpan -g ' LArgFile ' --sm ' LSkelMesh ' -n ' strname ' --sn ' subjId ' -o ' Outdir filesep subjId '_L --mt 12 --md 0.75' ];
        system(cad);
    end
    %% %%%%%%%%%%%%%%%%%%%%%  End of Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%  Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Reading ArgFile
    RArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep subjId '_Rgyri_default_session_auto.arg' ];
    [Lines, StNames] = Read_GyriArgFiles(RArgFile);
    StNames(ismember(StNames(:,1:11),'medial_wall','rows'),:) = []; % Removing Meadial Wall from the structures list
    
    % ---- Skeleton Mesh
    RSkelMesh = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh' filesep 'R' subjId '_skel.mesh'];
    
    % ---- Computing Gyral White Matter Spams
    Nstruct = size(StNames,1);
    parfor j = 1:Nstruct
        strname = deblank(StNames(j,:));
        cad = ['RicGyralSpan -g ' RArgFile ' --sm ' RSkelMesh ' -n ' strname ' --sn ' subjId ' -o ' Outdir filesep subjId '_R --mt 12 --md 0.75' ];
        system(cad);
    end
   %% %%%%%%%%%%%%%%%%%%%%%  End of Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%
end
delete(gcp);
return;



