function My_WhiteMatter_Skeleton_From_BrainVisa(InputImage, OutputDirectory);

% FreeSurferDatabaseDir = '/media/Data/PROCESSING_RESULTS/URNC/Longitudinal_analysis/Structural/5-freesurfer_processing/New_annot/subjects_long_temp/';
% BrainVisaDatabaseDir = '/media/COSAS/URNC-BrainVISADataBase/';
BrainVisaDatabaseDir = '/media/COSAS/8-BrainVISADataBase-HCP/';

% IdFile = '/media/Data/PROCESSING_RESULTS/URNC/Longitudinal_analysis/Structural/5-freesurfer_processing/New_annot/subjects_long_temp/URNC_Long_IDs.txt';
IdFile= '100206';
if exist(IdFile,'file')
Ids = char(textread(IdFile,'%s'));
else
    Ids = IdFile;
end
clust = parcluster('local');
numWorkers = clust.NumWorkers;

% try
%     parpool(numWorkers-2);
% catch
%     parpool(numWorkers-2);
% end

%Ids = Ids(1:2,:);
Ns = size(Ids,1);
for i = 1:Ns
    subjId = deblank(Ids(i,:));
%     if ~exist([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'spams' filesep subjId '_R_temporallateral_gyvec.mesh'],'file')&~exist([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'spams' filesep subjId '_L_temporallateral_gyvec.mesh'],'file')
% % % % % %     if ~exist([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh' filesep 'L' subjId '_skel.mesh'],'file')|~exist([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh' filesep 'R' subjId '_skel.mesh'],'file')

        %% Creating Directories
        if ~exist([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh'],'dir')
            mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII']);
            mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh']);
        end
        OutputDirectory = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh'];
        files2delete= '';
        
        %% Detecting CRAS
%         Taltransf = [FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
%         cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);
%         cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
        
        
        %% ===================== Left Hemisphere =============================== %%
        InputImage = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Lgrey_white_' subjId '.nii.gz'];
        system(['gunzip -d ' InputImage]);
        InputImage = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Lgrey_white_' subjId '.nii'];
        files2delete= strvcat(files2delete,InputImage);
        
        [OutImage,OutMesh] = Gyral_WM_Skeleton_Computation(InputImage, OutputDirectory);
        try
            movefile([OutputDirectory filesep 'final_skel.nii'], [OutputDirectory filesep 'L' subjId '_skel.nii']);
        end
        try
            movefile(OutMesh, [OutputDirectory filesep 'L' subjId '_skel.mesh']);
        end
        try
            movefile([OutputDirectory filesep 'final_skel.nii.minf'], [OutputDirectory filesep 'L' subjId '_skel.nii.minf']);
        end
        OutImageL = [OutputDirectory filesep 'L' subjId '_skel.nii'];
        %--- Adding RAs center
        VL = spm_vol(OutImageL);
        VL.mat(1:3,4) = VL.mat(1:3,4) + cras';
        VL = spm_create_vol(VL);
        
        %% %%%%%%%%%%%%%%%%%%%%%%%  Computing Gyral Spam %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ---- Creating Output Directory
        mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII']);
        mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'spams']);
        Outdir = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'spams'];
        
        
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
        
        
        
        
        %% =====================End of Left Hemisphere ========================= %%
        
        %% ===================== Right Hemisphere ============================== %%
% % % % % %         InputImage = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Rgrey_white_' subjId '.nii.gz'];
% % % % % %         system(['gunzip -d ' InputImage]);
% % % % % %         InputImage = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Rgrey_white_' subjId '.nii'];
% % % % % %         files2delete= strvcat(files2delete,InputImage);
% % % % % %         
% % % % % %         [OutImage,OutMesh] = Gyral_WM_Skeleton_Computation(InputImage, OutputDirectory);
% % % % % %         
% % % % % %         try
% % % % % %             movefile([OutputDirectory filesep 'final_skel.nii'], [OutputDirectory filesep 'R' subjId '_skel.nii']);
% % % % % %         end
% % % % % %         
% % % % % %         try
% % % % % %             movefile(OutMesh, [OutputDirectory filesep 'R' subjId '_skel.mesh']);
% % % % % %         end
% % % % % %         
% % % % % %         try
% % % % % %             movefile([OutputDirectory filesep 'final_skel.nii.minf'], [OutputDirectory filesep 'R' subjId '_skel.nii.minf']);
% % % % % %         end
% % % % % %         OutImageR = [OutputDirectory filesep 'R' subjId '_skel.nii'];
% % % % % %         %--- Adding RAs center
% % % % % %         VR = spm_vol(OutImageR);
% % % % % %         VR.mat(1:3,4) = VR.mat(1:3,4) + cras';
% % % % % %         VR = spm_create_vol(VR);
% % % % % %         
% % % % % %         
% % % % % %         %% %%%%%%%%%%%%%%%%%%%%%%%  Computing Gyral Spam %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %         % ---- Creating Output Directory
% % % % % %         mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII']);
% % % % % %         mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'spams']);
        Outdir = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'spams'];
        
        
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
        
        
        
        %% =====================End of Right Hemisphere ======================== %%
        
        % % % %     Vj = VR;
        % % % %     IL = spm_read_vols(VL);
        % % % %     IR = spm_read_vols(VR);
        % % % %     Vj.fname = [OutputDirectory filesep subjId '_skel.nii'];
        % % % %     spm_write_vol(Vj,logical(IL+IR));
        % % % %
        % % % %
        % % % %     cad = ['mri_convert ' FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep  'wmparc.mgz ' FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep  'wmparc.nii'];
        % % % %     system(cad);
        % % % %     AtlasFile = [FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep  'wmparc.nii'];
        % % % %     files2delete= strvcat(files2delete,AtlasFile);
        % % % %
        % % % %     cad = ['mri_convert -c -i ' OutputDirectory filesep subjId '_skel.nii -o '  FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'wmskel.nii'];
        % % % %     system(cad);
        % % % %     files2delete= strvcat(files2delete,[FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'wmskel.nii']);
        % % % %
        % % % %     V = spm_vol([FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'wmskel.nii']);
        % % % %     I = spm_read_vols(V);
        % % % %     VA = spm_vol(AtlasFile);
        % % % %     IA = spm_read_vols(VA);
        % % % %     IA = logical(I).*IA;
        % % % %     ind = find(IA<3000);
        % % % %     IA(ind) = 0;
        % % % %     ind = find(IA>5000);
        % % % %     IA(ind) = 0;
        % % % %     VA.fname =  [FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep subjId '_skelwmatlas.nii'];
        % % % %     spm_write_vol(VA,IA);
        % % % %
        % % % %     files2delete= strvcat(files2delete,VA.fname);
        % % % %
        % % % %     OutFile = [FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep subjId '_skelwmatlas.nii'];
        % % % %     cad = ['mri_convert ' OutFile ' ' FreeSurferDatabaseDir filesep subjId filesOut_meshep 'mri' filesep 'wm.skelatlas.BV.mgz'];
        % % % %     system(cad);
        % % % %
        % % % %     for i = 1:size(files2delete,1)
        % % % %         cad = ['rm -r ' deblank(files2delete(i,:))];
        % % % %         system(cad);
        % % % %     end
% % % % % %     end
end
return

