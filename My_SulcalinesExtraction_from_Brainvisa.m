function My_SulcallinesExtraction_from_Brainvisa(BrainVisaDatabaseDir,IdFile);

FreeSurferDatabaseDir = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase';
IdFile = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing/HCPwmPR1_Ids.txt';

Ids = char(textread(IdFile,'%s'));
% Ids = '1Test_HCP_899885-20140807-T1wMPR1';

Files2Delete = '';

Ns = size(Ids, 1);
for i = 1:Ns
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    
    tic;
    % ---- Creating Output Directories
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 'surface' ]);
    Outdir = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 'surface' ];
    
    % Moving FreeSurfer Surfaces into GIS surfaces
    OutSurfFiles = FreeS_Surface_to_Brainvisa_Mesh(FreeSurferDatabaseDir, BrainVisaDatabaseDir, subjId, 'white');
    
    % General Mandatory Inputs
    nobiasFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'nobias_' subjId '.nii.gz'];
    translFile = ['/usr/local/brainvisa/share/brainvisa-share-4.4/nomenclature/translation/sulci_model_2008.trl'];
    correspoFile = '/usr/local/brainvisa/share/brainvisa-share-4.4/nomenclature/surfaceanalysis/constraint_correspondance_2012.txt';
    
    %% %%%%%%%%%%%%%%%%%%%%%%%  Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Reading Lgrey_white image
    greywhiteFileNii = [ BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Lgrey_white_' subjId '.nii.gz'];
    greywhiteFile = [ BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Lgrey_white_' subjId '.ima'];
    
    if ~exist(greywhiteFile,'file')&exist(greywhiteFileNii,'file')
        cad = ['AimsFileConvert -i ' greywhiteFileNii ' -o ' greywhiteFile];
        system(cad);
    elseif ~exist(greywhiteFile,'file')&~exist(greywhiteFileNii,'file')
        %% Aqui hay que poner lo de convertir el ribbon en LGrey_white
    end
    
    % Reading FreeSurfer Curvature file
    [txt] = read_cfiles([FreeSurferDatabaseDir filesep subjId filesep  'surf' filesep 'lh.curv']);
    Text.Datatype = 'float32';
    Text.Values = -1*txt; % Negative values in sulci (FreeSurfer have the opposite)
    
    % Saving FreeSurfer Curvature
    [curvFile] = save_texBrainvisa(Text,[Outdir filesep subjId '_LFScurv.tex']);
    scurvFile = [Outdir filesep subjId '_LsmoothedFScurv.tex'];
    whiteSurf = deblank(OutSurfFiles(1,:));
    
    % Texture Smoothing
    cad = [ 'AimsTextureSmoothing -i '  curvFile ' -o ' scurvFile ' -m ' whiteSurf ' -s 2  -t 0.01'];
    system(cad);
    
    % Computing Geodesic Depth
    depthFile = [Outdir filesep subjId '_LgeoDepth.tex'];
    cad = [ 'AimsMeshGeodesicDepth -v ' greywhiteFile ' -i ' whiteSurf ' -c 10 -e 8 -o ' depthFile];
    system(cad);
    
    % Sulcal Arg File
    ArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.arg' ];
    graphbasFile = [Outdir filesep subjId '_LgraphLabelBasins.txt'];
    volbasFile = [Outdir filesep subjId '_LvolumeBasins.nii.gz'];
    
    cad = [ 'siGraph2Label -g ' ArgFile ...
        ' -a label' ...
        ' -tv ' nobiasFile ...
        ' -tr ' translFile ...
        ' -o ' volbasFile ...
        ' -b aims_junction -b aims_bottom -b aims_ss -b aims_other' ...
        ' -ot ' graphbasFile];
    system(cad);
    
    % Sulcal Lines File
    sulclinesFile = [Outdir filesep subjId '_Lwhite_sulcalines_GPDM.tex'];
    
    % Sulcal Lines Computation
    cad  = [ 'AimsSulcalLines' ...
        ' -d '  depthFile ...
        ' -c '  scurvFile ...
        ' -i ' whiteSurf ...
        ' -b ' volbasFile ...
        ' -lb ' graphbasFile ...
        ' -ls ' correspoFile ...
        ' -m 4 -t 1 -st 15'...
        ' -o ' sulclinesFile ...
        ' -si left'...
        ' -sb 100.0'...
        ' -cv 1'];
    system(cad);
    try;copyfile('/media/COSAS/scripts/sulcalines_end.tex',sulclinesFile);end
    %% %%%%%%%%%%%%%%%%%%%%%%%  End of Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%
    
    %% %%%%%%%%%%%%%%%%%%%%%%%  Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Reading Lgrey_white image
    greywhiteFileNii = [ BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Rgrey_white_' subjId '.nii.gz'];
    greywhiteFile = [ BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Rgrey_white_' subjId '.ima'];
    
    if ~exist(greywhiteFile,'file')&exist(greywhiteFileNii,'file')
        cad = ['AimsFileConvert -i ' greywhiteFileNii ' -o ' greywhiteFile];
        system(cad);
    elseif ~exist(greywhiteFile,'file')&~exist(greywhiteFileNii,'file')
        %% Aqui hay que poner lo de convertir el ribbon en RGrey_white
    end
    
    % Reading FreeSurfer Curvature file
    [txt] = read_cfiles([FreeSurferDatabaseDir filesep subjId filesep  'surf' filesep 'rh.curv']);
    Text.Datatype = 'float32';
    Text.Values = -1*txt; % Negative values in sulci (FreeSurfer have the opposite)
    
    % Saving FreeSurfer Curvature
    [curvFile] = save_texBrainvisa(Text,[Outdir filesep subjId '_RFScurv.tex']);
    scurvFile = [Outdir filesep subjId '_RsmoothedFScurv.tex'];
    whiteSurf = deblank(OutSurfFiles(2,:));
    
    % Texture Smoothing
    cad = [ 'AimsTextureSmoothing -i '  curvFile ' -o ' scurvFile ' -m ' whiteSurf ' -s 2  -t 0.01'];
    system(cad);
    
    % Computing Geodesic Depth
    depthFile = [Outdir filesep subjId '_RgeoDepth.tex'];
    cad = [ 'AimsMeshGeodesicDepth -v ' greywhiteFile ' -i ' whiteSurf ' -c 10 -e 8 -o ' depthFile];
    system(cad);
    
    % Sulcal Arg File
    ArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.arg' ];
    graphbasFile = [Outdir filesep subjId '_RgraphLabelBasins.txt'];
    volbasFile = [Outdir filesep subjId '_RvolumeBasins.nii.gz'];
    
    cad = [ 'siGraph2Label -g ' ArgFile ...
        ' -a label' ...
        ' -tv ' nobiasFile ...
        ' -tr ' translFile ...
        ' -o ' volbasFile ...
        ' -b aims_junction -b aims_bottom -b aims_ss -b aims_other' ...
        ' -ot ' graphbasFile];
    system(cad);
    
    % Sulcal Lines File
    sulclinesFile = [Outdir filesep subjId '_Rwhite_sulcalines_GPDM.tex'];
    
    % Sulcal Lines Computation
    cad  = [ 'AimsSulcalLines' ...
        ' -d '  depthFile ...
        ' -c '  scurvFile ...
        ' -i ' whiteSurf ...
        ' -b ' volbasFile ...
        ' -lb ' graphbasFile ...
        ' -ls ' correspoFile ...
        ' -m 4 -t 1 -st 15'...
        ' -o ' sulclinesFile ...
        ' -si right'...
        ' -sb 100.0'...
        ' -cv 1'];
    system(cad);
    try;copyfile('/media/COSAS/scripts/sulcalines_end.tex',sulclinesFile);end
    %% %%%%%%%%%%%%%%%%%%%%%%%  End of Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%
    
    
    
end
return;