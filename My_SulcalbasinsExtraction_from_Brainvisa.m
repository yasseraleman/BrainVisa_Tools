function My_SulcalbasinsExtraction_from_Brainvisa(BrainVisaDatabaseDir,IdFile);

FreeSurferDatabaseDir = '/media/HCPData/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/COSAS/8-BrainVISADataBase-HCP';
IdFile = '/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';

Ids = char(textread(IdFile,'%s'));
% Ids = '1Test_HCP_899885-20140807-T1wMPR1';

Files2Delete = '';

Nsubj = size(Ids, 1);
for i = 2:Nsubj
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Nsubj)]);
    
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

    %% %%%%%%%%%%%%%%%%%%%%%%%  Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Reading ArgFile
    ArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.arg' ];
    
% %     cad = ['cp -r ' BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data ' BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto_lobar.data' ];
% %     system(cad);
    
    % ---- Reading Tmtktri File
    LSulcTmtktri = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep  'aims_Tmtktri.gii' ];
    LSulcTmtktriFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'aims_Tmtktri.mesh' ];
    
    % Converting Sulci Surface to Mesh format
    cad = ['AimsFileConvert -i ' LSulcTmtktri ' -f GIS -o ' LSulcTmtktriFile ' --verbose 0'];
    system(cad);
    
    % Loading Mesh Surface
    Surf = load_mesh(LSulcTmtktriFile);
    contdelsurf = 0;
    for j = 1:length(Surf)
        if isempty(Surf(j).Tri)
            contdelsurf = contdelsurf + 1;
            indel(contdelsurf) = j;
        end
    end
    
    % Removing empty Surfaces
    Surf(indel) = [];
    
    % Saving the new mesh surface
    save_mesh(Surf, LSulcTmtktriFile);
    Files2Delete = strvcat(Files2Delete,LSulcTmtktriFile);
    
    % Reading Sulci attributes form Arg file
    [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(ArgFile);
    NodeNpoint = NodeNpoint(1:length(TmtktriIds));
    
    
    
    % Unifying sulcal labels
    SulcStr = unique(SulcLabels,'rows');
    Ns = size(SulcStr,1);
    
    % All labels in the same line
    SulcNames = '';
    for j = 1:Ns
        SulcNames = [SulcNames '''' deblank(SulcStr(j,:)) '''' ' '];
    end
    SulcNames(end) = [];
    
     %%  Preparing for Sulcal Length and Depth
    
    BottonFullPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'bottom.ima'];
    cmd = [ 'siGraph2Label -g ' ArgFile ' -o ' BottonFullPath ' -l ' SulcNames ' -b aims_bottom' ];
    system(cmd)
    
    HullJunctPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'hull_junction.ima'];
    cmd = [ 'siGraph2Label -g ' ArgFile ' -o ' HullJunctPath ' -l ' SulcNames ' -b aims_junction -s hull_junction' ];
    system(cmd)
    
    
    
    % Sulcal Arg File
    
    graphbasFile = [Outdir filesep subjId '_LvolumeBasins.txt'];
    volbasFile = [Outdir filesep subjId '_LvolumeBasins.nii.gz'];
    
    cad = [ 'siGraph2Label -g ' ArgFile ...
        ' -a label' ...
        ' -tv ' nobiasFile ...
        ' -tr ' translFile ...
        ' -o ' volbasFile ...
        ' -b aims_junction -b aims_bottom -b aims_ss -b aims_other' ...
        ' -ot ' graphbasFile];
    system(cad);
    
  
    %% %%%%%%%%%%%%%%%%%%%%%%%  End of Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%
    
    %% %%%%%%%%%%%%%%%%%%%%%%%  Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
       % ---- Reading ArgFile
    ArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto_lobar.arg' ];
    
    cad = ['cp -r ' BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data ' BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto_lobar.data' ];
    system(cad);
    
    % ---- Reading Tmtktri File
    LSulcTmtktri = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto_lobar.data' filesep  'aims_Tmtktri.gii' ];
    LSulcTmtktriFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto_lobar.data' filesep 'aims_Tmtktri.mesh' ];
    
    % Converting Sulci Surface to Mesh format
    cad = ['AimsFileConvert -i ' LSulcTmtktri ' -f GIS -o ' LSulcTmtktriFile ' --verbose 0'];
    system(cad);
    
    % Loading Mesh Surface
    Surf = load_mesh(LSulcTmtktriFile);
    contdelsurf = 0;
    for j = 1:length(Surf)
        if isempty(Surf(j).Tri)
            contdelsurf = contdelsurf + 1;
            indel(contdelsurf) = j;
        end
    end
    
    % Removing empty Surfaces
    Surf(indel) = [];
    
    % Saving the new mesh surface
    save_mesh(Surf, LSulcTmtktriFile);
    Files2Delete = strvcat(Files2Delete,LSulcTmtktriFile);
    
    % Reading Sulci attributes form Arg file
    [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(ArgFile);
    NodeNpoint = NodeNpoint(1:length(TmtktriIds));
    
    
    
    % Unifying sulcal labels
    SulcStr = unique(SulcLabels,'rows');
    Ns = size(SulcStr,1);
    
    % All labels in the same line
    SulcNames = '';
    for j = 1:Ns
        SulcNames = [SulcNames '''' deblank(SulcStr(j,:)) '''' ' '];
    end
    SulcNames(end) = [];
    
     %%  Preparing for Sulcal Length and Depth
    
    BottonFullPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto_lobar.data' filesep 'bottom.ima'];
    cmd = [ 'siGraph2Label -g ' ArgFile ' -o ' BottonFullPath ' -l ' SulcNames ' -b aims_bottom' ];
    system(cmd)
    
    HullJunctPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto_lobar.data' filesep 'hull_junction.ima'];
    cmd = [ 'siGraph2Label -g ' ArgFile ' -o ' HullJunctPath ' -l ' SulcNames ' -b aims_junction -s hull_junction' ];
    system(cmd)
    
    
    
    % Sulcal Arg File
    
    graphbasFile = [Outdir filesep subjId '_RlobarvolumeBasins.txt'];
    volbasFile = [Outdir filesep subjId '_RlobarvolumeBasins.nii.gz'];
    
    cad = [ 'siGraph2Label -g ' ArgFile ...
        ' -a label' ...
        ' -tv ' nobiasFile ...
        ' -tr ' translFile ...
        ' -o ' volbasFile ...
        ' -b aims_junction -b aims_bottom -b aims_ss -b aims_other' ...
        ' -ot ' graphbasFile];
    system(cad);
    %% %%%%%%%%%%%%%%%%%%%%%%%  End of Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%
    
    
    
end
return;