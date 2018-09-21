function varargout = My_SulcallinesExtraction_from_Brainvisa(varargin);
%
% Syntax :
%    OutSulcLines = My_SulcallinesExtraction_from_Brainvisa(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% This script computes sulcal lines along cortical white matter surface
% using the curvature map created by FreeSurfer.
%
% Input Parameters:
%     FreeSurferDatabaseDir          : FreeSurfer Database directory.
%     BrainVisaDatabaseDir           : Brainvisa Database Directory.
%     IdFile                         : Ids File.
%
% Output Parameters:
%    OutSulcLines                    : Sulcal Lines Files
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



FreeSurferDatabaseDir = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/Data/PROCESSING_RESULTS/HCP/8-BrainVisaDataBase';
IdFile = '/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';

IdFile = '1Test_HCP_899885-20140807-T1wMPR1';

% Locating Labels Translation Map
[A,B]= system('locate sulci_model_2008.trl');
ind = strfind(B,'.trl');
translFile = deblank(B(1:ind(1)+3));

% Locating File Correspondance Constraint
[A,B]= system('locate constraint_correspondance_2012.txt');
ind = strfind(B,'.txt');
correspoFile = deblank(B(1:ind(1)+3));


% %% ================== Checking Input parameters ========================= %
% if nargin <3
%     error('Please enter a correct number of inputs');
%     return;
% end
% FreeSurferDatabaseDir = varargin{1};
% BrainVisaDatabaseDir = varargin{2};
% IdFile = varargin{3};
% 
% %% ================== End of Checking Input parameters ================== %

if exist(IdFile, 'file')
    Ids = char(textread(IdFile,'%s'));
else
    Ids = IdFile;
end
OutSulcLines = '';
Files2Delete = '';

%% ============================ Main Program ============================ %
Ns = size(Ids, 1);
for i = 1:Ns
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Ns)]);

    % ---- Creating Output Directories 
    % Outdir = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'SulcalLines'];
    
    Outdir = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 'surface'];
    if ~exist(Outdir,'dir')
        mkdir(Outdir);
    end
    
    % Moving FreeSurfer Surfaces into GIS surfaces
    OutSurfFiles = FreeS_Surface_to_Brainvisa_Mesh(FreeSurferDatabaseDir, BrainVisaDatabaseDir, subjId, 'white');
        
    % General Mandatory Inputs
    nobiasFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'nobias_' subjId '.nii.gz'];
   
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
    Files2Delete = strvcat(Files2Delete, curvFile); % Deleting Temporary files
    
    scurvFile = [Outdir filesep subjId '_LsmoothedFScurv.tex'];
    whiteSurf = deblank(OutSurfFiles(1,:));
    
    % Texture Smoothing
    cad = [ 'AimsTextureSmoothing -i '  curvFile ' -o ' scurvFile ' -m ' whiteSurf ' -s 2  -t 0.01'];
    system(cad);
    Files2Delete = strvcat(Files2Delete, scurvFile); % Deleting Temporary files
    
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
    Files2Delete = strvcat(Files2Delete, volbasFile); % Deleting Temporary files
    
    
    % Sulcal Lines File
    sulclinesFile = [Outdir filesep subjId '_Lwhite_sulcalines.tex'];
    
    % Sulcal Lines Computation
    cad  = [ 'AimsSulcalLines' ...
        ' -d '  depthFile ...
        ' -c '  scurvFile ...
        ' -i ' whiteSurf ...
        ' -b ' volbasFile ...
        ' -lb ' graphbasFile ...
        ' -ls ' correspoFile ...
        ' -m 4 -t 2 -st 15'...
        ' -o ' sulclinesFile ...
        ' -si left'...
        ' -sb 50.0'...
        ' -cv 1'];
    system(cad);
    
    OutSulcLines = strvcat(OutSulcLines,sulclinesFile);
    
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
    Files2Delete = strvcat(Files2Delete, curvFile); % Deleting Temporary files
    
    scurvFile = [Outdir filesep subjId '_RsmoothedFScurv.tex'];
    whiteSurf = deblank(OutSurfFiles(2,:));
    
    % Texture Smoothing
    cad = [ 'AimsTextureSmoothing -i '  curvFile ' -o ' scurvFile ' -m ' whiteSurf ' -s 2  -t 0.01'];
    system(cad);
    Files2Delete = strvcat(Files2Delete, curvFile); % Deleting Temporary files
     
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
    Files2Delete = strvcat(Files2Delete, volbasFile); % Deleting Temporary files
    
    % Sulcal Lines File
    sulclinesFile = [Outdir filesep subjId '_Rwhite_sulcalines.tex'];
    
    % Sulcal Lines Computation
    cad  = [ 'AimsSulcalLines' ...
        ' -d '  depthFile ...
        ' -c '  scurvFile ...
        ' -i ' whiteSurf ...
        ' -b ' volbasFile ...
        ' -lb ' graphbasFile ...
        ' -ls ' correspoFile ...
        ' -m 4 -t 2 -st 15'...
        ' -o ' sulclinesFile ...
        ' -si right'...
        ' -sb 50.0'...
        ' -cv 1'];
    system(cad);
    OutSulcLines = strvcat(OutSulcLines,sulclinesFile);
    %% %%%%%%%%%%%%%%%%%%%%%%%  End of Right Hemisphere %%%%%%%%%%%%%%%%%%%
    
    %% ================== Removing Temporary Files ====================== %
    for j = 1:size(Files2Delete)
        system(['rm ' deblank(Files2Delete(j,:))]);
    end
    %% ================ End of Removing Temporary Files ================= %
end

varargout{1} = OutSulcLines;

return;