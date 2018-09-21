function varargout = FreeS_Surface_to_Brainvisa_Mesh(varargin);
%
% Syntax :
%    OutSurfFiles = FreeS_Surface_to_Brainvisa_Mesh(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile, surfType);
%
% This script converts freesurfer surfaces coordinates into GIS coordinate
% system. Brainvisa uses GIS coordinates system.
%
% Input Parameters:
%     FreeSurferDatabaseDir          : FreeSurfer Database directory.
%     BrainVisaDatabaseDir           : Brainvisa Database Directory.
%     IdFile                         : Ids File.
%     surfType                       : Freesurfer Surface type (inflated,
%                                      sphere, pial or white).
%
% Output Parameters:
%    OutSurfFiles                    : FreeSurfer Surfaces in Brainvisa
%                                      Space
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
% FreeSurferDatabaseDir = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing';
% BrainVisaDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase/BrainVisaDataBase';
% IdFile = '/home/yaleman/BrainVisaDataBase/Ids2Process.txt';
% surfType = 'white';

%% ================== Checking Input parameters ========================= %
if nargin <3
    error('Please enter a correct number of inputs');
    return;
end
FreeSurferDatabaseDir = varargin{1};
BrainVisaDatabaseDir = varargin{2};
IdFile = varargin{3};
if nargin == 4
    surfType = varargin{4};
else
    surfType = 'white';
end

%% ================== End of Checking Input parameters ================== %

if exist(IdFile, 'file')
    Ids = char(textread(IdFile,'%s'));
else
    Ids = IdFile;
end

OutSurfFiles = '';
Files2Delete = '';
Ns = size(Ids, 1); % Number of IDs
for i = 1:Ns
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    %subjId = '1Test_HCP_899885-20140807-T1wMPR1';
    tic;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%  Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Reading White Surface File
    LHemiMesh = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep subjId '_L' surfType '.mesh'];
    
    Imfile = [FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'T1.nii'];
    lhsurf = [FreeSurferDatabaseDir filesep subjId filesep 'surf' filesep 'lh.' surfType];
    Taltransf = [FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
    
    % Converting T1 to nifti format
    cad = ['mri_convert -i ' FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep 'T1.mgz -o ' Imfile] ;
    system(cad);
    
    % Reading freesurfer white surface
    [OutFiles, SurfF] = Exp_Surf(lhsurf, '0', '','', 'imp','n');
    Surfwl= SurfF{1};
    Surfwl.Name = 'LH.WHITE';
    
    % Reading Talairach Transformation
    cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
    Surfwl.SurfData.vertices =Surfwl.SurfData.vertices+repmat(cras,[size(Surfwl.SurfData.vertices,1) 1]); % adding RAS center
    Surfwl = freesCS2brainvisaCS(Surfwl,Imfile,'f2b'); % Converting to Brainvisa Coordinate system
    Out_mesh = save_mesh(Surfwl,LHemiMesh); %Saving the mesh
    OutSurfFiles = strvcat(OutSurfFiles,Out_mesh);
    %% %%%%%%%%%%%%%%%%%%%%%  End of Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%  Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Reading White Surface File
    RHemiMesh = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep subjId '_R' surfType '.mesh'];
    
    Imfile = [FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'T1.nii'];
    rhsurf = [FreeSurferDatabaseDir filesep subjId filesep 'surf' filesep 'rh.' surfType];
    
    
    % Reading freesurfer white surface
    [OutFiles, SurfF] = Exp_Surf(rhsurf, '0', '','', 'imp','n');
    Surfwr= SurfF{1};
    Surfwr.Name = 'RH.WHITE';
    
    % Reading Talairach Transformation
    Surfwr.SurfData.vertices =Surfwr.SurfData.vertices+repmat(cras,[size(Surfwr.SurfData.vertices,1) 1]); % adding RAS center
    Surfwr = freesCS2brainvisaCS(Surfwr,Imfile,'f2b'); % Converting to Brainvisa Coordinate system
    Out_mesh = save_mesh(Surfwr,RHemiMesh); %Saving the mesh
    OutSurfFiles = strvcat(OutSurfFiles,Out_mesh);
    %% %%%%%%%%%%%%%%%%%%%%%  End of Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%
    toc;
    delete(Imfile);
end
varargout{1} = OutSurfFiles;
return;
