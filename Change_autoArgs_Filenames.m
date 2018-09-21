function varagout = Change_autoArgs_Filenames(varagin);
%
% Syntax :
%   outArgFilenames = Change_autoArgs_Filenames(brainvisaDir, IdFile, lookForcad, rePlacecad);
%
% This function renames arg files for a specified subjects.
%
% Input Parameters:
%   brainvisaDir       : BrainVisa Database Directory.
%   IdFile             : Subjects Ids.
%   lookForcad         : String of characters that will be replaced.
%   rePlacecad         : String of characters that will be puted in place of the lookForcad.
%
% Output Parameters:
%   outArgFilenames    : Replaced Arg Filenames
%
% Related references:
%
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 30th 2006
% Version $1.0

%% ====================== Checking input parameters ===================== %
if nargin<4 % the indispensable input arguments are not provided
    error('Four inputs are mandatory');
else
    brainvisaDir = varagin{1};
    IdFile = varagin{2};
    lookForcad = varagin{3};
    rePlacecad = varagin{4};
end
% brainvisaDir  =  '/media/COSAS/8-BrainVISADataBase-HCP';
% IdFile = '/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';

%% =================== End of checking input parameters ================= %

%% ======================== Main Program ================================ %
if exist(deblank(IdFile),'file');
    Ids = char(textread(IdFile,'%s'));
else
    Ids = IdFile;
end

Ns = size(Ids, 1);
outArgFilenames = '';
for i = 1:Ns
    opts.subjid = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' opts.subjid '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    argDir = [brainvisaDir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' ];
    tempFilenames = dir([argDir filesep '*' lookForcad '.arg']);
    for j = 1:length(tempFilenames)
        origName = [argDir filesep tempFilenames(j).name];
        
        ind = strfind(tempFilenames(j).name,'default_session');
        rawcad = tempFilenames(j).name(1:ind+19);
        if ~strcmp(rePlacecad,'auto')
            replaceName = [argDir filesep rawcad '_' rePlacecad '.arg'];
        else
            replaceName = [argDir filesep rawcad '.arg'];
        end
        try
            movefile(origName,replaceName);
            outArgFilenames = strvcat(outArgFilenames, replaceName);
        end
    end
end
%% ==================== End of Main Program ============================= %
% -- Outputs
varagout{1} = outArgFilenames;
return;