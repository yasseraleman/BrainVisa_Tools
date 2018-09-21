function Copy_FreeSurfer_toImport_into_BrainVisa(FreesDir, FreesDatabaseDir, IdFile);
%
% Syntax :
%    Import_FreeSurfer_into_BrainVisa(FreesDir, FreesDatabaseDir, IdFile);
%
% This script replicates the freesurfer structure for BrainVisa Importation. 
%
% Input Parameters:
%     FreesDir            : Original FreeSurfer Directory.
%     FreesDatabaseDir    : Database FreeSurfer Directory.
%     IdFile              :Ids File.
%
% Output Parameters:
%
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

Ids = char(textread(IdFile,'%s'));
Ns = size(Ids,1);
for i = 1:Ns
    Id = deblank(Ids(i,:));
    disp(['Copying Subject: ' Id '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    mkdir([FreesDatabaseDir filesep Id filesep 'bem']);
    mkdir([FreesDatabaseDir filesep Id filesep 'label']);
    mkdir([FreesDatabaseDir filesep Id filesep 'mri']);
    mkdir([FreesDatabaseDir filesep Id filesep 'mri' filesep 'orig']);
    mkdir([FreesDatabaseDir filesep Id filesep 'mri' filesep 'transforms']);
    mkdir([FreesDatabaseDir filesep Id filesep 'scripts']);
    mkdir([FreesDatabaseDir filesep Id filesep 'src']);
    mkdir([FreesDatabaseDir filesep Id filesep 'stats']);
    mkdir([FreesDatabaseDir filesep Id filesep 'surf']);
    mkdir([FreesDatabaseDir filesep Id filesep 'tmp']);
    mkdir([FreesDatabaseDir filesep Id filesep 'touch']);
    mkdir([FreesDatabaseDir filesep Id filesep 'trash']);
    mkdir([FreesDatabaseDir filesep Id filesep 'mri']);
    try
        copyfile([FreesDir filesep Id filesep 'mri' filesep 'orig.mgz'], [FreesDatabaseDir filesep Id filesep 'mri' filesep 'orig.mgz']);
    end
    try
        copyfile([FreesDir filesep Id filesep 'mri' filesep 'ribbon.mgz'], [FreesDatabaseDir filesep Id filesep 'mri' filesep 'ribbon.mgz']);
    end
    try
        copyfile([FreesDir filesep Id filesep 'mri' filesep 'transforms' filesep 'talairach.auto.xfm'],[FreesDatabaseDir filesep Id filesep 'mri' filesep 'transforms' filesep 'talairach.auto.xfm']);
    end
end