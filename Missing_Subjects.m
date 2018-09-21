function Idf = Missing_Subjects(BrainVISADataBase,IdFile);
%
% Syntax :
%   Idf = Missing_Subjects(BrainVISADataBase,IdFile);
%
% This script detects subjects with uncompleted BrainVISA processing.
%
%
% Input Parameters:
%       BrainVISADataBase     : BrainVISA database
%
%       IdFile                : Ids File
%
% Output Parameters:
%      Idf                    : Subjects with uncomplete results
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% July 22th 2014
% Version $1.0

if ~exist('BrainVISADataBase','var')
    error('Wrong BrainVISA Database folder');
    return;
else
    if ~exist([BrainVISADataBase filesep 'subjects'],'dir')
        error('Wrong BrainVISA Database folder');
        return;
    end
end

if ~exist('IdFile','var')
    error('Wrong Ids File ');
    return;
end

% BrainVISADataBase  =  '/media/COSAS/8-BrainVISADataBase-HCP';
% IdFile = '/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';

if exist(deblank(IdFile),'file');
    Ids = char(textread(IdFile,'%s'));
else
    Ids = IdFile;
end
Ns = size(Ids, 1);
Idf = '';
cont = 0;
for i = 1:Ns
    opts.subjid = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' opts.subjid '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    LArgfile = [BrainVISADataBase filesep 'subjects' filesep opts.subjid filesep 't1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/L' opts.subjid '_default_session_auto.arg'];
    RArgfile = [BrainVISADataBase filesep 'subjects' filesep opts.subjid filesep 't1mri/default_acquisition/default_analysis/folds/3.1/default_session_auto/R' opts.subjid '_default_session_auto.arg'];
    if ~exist(LArgfile,'file')|~exist(RArgfile,'file')
        Idf = strvcat(Idf,opts.subjid);
    end
end
return;