function  OutTextFile = Import_into_MorphologistPipeline(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% Syntax :
%    OutTextFile = Import_into_MorphologistPipeline(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% This script generates the inputs for morphologist pipeline
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

FreeSurferDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase';
IdFile = '/media/MyDisk/PROCESSING_RESULTS/Total_Ids.txt';
OutTextFile = [BrainVisaDatabaseDir filesep 'Config_Import_File.txt'];


Ids = char(textread(IdFile,'%s'));

% Ids = strvcat('HCP_127630-20140807-T1wMPR2',...
% 'HCP_135225-20140807-T1wMPR1');


Ns = size(Ids,1);

cadt1out= '';

for i = 1:Ns
    Id = deblank(Ids(i,:));
    disp(['Preparing Subject: ' Id '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    
    
    %% Inputs
    
  
%   T1_input
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep Id '.nii.gz' '''' ' '];
    cadt1out= [cadt1out ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/HCP_100307-20140807-T1wMPR1.nii.gz
%     
end
disp('  ');
disp('Input Filenames');
cadt1out(end) = [];
disp(cadt1out);

return;