function  OutTextFile = Import_into_HIPHOP_Parcellation(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% Syntax :
%    OutTextFile = Import_into_HIPHOP_Parcellation(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% This script generates the inputs for HIP HOP Cortical Parcellation pipeline
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
OutTextFile = [BrainVisaDatabaseDir filesep 'Gyral_HIPHOP_Parcellation_Config_File.txt'];


Ids = char(textread(IdFile,'%s'));

% Ids = strvcat('HCP_127630-20140807-T1wMPR2',...
% 'HCP_135225-20140807-T1wMPR1');


Ns = size(Ids,1);

cadt1left= '';
cadt1right = '';
for i = 1:Ns
    Id = deblank(Ids(i,:));
    disp(['Preparing Subject: ' Id '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    
    
    %% Inputs
    
  
%   T1_input

    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' Id '_default_session_auto.arg' '''' ' '];
    cadt1left= [cadt1left ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' Id '_default_session_auto.arg' '''' ' '];
    cadt1right= [cadt1right ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/HCP_100307-20140807-T1wMPR1.nii.gz
%     
end
disp('  ');
disp('Left Automatic Recognition Arg Filenames');
cadt1left(end) = [];
disp(cadt1left);

disp('  ');
disp('Right Automatic Recognition Arg Filenames');
cadt1right(end) = [];
disp(cadt1right);

cads = strvcat('## =========== Inputs ===============','# Left Automatic Recognition Arg Filenames', cadt1left, ' ' , '# Right Automatic Recognition Arg Filenames' ,cadt1right);
fid = fopen(OutTextFile,'wt');
for i = 1:size(cads,1)
    fprintf(fid,'%s\n',deblank(cads(i,:)));
end
fclose(fid);


return;