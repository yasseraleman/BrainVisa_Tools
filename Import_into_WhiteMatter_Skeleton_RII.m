function  OutTextFile = Import_into_WhiteMatter_Skeleton_RII(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% Syntax :
%    OutTextFile = Import_into_WhiteMatter_Skeleton_RII(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% This script generates the inputs for WhiteMatter_Skeleton_RII
% (UTHSCSA/Gyral Skeleton Calculation).
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
BrainVisaDatabaseDir = '/media/MyDisk/Test';
IdFile = '/media/MyDisk/PROCESSING_RESULTS/Total_Ids.txt';
OutTextFile = [BrainVisaDatabaseDir filesep 'WhiteMatter_Skeleton_RII_Config_File.txt'];


Ids = char(textread(IdFile,'%s'));
Ids = Ids(1:2,:);
% Ids = strvcat('HCP_127630-20140807-T1wMPR2',...
% 'HCP_135225-20140807-T1wMPR1');


Ns = size(Ids,1);
cadt1in = '';
cadt1out= '';
cadtId = '';
for i = 1:Ns
    Id = deblank(Ids(i,:));
    disp(['Preparing Subject: ' Id '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'RII']);
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh']);
   % system(['rm -r ' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh' filesep '*']);
    
    %% Inputs
%   T1_input

    % Grey white Input Filenames for Both Hemispheres
%     ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Lgrey_white_' Id '.ima' '''' ' '];
%     cadt1in= [cadt1in ncad];
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Rgrey_white_' Id '.ima' '''' ' '];
    cadt1in= [cadt1in ncad];
    
    
    %% Outputs
    % Output Names
%     ncad = ['''' Id '''' ' '];
%     cadtId= [cadtId ncad];
    ncad = ['''' Id '''' ' '];
    cadtId= [cadtId ncad];
    
    
     % Output Directories
%     ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh' '''' ' '];
%     cadt1out= [cadt1out ncad];
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh' '''' ' '];
    cadt1out= [cadt1out ncad];
%     
end
cadt1in(end) = [];
cadt1out(end) = [];
cadtId(end) = [];

cads = strvcat('## =========== Inputs ===============','# Grey/White Input Images (GIS Format)', cadt1in, ' ' ,'# Subjects IDs', cadtId ,' ' ,' ' , '## =========== Outputs ===============', '# Grey/White Output Images' ,cadt1out);
fid = fopen(OutTextFile,'wt');
for i = 1:size(cads,1)
    fprintf(fid,'%s\n',deblank(cads(i,:)));
end
fclose(fid);

return;