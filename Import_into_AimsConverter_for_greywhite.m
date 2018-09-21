function  OutTextFile = Import_into_AimsConverter_for_greywhite(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% Syntax :
%    OutTextFile = Import_into_AimsConverter(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% This script generates the inputs for AimsConverter (Brainvisa/Tools)
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

FreeSurferDatabaseDir = '/home/yaleman/FreeSurferDataBase';
BrainVisaDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase/BrainVisaDataBase';
IdFile = '/home/yaleman/BrainVisaDataBase/Ids2Process.txt';
OutTextFile = [BrainVisaDatabaseDir filesep 'Grey_White_AimsConverter_Config_File.txt'];


Ids = char(textread(IdFile,'%s'));

% Ids = strvcat('HCP_127630-20140807-T1wMPR2',...
% 'HCP_135225-20140807-T1wMPR1');


Ns = size(Ids,1);
cadt1in = '';
cadt1out= '';
cadtId = '';
for i = 1:Ns
    Id = deblank(Ids(i,:));
    disp(['Preparing Subject: ' Id '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    
    
    %% Inputs
%   T1_input
    % Grey white Input Filenames for Both Hemispheres
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Lgrey_white_' Id '.nii.gz' '''' ' '];
    cadt1in= [cadt1in ncad];
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Rgrey_white_' Id '.nii.gz' '''' ' '];
    cadt1in= [cadt1in ncad];
    
    
    % Grey white Output Filenames for Both Hemispheres
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Lgrey_white_' Id '.ima' '''' ' '];
    cadt1out= [cadt1out ncad];
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Rgrey_white_' Id '.ima' '''' ' '];
    cadt1out= [cadt1out ncad];
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/HCP_100307-20140807-T1wMPR1.nii.gz
%     
end
cadt1in(end) = [];
cadt1out(end) = [];

cads = strvcat('## =========== Inputs ===============','# Grey/White Input Images', cadt1in, ' ' ,' ' , '## =========== Outputs ===============', '# Grey/White Output Images' ,cadt1out);
fid = fopen(OutTextFile,'wt');
for i = 1:size(cads,1)
    fprintf(fid,'%s\n',deblank(cads(i,:)));
end
fclose(fid);

return;