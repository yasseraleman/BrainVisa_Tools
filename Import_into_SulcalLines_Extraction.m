function  OutTextFile = Import_into_SulcalLines_Extraction(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% Syntax :
%    OutTextFile = Import_into_SulcalLines_Extraction(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% This script generates the inputs for Sulcal Lines Extraction Pipeline.
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

FreeSurferDatabaseDir = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase/BrainVisaDataBase';
IdFile = '/home/yaleman/BrainVisaDataBase/Ids2Process.txt';
OutTextFile = [BrainVisaDatabaseDir filesep 'Sulcal_Lines_Extraction_Config_File.txt'];


Ids = char(textread(IdFile,'%s'));

% Ids = strvcat('HCP_127630-20140807-T1wMPR2',...
% 'HCP_135225-20140807-T1wMPR1');
%Ids = strvcat('HCP_116524-20140807-T1wMPR1','HCP_125525-20140807-T1wMPR1');

Ns = size(Ids,1);

cadt1left= '';
cadwmeshleft = '';
cadgreywhiteleft = '';
cadt1right = '';
cadwmeshright = '';
cadgreywhiteright = '';
cadsulctrl = '';
cadcorresp = '';
for i = 1:Ns
    Id = deblank(Ids(i,:));
    disp(['Preparing Subject: ' Id '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    
    
    %% Inputs
    
  
%   T1_input

    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' Id '_default_session_auto.arg' '''' ' '];
    cadt1left= [cadt1left ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep  Id '_Lwhite.mesh' '''' ' '];
    cadwmeshleft= [cadwmeshleft ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Lgrey_white_' Id '.ima' '''' ' '];
    cadgreywhiteleft= [cadgreywhiteleft ncad];
    
    
    ncad = ['''' '/usr/local/brainvisa/share/brainvisa-share-4.4/nomenclature/translation/sulci_model_2008.trl' '''' ' '];
    cadsulctrl= [cadsulctrl ncad];
    
    ncad = ['''' '/usr/local/brainvisa/share/brainvisa-share-4.4/nomenclature/surfaceanalysis/constraint_correspondance_2012.txt' '''' ' '];
    cadcorresp= [cadcorresp ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' Id '_default_session_auto.arg' '''' ' '];
    cadt1right= [cadt1right ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep  Id '_Rwhite.mesh' '''' ' '];
    cadwmeshright= [cadwmeshright ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'Rgrey_white_' Id '.ima' '''' ' '];
    cadgreywhiteright= [cadgreywhiteright ncad];
    
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

cads = strvcat('## =========== Inputs ===============','# Left Automatic Recognition Arg Filenames', cadt1left, ' ' ... 
                                                      ,'# Left White Mesh', cadwmeshleft ,' ' ...
                                                      ,'# Left Grey-White Image', cadgreywhiteleft ,' ' ...
                                                      ,'# Labels Translation Map', cadsulctrl ,' ' ...
                                                      ,'# Correspondance Constrain', cadcorresp ,' ' ...
                                                      ,'# Right Automatic Recognition Arg Filenames' ,cadt1right,' '...
                                                      ,'# Right White Mesh', cadwmeshright ,' ' ...
                                                      ,'# Right Grey-White Image', cadgreywhiteright, ' '...
                                                      ,'# Labels Translation Map', cadsulctrl,' '...
                                                      ,'# Correspondance Constrain', cadcorresp);
fid = fopen(OutTextFile,'wt');                          
for i = 1:size(cads,1)
    fprintf(fid,'%s\n',deblank(cads(i,:)));
end
fclose(fid);


return;