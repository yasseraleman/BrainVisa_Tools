function Rename_WmSkeletons(BrainVisaDatabaseDir, IdFile);
%
% Syntax :
%    OutTextFile = Import_into_WhiteMatter_Skeleton_RII(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% This script renames the skeleton meshes produced during the Gyral
% Skeleton Calculation( Brainvisa). All the skeleton meshes have
% NRIC_Skel_skel.mesh as name. This occurs only using the iterate option.
%
%
% Input Parameters:
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

BrainVisaDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase';
%BrainVisaDatabaseDir = '/media/MyDisk/Test';
IdFile = '/media/MyDisk/PROCESSING_RESULTS/Total_Ids.txt';
OutTextFile = [BrainVisaDatabaseDir filesep 'WhiteMatter_Skeleton_RII_Config_File.txt'];


Ids = char(textread(IdFile,'%s'));
%Ids = Ids(1:2,:);
% Ids = strvcat('HCP_127630-20140807-T1wMPR2',...
% 'HCP_135225-20140807-T1wMPR1');


Ns = size(Ids,1);
cadt1in = '';
cadt1out= '';
cadtId = '';
for i = 1:Ns
    Id = deblank(Ids(i,:));
    disp(['Preparing Subject: ' Id '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    if exist([BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'spams'],'dir')
        movefile([BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'spams'],[BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'wmspams']);
    end
        
     
end
return;