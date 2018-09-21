function  OutTextFile = Import_into_Mesh_Watershed(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
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

BrainVisaDatabaseDir = '/media/COSAS/8-BrainVISADataBase-HCP';
IdFile = '/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';
OutTextFile = [BrainVisaDatabaseDir filesep 'Mesh_Watershed_Extraction_Config_File.txt'];


Ids = char(textread(IdFile,'%s'));

% Ids = strvcat('HCP_127630-20140807-T1wMPR2',...
% 'HCP_135225-20140807-T1wMPR1');
%Ids = strvcat('HCP_116524-20140807-T1wMPR1','HCP_125525-20140807-T1wMPR1');

Ns = size(Ids,1);


cadwmeshleft = '';
cadwhitedpfleft = '';
cadwhitepitsleft = '';
cadwhitenoysipitsleft = '';
cadwhiteridgesleft = '';
cadwhitebasinsleft = '';

cadwmeshright = '';
cadwhitedpfright = '';
cadwhitepitsright = '';
cadwhitenoysipitsright = '';
cadwhiteridgesright = '';
cadwhitebasinsright = '';

for i = 2:Ns
    Id = deblank(Ids(i,:));
    disp(['Preparing Subject: ' Id '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    
    
    %% Inputs
    
  
%   T1_input

    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep Id '_Lwhite.gii' '''' ' '];
    cadwmeshleft= [cadwmeshleft ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep 'surface_analysis' filesep  Id '_Lwhite_DPF.gii' '''' ' '];
    cadwhitedpfleft= [cadwhitedpfleft ncad];
    
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep 'surface_analysis' filesep  Id '_Lwhite_pits.gii' '''' ' '];
    cadwhitepitsleft= [cadwhitepitsleft ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep 'surface_analysis' filesep  Id '_Lwhite_noisy_pits.gii' '''' ' '];
    cadwhitenoysipitsleft= [cadwhitenoysipitsleft ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep 'surface_analysis' filesep  Id '_Lwhite_ridges.gii' '''' ' '];
    cadwhiteridgesleft= [cadwhiteridgesleft ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep 'surface_analysis' filesep  Id '_Lwhite_basins.gii' '''' ' '];
    cadwhitebasinsleft= [cadwhitebasinsleft ncad];
    
   

    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep Id '_Rwhite.gii' '''' ' '];
    cadwmeshright= [cadwmeshright ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep 'surface_analysis' filesep  Id '_Rwhite_DPF.gii' '''' ' '];
    cadwhitedpfright= [cadwhitedpfright ncad];
    
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep 'surface_analysis' filesep  Id '_Rwhite_pits.gii' '''' ' '];
    cadwhitepitsright= [cadwhitepitsright ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep 'surface_analysis' filesep  Id '_Rwhite_noisy_pits.gii' '''' ' '];
    cadwhitenoysipitsright= [cadwhitenoysipitsright ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep 'surface_analysis' filesep  Id '_Rwhite_ridges.gii' '''' ' '];
    cadwhiteridgesright= [cadwhiteridgesright ncad];
    
    ncad = ['''' BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep 'surface_analysis' filesep  Id '_Rwhite_basins.gii' '''' ' '];
    cadwhitebasinsright= [cadwhitebasinsright ncad];
    
    
    
    
    
   
%     /home/yaleman/BrainVisaDataBase/subjects/HCP_100307-20140807-T1wMPR1/t1mri/default_acquisition/HCP_100307-20140807-T1wMPR1.nii.gz
%     
end
disp('  ');

cads = strvcat('## =========== Inputs Left Hemisphere ===============',' ','# Input Mesh', cadwmeshleft, ' ' ... 
                                                      ,'# DPF Texture', cadwhitedpfleft ,' ' ...
                                                      ,'# Pits Texture', cadwhitepitsleft ,' ' ...
                                                      ,'# Noysi Pits Texture', cadwhitenoysipitsleft ,' ' ...
                                                      ,'# Ridges Texture', cadwhiteridgesleft ,' ' ...
                                                      ,'# Basins Texture' ,cadwhitebasinsleft,' '...
              ,'## =========== Inputs Right Hemisphere ===============',' ','# Input Mesh', cadwmeshright, ' ' ... 
                                                      ,'# DPF Texture', cadwhitedpfright ,' ' ...
                                                      ,'# Pits Texture', cadwhitepitsright ,' ' ...
                                                      ,'# Noysi Pits Texture', cadwhitenoysipitsright ,' ' ...
                                                      ,'# Ridges Texture', cadwhiteridgesright ,' ' ...
                                                      ,'# Basins Texture' ,cadwhitebasinsright);
           
                                                  
fid = fopen(OutTextFile,'wt');                          
for i = 1:size(cads,1)
    fprintf(fid,'%s\n',deblank(cads(i,:)));
end
fclose(fid);


return;