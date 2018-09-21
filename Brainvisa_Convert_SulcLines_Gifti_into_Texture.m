function OutTextFiles = Brainvisa_Convert_SulcLines_Gifti_into_Texture(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%   
% Syntax :
%    OutTextFiles = Brainvisa_Convert_Gifti_into_Texture(FreeSurferDatabaseDir, BrainVisaDatabaseDir, IdFile);
%
% This script generates the inputs for AimsConverter (Brainvisa/Tools)
%
% Input Parameters:
%     FreeSurferDatabaseDir          : FreeSurfer Database directory.
%     BrainVisaDatabaseDir           : Brainvisa Database Directory.
%     IdFile                         :Ids File.
%
% Output Parameters:
%    OutTextFiles                    : Output Texture File.
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

% % FreeSurferDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing';
% % BrainVisaDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase';
% % IdFile = '/media/MyDisk/PROCESSING_RESULTS/Total_Ids.txt';


Ids = char(textread(IdFile,'%s'));

% Ids = strvcat('HCP_127630-20140807-T1wMPR2',...
% 'HCP_135225-20140807-T1wMPR1');

OutTextFiles = '';
Ns = size(Ids,1);
for i = 1:Ns
    Id = deblank(Ids(i,:));
    disp(['Converting Subject: ' Id '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    
    %% ====================== Left Hemisphere =========================== %
    InputFile = [BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 'surface' filesep Id '_Lwhite_sulcalines.gii' ];
    OutputFile = [BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 'surface' filesep Id '_Lwhite_sulcalines.tex' ];
    cad = ['AimsFileConvert -i ' InputFile  ' -o ' OutputFile];
    if exist(InputFile,'file')
        system(cad);
    else
        OutputFile = '';
    end 
    OutTextFiles = strvcat(OutTextFiles, OutputFile);
    
    %% ====================== Right Hemisphere ========================== %
    InputFile = [BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 'surface' filesep Id '_Rwhite_sulcalines.gii' ];
    OutputFile = [BrainVisaDatabaseDir filesep 'subjects' filesep Id filesep 'surface' filesep Id '_Rwhite_sulcalines.tex' ];
    cad = ['AimsFileConvert -i ' InputFile  ' -o ' OutputFile];
    if exist(InputFile,'file')
        system(cad);
    else
        OutputFile = '';
    end
    OutTextFiles = strvcat(OutTextFiles, OutputFile);
end
return;