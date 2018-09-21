function My_Sulcal_Length_and_Depth_From_Brainvisa(BrainVisaDatabaseDir,IdFile);

FreeSurferDatabaseDir = '/media/HCPData/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/COSAS/8-BrainVISADataBase-HCP';

%IdFile = '/media/MyDisk/PROCESSING_RESULTS/Total_Ids.txt';
IdFile ='/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';

Ids = char(textread(IdFile,'%s'));
Files2Delete = '';

Ns = size(Ids, 1);
for i = 1:Ns
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    %subjId = '1Test_HCP_899885-20140807-T1wMPR1';
    tic;
    
    % ---- Creating Output Directories
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII']);
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcalspams']);
    OutdirSPAM = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcalspams'];
    
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcallength']);
    OutdirLENGTH = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcallength'];
    
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcaldepth']);
    OutdirDEPTH = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcaldepth'];
    
    %% %%%%%%%%%%%%%%%%%%%%%%%  Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Reading ArgFile
    LArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.arg' ];
    
    % ---- Reading Tmtktri File
    LSulcTmtktri = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'aims_Tmtktri.gii' ];
    LSulcTmtktriFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'aims_Tmtktri.mesh' ];
    
    % Converting Sulci Surface to Mesh format
    cad = ['AimsFileConvert -i ' LSulcTmtktri ' -f GIS -o ' LSulcTmtktriFile ' --verbose 0'];
    system(cad);
    
    % Loading Mesh Surface
    Surf = load_mesh(LSulcTmtktriFile);
    contdelsurf = 0;
    for j = 1:length(Surf)
        if isempty(Surf(j).Tri)
            contdelsurf = contdelsurf + 1;
            indel(contdelsurf) = j;
        end
    end
    
    % Removing empty Surfaces
    Surf(indel) = [];
    
    % Saving the new mesh surface
    save_mesh(Surf, LSulcTmtktriFile);
    Files2Delete = strvcat(Files2Delete,LSulcTmtktriFile);
    
    % Reading Sulci attributes form Arg file
    [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(LArgFile);
    NodeNpoint = NodeNpoint(1:length(TmtktriIds));
    
    
    
    % Unifying sulcal labels
    SulcStr = unique(SulcLabels,'rows');
    Ns = size(SulcStr,1);
    
    % All labels in the same line
    SulcNames = '';
    for j = 1:Ns
        SulcNames = [SulcNames '''' deblank(SulcStr(j,:)) '''' ' '];
    end
    SulcNames(end) = [];
    
    %%  Preparing for Sulcal Length and Depth
    
    BottonFullPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'bottom.ima'];
    cmd = [ 'siGraph2Label -g ' LArgFile ' -o ' BottonFullPath ' -l ' SulcNames ' -b aims_bottom' ];
    system(cmd)
    
    HullJunctPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'hull_junction.ima'];
    cmd = [ 'siGraph2Label -g ' LArgFile ' -o ' HullJunctPath ' -l ' SulcNames ' -b aims_junction -s hull_junction' ];
    system(cmd)
    
    
    %%  End of Preparing for Sulcal Length and Depth
    
    
    
    
    
%     SulcName = ['''F.C.M.post._left''' ' ' '''F.C.L.r.ant._left'''];
%     subjId = '1Test_HCP_899885-20140807-T1wMPR1';
    
    
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcallength']);
    OutdirLENGTH = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcallength'];
    
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcaldepth']);
    OutdirDEPTH = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcaldepth'];
    
    
    %1Test_HCP_899885-20140807-T1wMPR1_L_bankssts
    
   
%      cmd = [ 'siGraph2Label -g ' LArgFile ' -o ' HullJunctPath ' -l ' SulcNames ' -b aims_junction -s hull_junction' ];
%      system(cmd)
%      cmd = [ 'siGraph2Label -g ' LArgFile ' -o ' BottonFullPath ' -l ' SulcNames ' -b aims_bottom' ];
%      system(cmd)
    
    for j = 1:size(SulcNames,1)
        sname = deblank(SulcNames(j,:));
        cmd = ['RiiSulcalLength --tvol '  HullJunctPath ' --bvol ' BottonFullPath ' --sm ' LSulcTmtktriFile	' --sn '  sname ...
            ' --sg ' LArgFile ...
            ' -o ' OutdirLENGTH filesep subjId '_L_' sname ...
            ' -s  ' subjId ...
            ' --bv'];
        system(cmd);
         Basename = [OutdirLENGTH filesep subjId '_L_' sname];
        cmd = ['RicSulcalDepth --sm ' Basename '.mesh ' ...
            ' --tm ' Basename '_top.mesh' ...
            ' --bm ' Basename '_bottom.mesh' ...
            ' -o '   OutdirDEPTH filesep subjId '_L_' sname ...
            ' --ns  40' ...
            ' --nd 3 ' ...
            ' --fs 0.2' ...
            ' -s '   subjId ...
            ' --bv '];
        system(cmd);
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        % ----- Sulcal Length
        TempBaseName = [OutdirLENGTH filesep subjId '_L_' SulcName];
        cad = ['RiiSulcalLength --tvol '  HullJunctPath ' --bvol ' BottonFullPath ' --sm ' LSulcTmtktriFile	' --sn '  SulcName ...
            ' --sg ' LArgFile ...
            ' -o ' TempBaseName...
            ' -s  ' subjId ...
            ' --bv'];
        system(cad);
        
        % ----- Sulcal Depth
        cad = ['RicSulcalDepth --sm ' TempBaseName '.mesh ' ...
            ' --tm ' TempBaseName '_top.mesh' ...
            ' --bm ' TempBaseName '_bottom.mesh' ...
            ' -o '   TempBaseName ...
            ' --ns  40' ...
            ' --nd 3 ' ...
            ' --fs 0.2' ...
            ' -s '   subjId ...
            ' --bv '];
        system(cad);
        
        
        
    end
    %matlabpool close;
    
    toc;
    %% %%%%%%%%%%%%%%%%%%%%%  End of Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%
    
    
end
return;


function [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(ArgFile);
%
% Syntax :
% [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(ArgFile);
%
% This function replace the sulcis names, contained in Sulclist, with a new
% sulc name in the Arg File (Brainvisa Sulci organization format)
%
% Input Parameters:
%   InArg             : Input Arg File
%   Sulclist          : List names for sulcis that will be renamed.
%   Nsulcname         : New Sulcis name
%   Outdir            : Output Directory.
%
% Output Parameters:
%   OutArg            : New/Output Arg file with the sulcis names changed
%
%
% Related references:
%
%
% See also: Multi_Replace_Sulc_name_Argfile
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 20th 2012
% Version $1.0


%%
[pth, nm, ext] = fileparts(ArgFile);

%=====================Checking Input Parameters===========================%
% if nargin<4
%     Outdir = pth;
% end
%=========================================================================%
%==========================  Reading Arg File  ===========================%
fio = fopen(ArgFile);
contid = 0;
contpo = 0;
conttm = 0;
cont = 0;
SulcNames = '';
SulcLabels = '';
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    % Finding Sulcus Node Id
    if ~isempty(strfind(lower(line),'*begin node fold'))
        contid = contid + 1;
        temp = strread(line,'%s');
        NodeIds(contid) = str2num(temp{4});
    end
    % Finding Sulcus name
    if ~isempty(strfind(lower(line),'name'))
        if strcmp(lower(line(1:4)),'name')
            temp = strread(line,'%s');
            sname = temp{2};
            SulcNames = strvcat(SulcNames,sname);
        end
    end
    % Finding Sulcus number of points
    if ~isempty(strfind(lower(line),'point_number'))
        if strcmp(lower(line(1:12)),'point_number')
            contpo = contpo + 1;
            temp = strread(line,'%s');
            NodeNpoint(contpo) = str2num(temp{2});
        end
    end
    % Finding Sulcus label
    if ~isempty(strfind(lower(line),'label'))
        if strcmp(lower(line(1:5)),'label')
            temp = strread(line,'%s');
            slabel = temp{2};
            SulcLabels = strvcat(SulcLabels,slabel);
        end
    end
    % Finding Sulcus Tmtktri Id
    if ~isempty(strfind(lower(line),'tmtktri_label'))
        if strcmp(lower(line(1:13)),'tmtktri_label')
            conttm = conttm + 1;
            temp = strread(line,'%s');
            TmtktriIds(conttm) = str2num(temp{2});
        end
    end
    
end
fclose(fio);
return

return
