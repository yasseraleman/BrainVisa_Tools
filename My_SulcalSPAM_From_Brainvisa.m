function My_SulcalSPAM_From_Brainvisa(BrainVisaDatabaseDir,IdFile);

FreeSurferDatabaseDir = '/media/HCPData/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/COSAS/8-BrainVISADataBase-HCP';

IdFile ='/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';

Ids = char(textread(IdFile,'%s'));

Files2Delete = '';

Ns = size(Ids, 1);
for i = 1:Ns
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    %subjId = '1Test_HCP_899885-20140807-T1wMPR1';
    tic;
    % ---- Creating Output Directory
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII']);
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcalspams']);
    Outdir = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'sulcalspams'];
    
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
            indel(contdelsurf) = j
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
    
     % ---- Reading White Surface File
    LHemiMeshGii = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep subjId '_Lhemi.gii'];
    LHemiMesh = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep subjId '_Lhemi.mesh'];
    Files2Delete = strvcat(Files2Delete,LHemiMesh);
    
     %% Uncomment just in case that we want to use the Brainvisa Pial
    % Converting Surface to Mesh format
%     cad = ['AimsFileConvert -i ' LHemiMeshGii ' -f GIS -o ' LHemiMesh ' --verbose 0'];
%     system(cad);
    %% End
    %%
    
    %% Uncomment just in case that we want to use the freesurfer Pial
    %% matter surface
    % Files
    Imfile = [FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'T1.nii'];
    lhsurf = [FreeSurferDatabaseDir filesep subjId filesep 'surf' filesep 'lh.pial'];
    Taltransf = [FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
    
    % Converting T1 to nifti format
    cad = ['mri_convert -i ' FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep 'T1.mgz -o ' Imfile] ;
    system(cad);
    
    % Reading freesurfer white surface
    [OutFiles, SurfF] = Exp_Surf(lhsurf, '0', '','', 'imp','n');
    Surfwl= SurfF{1};
    Surfwl.Name = 'LH.WHITE';
    
    % Reading Talairach Transformation
    cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
    Surfwl.SurfData.vertices =Surfwl.SurfData.vertices+repmat(cras,[size(Surfwl.SurfData.vertices,1) 1]); % adding RAS center
    Surfwl = freesCS2brainvisaCS(Surfwl,Imfile,'f2b'); % Converting to Brainvisa Coordinate system
    Out_mesh = save_mesh(Surfwl,LHemiMesh); %Saving the mesh
    
    %% End 
    %%
    
    
    % ---- Computing Gyral White Matter Spams
    %Nstruct = size(SulcNames,1);
    SulcStr = unique(SulcLabels,'rows');
    Ns = size(SulcStr,1);
    %matlabpool open;
    for j = 1:Ns
        SulcName = deblank(SulcStr(j,:));
        ind = find(ismember(SulcLabels,SulcName,'rows')); % Same Sulci
        pointdistr = NodeNpoint(ind);
        cad = ['RicSulcalSpan'...
            ' --sg ' LArgFile ...
            ' -o ' Outdir filesep subjId '_L' ...
            ' -n ' subjId ...
            ' --bv' ...
            ' --sm ' LSulcTmtktriFile ...
            ' --gm ' LHemiMesh ...
            ' --sn ' SulcName ...
            ' --maxa 1' ...
            ' --mina  0.9' ...
            ' --maxd 10' ...
            ' --mind 0.1' ...
            ' --mins ' num2str(floor(mean(pointdistr) - 0.5*mean(pointdistr)))...
            ' --sr 4' ...
            ' --nsd 3'];
        
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
