function My_Sulcal_Processing_From_Brainvisa(BrainVisaDatabaseDir,IdFile);

FreeSurferDatabaseDir = '/media/HCPData/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/COSAS/8-BrainVISADataBase-HCP';
IdFile ='/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';

Ids = char(textread(IdFile,'%s'));

Ns = size(Ids, 1);
%matlabpool('open',3);
for i = 82:101
    Files2Delete = '';
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    %subjId = '1Test_HCP_899885-20140807-T1wMPR1';
    tic;
    
    % ---- Creating Output Directories
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII1']);
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII1' filesep 'sulcalspams']);
    OutdirSPAM = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII1' filesep 'sulcalspams'];
    
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII1' filesep 'sulcallength']);
    OutdirLENGTH = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII1' filesep 'sulcallength'];
    
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII1' filesep 'sulcaldepth']);
    OutdirDEPTH = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII1' filesep 'sulcaldepth'];
    
    % Reading Talairach Transformation
    Taltransf = [FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
    cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
    
    % Converting T1 to nifti format
    Imfile = [FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'T1.nii'];
    cad = ['mri_convert -i ' FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep 'T1.mgz -o ' Imfile] ;
    system(cad);
    
    
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
    % % % %     contdelsurf = 0;
    % % % %     indel = 0;
    % % % %     for j = 1:length(Surf)
    % % % %         if isempty(Surf(j).Tri)
    % % % %             contdelsurf = contdelsurf + 1;
    % % % %             indel(contdelsurf) = j;
    % % % %         else
    % % % %             try
    % % % %                 Surft = Sulcal_Face_Labelling(Surf(j));
    % % % %                 Normals = Surft.SurfData.VertexNormals;
    % % % %                 [Trip] = Vert_Neibp(double(Surft.SurfData.faces),size(Surft.SurfData.vertices,1),size(Surft.SurfData.faces,1));
    % % % %                 Temp = sum(Trip);
    % % % %                 Trip(:,Temp==0) = [];
    % % % %                 temp = Trip(:,3:end);
    % % % %                 indz = find(temp == 0);
    % % % %
    % % % %                 Coord = [repmat(Trip(:,1),[size(temp,2) 1]) temp(:)];
    % % % %                 ind = find(sum(logical(Coord),2) == 2);
    % % % %                 Coord = Coord(ind,:);
    % % % %                 temp = sort(Coord')';
    % % % %                 Coord = unique(temp,'rows');
    % % % %
    % % % %                 T = Surft.Is(Coord);
    % % % %                 ind = find(T(:,1) - T(:,2) ~=0); % Crossing edges between walls
    % % % %                 Coord = Coord(ind,:);
    % % % %
    % % % %                 indd = accumarray(Coord(:),Coord(:)*0+1); % Points that only appears in a single edge
    % % % %                 ind = find(indd == 1);
    % % % %                 edge2rem = find(sum(ismember(Coord,ind),2) ~=0);
    % % % %                 Coord(edge2rem,:) = []; % Removing unconnected edges
    % % % %                 Surf(j).SurfData.vertices(Coord(:),:) = Surf(j).SurfData.vertices(Coord(:),:) + 2*Normals(Coord(:),:);
    % % % %             catch
    % % % %                 contdelsurf = contdelsurf + 1;
    % % % %                 indel(contdelsurf) = j;
    % % % %             end
    % % % %         end
    % % % %     end
    
    % % % %     if sum(indel) > 0
    % % % %         % Removing empty Surfaces
    % % % %         Surf(indel) = [];
    % % % %     end
    
    
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
    
    
    % % %     % Saving the new mesh surface
    % % %     save_mesh(Surf, LSulcTmtktriFile);
    % % %     Files2Delete = strvcat(Files2Delete,LSulcTmtktriFile);
    
    % Reading Sulci attributes form Arg file
    [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktRII1ds] = Read_SulcalArgFiles(LArgFile);
    NodeNpoint = NodeNpoint(1:length(TmtktRII1ds));
    
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
    lhsurf = [FreeSurferDatabaseDir filesep subjId filesep 'surf' filesep 'lh.pial'];
    
    
    % Reading freesurfer white surface
    [OutFiles, SurfF] = Exp_Surf(lhsurf, '0', '','', 'imp','n');
    Surfwl= SurfF{1};
    Surfwl.Name = 'LH.WHITE';
    
    Surfwl.SurfData.vertices =Surfwl.SurfData.vertices+repmat(cras,[size(Surfwl.SurfData.vertices,1) 1]); % adding RAS center
    Surfwl = freesCS2brainvisaCS(Surfwl,Imfile,'f2b'); % Converting to Brainvisa Coordinate system
    Out_mesh = save_mesh(Surfwl,LHemiMesh); %Saving the mesh
    
    %% End
    %%
    
    % Unifying sulcal labels
    SulcStr = unique(SulcLabels,'rows');
    
    %     [SulcStr, colHeaders] = Detecting_BrainVisa_Groups('/media/COSAS/scripts/BrainVisa_myTools/Brainvisa_nomenclature_sulci.txt', 'left');
    %
    [groupsIds, groupsNames, SulcStr, labelsIds] = Detecting_BrainVisa_Groups('/media/COSAS/scripts/BrainVisa_myTools/Brainvisa_nomenclature_sulci+STS.txt', 'left');
    
    Nstruct = size(SulcStr,1);
    
% % % %     % All labels in the same line
% % % %     SulcNames = '';
% % % %     for j = 1:Nstruct
% % % %         SulcNames = [SulcNames '''' deblank(SulcStr(j,:)) '''' ' '];
% % % %     end
% % % %     SulcNames(end) = [];
    
    %%  Preparing for Sulcal Length and Depth
    
% % % %         BottonFullPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'bottom.ima'];
% % % %         cmd = [ 'siGraph2Label -g ' LArgFile ' -o ' BottonFullPath ' -l ' SulcNames ' -b aims_bottom' ];
% % % %         system(cmd)
% % % %     
% % % %         HullJunctPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'hull_junction.ima'];
% % % %         cmd = [ 'siGraph2Label -g ' LArgFile ' -o ' HullJunctPath ' -l ' SulcNames ' -b aims_junction -s hull_junction' ];
% % % %         system(cmd)
    % % % %
    % % % %
    % % % %
    % % % %     HullJuncttemp = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'hull_junction.nii'];
    % % % %     cad = ['AimsFileConvert -i ' HullJunctPath ' -o ' HullJuncttemp];
    % % % %     system(cad);
    % % % %     Vt = spm_vol(HullJuncttemp);
    % % % %     It = spm_read_vols(Vt);
    % % % %     for i = 1:1
    % % % %         It = imdilate(It,strel(ones(3,3,3)));
    % % % %     end
    % % % %     spm_write_vol(Vt,It);
    % % % %     cad = ['AimsFileConvert -i ' Vt.fname ' -o ' HullJunctPath];
    % % % %     system(cad);
    % % % %     %
    % % % %     BottonJuncttemp = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'bottom.nii'];
    % % % %     cad = ['AimsFileConvert -i ' BottonFullPath ' -o ' BottonJuncttemp];
    % % % %     system(cad);
    % % % %     Vt = spm_vol(BottonJuncttemp);
    % % % %     It = spm_read_vols(Vt);
    % % % %     for i = 1:1
    % % % %         It = imdilate(It,strel(ones(3,3,3)));
    % % % %     end
    % % % %     spm_write_vol(Vt,It);
    % % % %     cad = ['AimsFileConvert -i ' Vt.fname ' -o ' BottonFullPath];
    % % % %     system(cad);
    
    
    %%  End of Preparing for Sulcal Length and Depth
    
    
    % ---- Computing Gyral White Matter Spams
    
    %matlabpool open;
    
    
    for j = 1:Nstruct
        SulcName = deblank(SulcStr(j,:));
        
        % ----- Sulcal SPAM
        ind = find(ismember(SulcLabels,SulcName,'rows')); % Same Sulci
        try
            pointdistr = NodeNpoint(ind);
            cad = ['RicSulcalSpan'...
                ' --sg ' LArgFile ...
                ' -o ' OutdirSPAM filesep subjId '_L' ...
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
            
            
            BottonFullPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'bottom.ima'];
            cmd = [ 'siGraph2Label -g ' LArgFile ' -o ' BottonFullPath ' -l ' SulcName ' -b aims_bottom' ];
            system(cmd);
            
            HullJunctPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'hull_junction.ima'];
            cmd = [ 'siGraph2Label -g ' LArgFile ' -o ' HullJunctPath ' -l ' SulcName ' -b aims_junction -s hull_junction' ];
            system(cmd);
            
            
            % %         HullJuncttemp = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'hull_junction.nii'];
            % %         cad = ['AimsFileConvert -i ' HullJunctPath ' -o ' HullJuncttemp];
            % %         system(cad);
            % %         Vt = spm_vol(HullJuncttemp);
            % %         It = spm_read_vols(Vt);
            % %         for i = 1:1
            % %             It = imdilate(It,strel(ones(3,3,3)));
            % %         end
            % %         spm_write_vol(Vt,It);
            % %         cad = ['AimsFileConvert -i ' Vt.fname ' -o ' HullJunctPath];
            % %         system(cad);
            % %         %
            % %         BottonJuncttemp = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'bottom.nii'];
            % %         cad = ['AimsFileConvert -i ' BottonFullPath ' -o ' BottonJuncttemp];
            % %         system(cad);
            % %         Vt = spm_vol(BottonJuncttemp);
            % %         It = spm_read_vols(Vt);
            % %         for i = 1:1
            % %             It = imdilate(It,strel(ones(3,3,3)));
            % %         end
            % %         spm_write_vol(Vt,It);
            % %         cad = ['AimsFileConvert -i ' Vt.fname ' -o ' BottonFullPath];
            % %         system(cad);
            
            
            
            % ----- Sulcal Length
            
            TempBaseNameL = [OutdirLENGTH filesep subjId '_L_' SulcName];
            cad = ['RiiSulcalLength --tvol '  HullJunctPath ' --bvol ' BottonFullPath ' --sm ' LSulcTmtktriFile	' --sn '  SulcName ...
                ' --sg ' LArgFile ...
                ' -o ' TempBaseNameL...
                ' -s  ' subjId ...
                ' --bv'];
            system(cad);
            
            
            TempBaseName =  [OutdirDEPTH filesep subjId '_L_' SulcName];
            % ----- Sulcal Depth
            cad = ['RicSulcalDepth --sm ' TempBaseNameL '.mesh ' ...
                ' --tm ' TempBaseNameL '_top.mesh' ...
                ' --bm ' TempBaseNameL '_bottom.mesh' ...
                ' -o '   TempBaseName ...
                ' --ns  40' ...
                ' --nd 3 ' ...
                ' --fs 0.2' ...
                ' -s '   subjId ...
                ' --bv '];
            system(cad);
        end
    end
    %matlabpool close;
    
    toc;
    %% %%%%%%%%%%%%%%%%%%%%%  End of Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%  Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Reading ArgFile
    RArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.arg' ];
    
    % ---- Reading Tmtktri File
    RSulcTmtktri = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'aims_Tmtktri.gii' ];
    RSulcTmtktriFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'aims_Tmtktri.mesh' ];
    
    % Converting Sulci Surface to Mesh format
    cad = ['AimsFileConvert -i ' RSulcTmtktri ' -f GIS -o ' RSulcTmtktriFile ' --verbose 0'];
    system(cad);
    
    % Loading Mesh Surface
    Surf = load_mesh(RSulcTmtktriFile);
    % % % %     contdelsurf = 0;
    % % % %     indel = 0;
    % % % %     for j = 1:length(Surf)
    % % % %         if isempty(Surf(j).Tri)
    % % % %             contdelsurf = contdelsurf + 1;
    % % % %             indel(contdelsurf) = j;
    % % % %         else
    % % % %             try
    % % % %                 Surft = Sulcal_Face_Labelling(Surf(j));
    % % % %                 Normals = Surft.SurfData.VertexNormals;
    % % % %                 [Trip] = Vert_Neibp(double(Surft.SurfData.faces),size(Surft.SurfData.vertices,1),size(Surft.SurfData.faces,1));
    % % % %                 Temp = sum(Trip);
    % % % %                 Trip(:,Temp==0) = [];
    % % % %                 temp = Trip(:,3:end);
    % % % %                 indz = find(temp == 0);
    % % % %
    % % % %                 Coord = [repmat(Trip(:,1),[size(temp,2) 1]) temp(:)];
    % % % %                 ind = find(sum(logical(Coord),2) == 2);
    % % % %                 Coord = Coord(ind,:);
    % % % %                 temp = sort(Coord')';
    % % % %                 Coord = unique(temp,'rows');
    % % % %
    % % % %                 T = Surft.Is(Coord);
    % % % %                 ind = find(T(:,1) - T(:,2) ~=0); % Crossing edges between walls
    % % % %                 Coord = Coord(ind,:);
    % % % %
    % % % %                 indd = accumarray(Coord(:),Coord(:)*0+1); % Points that only appears in a single edge
    % % % %                 ind = find(indd == 1);
    % % % %                 edge2rem = find(sum(ismember(Coord,ind),2) ~=0);
    % % % %                 Coord(edge2rem,:) = []; % Removing unconnected edges
    % % % %                 Surf(j).SurfData.vertices(Coord(:),:) = Surf(j).SurfData.vertices(Coord(:),:) + 2*Normals(Coord(:),:);
    % % % %             catch
    % % % %                 contdelsurf = contdelsurf + 1;
    % % % %                 indel(contdelsurf) = j;
    % % % %             end
    % % % %         end
    % % % %     end
    % % % %
    % % % %     if sum(indel) > 0
    % % % %         % Removing empty Surfaces
    % % % %         Surf(indel) = [];
    % % % %     end
    % % % %
    % % % %     % Saving the new mesh surface
    % % % %     save_mesh(Surf, RSulcTmtktriFile);
    % % % %     Files2Delete = strvcat(Files2Delete,RSulcTmtktriFile);
    
    
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
    save_mesh(Surf, RSulcTmtktriFile);
    Files2Delete = strvcat(Files2Delete,RSulcTmtktriFile);
    
    
    % Reading Sulci attributes form Arg file
    [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktRII1ds] = Read_SulcalArgFiles(RArgFile);
    NodeNpoint = NodeNpoint(1:length(TmtktRII1ds));
    
    % ---- Reading White Surface File
    RHemiMeshGii = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep subjId '_Rhemi.gii'];
    RHemiMesh = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep subjId '_Rhemi.mesh'];
    Files2Delete = strvcat(Files2Delete,RHemiMesh);
    
    %% Uncomment just in case that we want to use the Brainvisa Pial
    % Converting Surface to Mesh format
    %     cad = ['AimsFileConvert -i ' RHemiMeshGii ' -f GIS -o ' RHemiMesh ' --verbose 0'];
    %     system(cad);
    %% End
    %%
    
    %% Uncomment just in case that we want to use the freesurfer Pial
    %% matter surface
    % Files
    rhsurf = [FreeSurferDatabaseDir filesep subjId filesep 'surf' filesep 'rh.pial'];
    
    
    % Reading freesurfer white surface
    [OutFiles, SurfF] = Exp_Surf(rhsurf, '0', '','', 'imp','n');
    Surfwr= SurfF{1};
    Surfwr.Name = 'RH.WHITE';
    
    % Reading Talairach Transformation
    Surfwr.SurfData.vertices =Surfwr.SurfData.vertices+repmat(cras,[size(Surfwr.SurfData.vertices,1) 1]); % adding RAS center
    Surfwr = freesCS2brainvisaCS(Surfwr,Imfile,'f2b'); % Converting to Brainvisa Coordinate system
    Out_mesh = save_mesh(Surfwr,RHemiMesh); %Saving the mesh
    
    %% End
    %%
    
    % Unifying sulcal labels
    SulcStr = unique(SulcLabels,'rows');
    
    [groupsIds, groupsNames, SulcStr, labelsIds] = Detecting_BrainVisa_Groups('/media/COSAS/scripts/BrainVisa_myTools/Brainvisa_nomenclature_sulci+STS.txt', 'right');
    
    Nstruct = size(SulcStr,1);
    
    % All labels in the same line
    SulcNames = '';
    for j = 1:Nstruct
        SulcNames = [SulcNames '''' deblank(SulcStr(j,:)) '''' ' '];
    end
    SulcNames(end) = [];
    
    %%  Preparing for Sulcal Length and Depth
    
% % %         BottonFullPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'bottom.ima'];
% % %         cmd = [ 'siGraph2Label -g ' RArgFile ' -o ' BottonFullPath ' -l ' SulcNames ' -b aims_bottom' ];
% % %         system(cmd)
% % %     
% % %         HullJunctPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'hull_junction.ima'];
% % %         cmd = [ 'siGraph2Label -g ' RArgFile ' -o ' HullJunctPath ' -l ' SulcNames ' -b aims_junction -s hull_junction' ];
% % %         system(cmd)
    % % %
    % % %
    % % %
    % % %     HullJuncttemp = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'hull_junction.nii'];
    % % %     cad = ['AimsFileConvert -i ' HullJunctPath ' -o ' HullJuncttemp];
    % % %     system(cad);
    % % %     Vt = spm_vol(HullJuncttemp);
    % % %     It = spm_read_vols(Vt);
    % % %     for i = 1:1
    % % %         It = imdilate(It,strel(ones(3,3,3)));
    % % %     end
    % % %     spm_write_vol(Vt,It);
    % % %     cad = ['AimsFileConvert -i ' Vt.fname ' -o ' HullJunctPath];
    % % %     system(cad);
    % % %     %
    % % %     BottonJuncttemp = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'bottom.nii'];
    % % %     cad = ['AimsFileConvert -i ' BottonFullPath ' -o ' BottonJuncttemp];
    % % %     system(cad);
    % % %     Vt = spm_vol(BottonJuncttemp);
    % % %     It = spm_read_vols(Vt);
    % % %     for i = 1:1
    % % %         It = imdilate(It,strel(ones(3,3,3)));
    % % %     end
    % % %     spm_write_vol(Vt,It);
    % % %     cad = ['AimsFileConvert -i ' Vt.fname ' -o ' BottonFullPath];
    % % %     system(cad);
    
    
    %%  End of Preparing for Sulcal Length and Depth
    
    
    % ---- Computing Gyral White Matter Spams
    
    %matlabpool open;
    for j = 1:Nstruct
        SulcName = deblank(SulcStr(j,:));
        
        % ----- Sulcal SPAM
        ind = find(ismember(SulcLabels,SulcName,'rows')); % Same Sulci
        try
            pointdistr = NodeNpoint(ind);
            pointdistr = 30;
            cad = ['RicSulcalSpan'...
                ' --sg ' RArgFile ...
                ' -o ' OutdirSPAM filesep subjId '_R' ...
                ' -n ' subjId ...
                ' --bv' ...
                ' --sm ' RSulcTmtktriFile ...
                ' --gm ' RHemiMesh ...
                ' --sn ' SulcName ...
                ' --maxa 1' ...
                ' --mina  0.9' ...
                ' --maxd 10' ...
                ' --mind 0.1' ...
                ' --mins ' num2str(floor(mean(pointdistr) - 0.5*mean(pointdistr)))...
                ' --sr 4' ...
                ' --nsd 3'];
            
            system(cad);
            
            
            BottonFullPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'bottom.ima'];
            cmd = [ 'siGraph2Label -g ' RArgFile ' -o ' BottonFullPath ' -l ' SulcName ' -b aims_bottom' ];
            system(cmd);
            
            HullJunctPath = 	[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'hull_junction.ima'];
            cmd = [ 'siGraph2Label -g ' RArgFile ' -o ' HullJunctPath ' -l ' SulcName ' -b aims_junction -s hull_junction' ];
            system(cmd);
            
            
            % %         HullJuncttemp = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'hull_junction.nii'];
            % %         cad = ['AimsFileConvert -i ' HullJunctPath ' -o ' HullJuncttemp];
            % %         system(cad);
            % %         Vt = spm_vol(HullJuncttemp);
            % %         It = spm_read_vols(Vt);
            % %         for i = 1:1
            % %             It = imdilate(It,strel(ones(3,3,3)));
            % %         end
            % %         spm_write_vol(Vt,It);
            % %         cad = ['AimsFileConvert -i ' Vt.fname ' -o ' HullJunctPath];
            % %         system(cad);
            % %         %
            % %         BottonJuncttemp = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'bottom.nii'];
            % %         cad = ['AimsFileConvert -i ' BottonFullPath ' -o ' BottonJuncttemp];
            % %         system(cad);
            % %         Vt = spm_vol(BottonJuncttemp);
            % %         It = spm_read_vols(Vt);
            % %         for i = 1:1
            % %             It = imdilate(It,strel(ones(3,3,3)));
            % %         end
            % %         spm_write_vol(Vt,It);
            % %         cad = ['AimsFileConvert -i ' Vt.fname ' -o ' BottonFullPath];
            % %         system(cad);
            
            
            
            % ----- Sulcal Length
            
            TempBaseNameL = [OutdirLENGTH filesep subjId '_R_' SulcName];
            cad = ['RiiSulcalLength --tvol '  HullJunctPath ' --bvol ' BottonFullPath ' --sm ' RSulcTmtktriFile	' --sn '  SulcName ...
                ' --sg ' RArgFile ...
                ' -o ' TempBaseNameL...
                ' -s  ' subjId ...
                ' --bv'];
            system(cad);
            
            TempBaseName = [OutdirDEPTH filesep subjId '_R_' SulcName];
            % ----- Sulcal Depth
            cad = ['RicSulcalDepth --sm ' TempBaseNameL '.mesh ' ...
                ' --tm ' TempBaseNameL '_top.mesh' ...
                ' --bm ' TempBaseNameL '_bottom.mesh' ...
                ' -o '   TempBaseName ...
                ' --ns  40' ...
                ' --nd 3 ' ...
                ' --fs 0.2' ...
                ' -s '   subjId ...
                ' --bv '];
            system(cad);
        end
        Files2Delete = strvcat(Files2Delete,BottonFullPath);
        Files2Delete = strvcat(Files2Delete,HullJunctPath);
    end
    %matlabpool close;
    for k = 1:size(Files2Delete,1)
        delete(deblank(Files2Delete(k,:)));
    end
    toc;
    %% %%%%%%%%%%%%%%%%%%%%%  End of Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%
    
    
end
%matlabpool('close');
return;


function [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktRII1ds] = Read_SulcalArgFiles(ArgFile);
%
% Syntax :
% [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktRII1ds] = Read_SulcalArgFiles(ArgFile);
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
            TmtktRII1ds(conttm) = str2num(temp{2});
        end
    end
    
end
fclose(fio);
return

return
