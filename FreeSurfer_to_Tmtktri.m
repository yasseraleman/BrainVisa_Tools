function FreeSurfer_to_Tmtktri(FreeSurf,Id,SulcArgFile);

% FreeSurferDatabaseDir = '/media/HCPData/5-freesurfer_processing';
% BrainVisaDatabaseDir = '/media/COSAS/8-BrainVISADataBase-HCP/';
% IdFile = '/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';
% aparcId = 'aparc';


FreeSurferDatabaseDir = '/media/Data/PROCESSING_RESULTS/URNC/Longitudinal_analysis/Structural/5-freesurfer_processing/New_annot/subjects_long_temp/';
BrainVisaDatabaseDir = '/media/COSAS/URNC-BrainVISADataBase/';
IdFile = '/media/Data/PROCESSING_RESULTS/URNC/Longitudinal_analysis/Structural/5-freesurfer_processing/New_annot/subjects_long_temp/URNC_Long_IDs.txt';
aparcId = 'aparc8regionsBV';



if exist(IdFile,'file')
    Ids = char(textread(IdFile,'%s'));
    
else
    Ids = IdFile;
end
%Ids = '1Test_HCP_899885-20140807-T1wMPR1';
Ns = size(Ids, 1);
for i = 1:Ns
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    %subjId = '1Test_HCP_899885-20140807-T1wMPR1';
    
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII']);
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'mesh']);
    
    %% ==================== Lobes Parcellation ============================== %
    % Selecting Frontal regions
    FroIds = sort([1028 1003 1027 1018 1019 1020 1012 1014 1024 1017 1032 1026 1002]);
    FroIdsR = FroIds + 1000;
    
    % Selecting Temporal regions
    TempIds = sort([1009 1015 1030 1001 1007 1034 1006 1033 1016]);
    TempIdsR = TempIds + 1000;
    
    % Selecting Parietal regions
    ParIds = sort([1029 1008 1031 1022 1025 1010 1023]);
    ParIdsR = ParIds + 1000;
    
    % Selecting Occipital regions
    OccIds = sort([1011 1013 1005 1021]);
    OccIdsR = OccIds + 1000;
    
    % Selecting Insula regions
    InsIds = [1035];
    InsIdsR = [2035];
    %% ================ End of Lobes Parcellation ======================= %
    
    %% ===================== Mandatory Files ============================ %
    lhsurf = [FreeSurferDatabaseDir filesep subjId filesep 'surf' filesep 'lh.white'];
    rhsurf = [FreeSurferDatabaseDir filesep subjId filesep 'surf' filesep 'rh.white'];
    
    lhannot = [FreeSurferDatabaseDir filesep subjId filesep 'label' filesep 'lh.' aparcId '.annot'];
    rhannot = [FreeSurferDatabaseDir filesep subjId filesep 'label' filesep 'rh.' aparcId '.annot'];
    Taltransf = [FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
    Imfile = [FreeSurferDatabaseDir filesep subjId filesep 'tmp' filesep 'T1.nii'];
    cad = ['mri_convert -i ' FreeSurferDatabaseDir filesep subjId filesep 'mri' filesep 'T1.mgz -o ' Imfile] ;
    system(cad);
    %% ================== End of Mandatory Files ======================== %
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%  Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
    LArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'L' subjId '.arg'];
    
    %% ===================== Reading Arg File =========================== %
    fio = fopen(LArgFile,'rt');lines = '';cont = 0;
    contmax = 0;
    contmin = 0;
    contvox = 0;
    
    while 1
        cont = cont + 1;
        line = fgetl(fio);
        if ~ischar(line),   break,   end
        linet= line;
        if strfind(linet,'boundingbox_max')&(contmax == 0);
            contmax = contmax + 1;
            bbmaxline = linet;
        end
        if strfind(linet,'boundingbox_min')&(contmin == 0);
            contmin = contmin + 1;
            bbminline = linet;
        end
        if strfind(linet,'voxel_size')&(contvox == 0);
            contvox = contvox + 1;
            voxsizeline = linet;
        end
    end
    fclose(fio);
    temp = strread(bbmaxline,'%s');
    bbmaxline = [char(temp(1)) ' ' char(temp(2)) ' ' char(temp(3)) ' ' char(temp(4))];
    temp = strread(bbminline,'%s');
    bbminline = [char(temp(1)) ' ' char(temp(2)) ' ' char(temp(3)) ' ' char(temp(4))];
    temp = strread(voxsizeline,'%s');
    voxsizeline = [char(temp(1))  repmat( ' ',[1 5]) ' ' char(temp(2)) ' ' char(temp(3)) ' ' char(temp(4))];
    %% =================End of Reading Arg File ========================= %
    
    
    %% Reading Data
    %----- Reading Surfaces
    [OutFiles, SurfF] = Exp_Surf(lhsurf, '0', '','', 'imp','n');
    Surfwl= SurfF{1};
    Surfwl.Name = 'LH.WHITE';
    cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
    Surfwl.SurfData.vertices =Surfwl.SurfData.vertices+repmat(cras,[size(Surfwl.SurfData.vertices,1) 1]);
    Surfwl = freesCS2brainvisaCS(Surfwl,Imfile,'f2b');
    %Out_mesh = save_mesh(Surfwl,[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep subjId '_Lwhite.mesh']);
    
    %----- Reading Annot Files
    [txtl,ctabl] = read_cfiles(lhannot);
    ctabl.table = [ctabl.table 1000+[0:size(ctabl.table,1)-1]' ];
    tempname = char(ctabl.struct_names);
    indu = find(ismember(tempname(:,1:7),'unknown','rows') == 1);
    indcc = find(ismember(tempname(:,1:14),'corpuscallosum','rows') == 1);
    lindwwall = [ctabl.table([indu indcc],5) ;0];          % Left Hemisphere Medial wall indexes
    indwall = find(ismember(txtl, lindwwall));
    txtl(indwall) = 0;
    ctabl.table([indu indcc],:) = [];
    ctabl.struct_names([indu indcc]) = [];
    ctabl.table = [255 255 255 0 0 0;ctabl.table];
    ctabl.struct_names= [{'medial_wall'};ctabl.struct_names];
    
    BaseName = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep subjId '_Lgyri_default_session_auto' ];
    mkdir([BaseName '.data']);
    
    % Saving Arg File
    fiw = fopen([BaseName '.arg'],'wt');
    fprintf(fiw,'%s\n','# graph 1.0');
    fprintf(fiw,'\n');
    fprintf(fiw,'%s\n','*BEGIN GRAPH RoiArg');
    fprintf(fiw,'%s\n', bbmaxline);
    fprintf(fiw,'%s\n', bbminline);
    fprintf(fiw,'%s\n', voxsizeline);
    fprintf(fiw,'%s\n','filename_base   *');
    fprintf(fiw,'%s\n',[ 'roi.global.tri  roi aims_Tmtktri.mesh Tmtktri_label']);
    fprintf(fiw,'%s\n','type.global.tri roi.global.tri');
    fprintf(fiw,'\n');
    
    Nstruct = size(ctabl.table,1);
    for j = 1:Nstruct
        ind2remove = find(txtl~=ctabl.table(j,5));
        
        ind2keep = find(ismember(txtl, lindwwall) == 0);
        Surft = Surfwl;
        indfac2remove = find(sum(ismember(Surft.SurfData.faces,ind2remove),2) == 3);
        
        Surft.SurfData.faces(indfac2remove,:) = []; % Medial Wall
        [Surft] = Reorg_Surf(Surft);
        Surfa(j) = Surft; % Medial Wall
        Surfa(j).Name = char(ctabl.struct_names(j));
        
        
        
        fprintf(fiw,'%s\n',['*BEGIN NODE roi ' num2str(j)]);
        fprintf(fiw,'%s\n',['name          ' Surfa(j).Name]);
        fprintf(fiw,'%s\n',['Tmtktri_label ' num2str(j)]);
        fprintf(fiw,'%s\n',['roi_label     ' num2str(j-1)]);
        fprintf(fiw,'%s\n',['surface_area  ' num2str(0)]);
        fprintf(fiw,'%s\n','*END');
        fprintf(fiw,'\n');
    end
    fprintf(fiw,'%s','*END');
    fclose(fiw);
    
    Out_mesh = save_mesh(Surfa,[BaseName '.data' filesep 'aims_Tmtktri.mesh' ]);
    
    %% %%%%%%%%%%%%%%%%%%%% End of  Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%
    
    %% %%%%%%%%%%%%%%%%%%%%%%%% Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    RArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'R' subjId '.arg'];
    
    %% ===================== Reading Arg File =========================== %
    fio = fopen(RArgFile,'rt');lines = '';cont = 0;
    contmax = 0;
    contmin = 0;
    contvox = 0;
    
    while 1
        cont = cont + 1;
        line = fgetl(fio);
        if ~ischar(line),   break,   end
        linet= line;
        if strfind(linet,'boundingbox_max')&(contmax == 0);
            contmax = contmax + 1;
            bbmaxline = linet;
        end
        if strfind(linet,'boundingbox_min')&(contmin == 0);
            contmin = contmin + 1;
            bbminline = linet;
        end
        if strfind(linet,'voxel_size')&(contvox == 0);
            contvox = contvox + 1;
            voxsizeline = linet;
        end
    end
    fclose(fio);
    temp = strread(bbmaxline,'%s');
    bbmaxline = [char(temp(1)) ' ' char(temp(2)) ' ' char(temp(3)) ' ' char(temp(4))];
    temp = strread(bbminline,'%s');
    bbminline = [char(temp(1)) ' ' char(temp(2)) ' ' char(temp(3)) ' ' char(temp(4))];
    temp = strread(voxsizeline,'%s');
    voxsizeline = [char(temp(1))  repmat( ' ',[1 5]) ' ' char(temp(2)) ' ' char(temp(3)) ' ' char(temp(4))];
    %% =================End of Reading Arg File ========================= %
    
    
    %% Reading Data
    %----- Reading Surfaces
    [OutFiles, SurfF] = Exp_Surf(rhsurf, '0', '','', 'imp','n');
    Surfwr= SurfF{1};
    Surfwr.Name = 'RH.WHITE';
    cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
    Surfwr.SurfData.vertices =Surfwr.SurfData.vertices+repmat(cras,[size(Surfwr.SurfData.vertices,1) 1]);
    Surfwr = freesCS2brainvisaCS(Surfwr,Imfile,'f2b');
    delete(Imfile);
    %Out_mesh = save_mesh(Surfwr,[BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep subjId '_Rwhite.mesh']);
    
    %----- Reading Annot Files
    [txtr,ctabr] = read_cfiles(rhannot);
    ctabr.table = [ctabr.table 2000+[0:size(ctabr.table,1)-1]' ];
    tempname = char(ctabr.struct_names);
    indu = find(ismember(tempname(:,1:7),'unknown','rows') == 1);
    indcc = find(ismember(tempname(:,1:14),'corpuscallosum','rows') == 1);
    lindwwall = [ctabr.table([indu indcc],5) ;0];          % Left Hemisphere Medial wall indexes
    indwall = find(ismember(txtr, lindwwall));
    txtr(indwall) = 0;
    ctabr.table([indu indcc],:) = [];
    ctabr.struct_names([indu indcc]) = [];
    ctabr.table = [255 255 255 0 0 0;ctabr.table];
    ctabr.struct_names= [{'medial_wall'};ctabr.struct_names];
    
    BaseName = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep subjId '_Rgyri_default_session_auto' ];
    mkdir([BaseName '.data']);
    
    % Saving Arg File
    fiw = fopen([BaseName '.arg'],'wt');
    fprintf(fiw,'%s\n','# graph 1.0');
    fprintf(fiw,'\n');
    fprintf(fiw,'%s\n','*BEGIN GRAPH RoiArg');
    fprintf(fiw,'%s\n', bbmaxline);
    fprintf(fiw,'%s\n', bbminline);
    fprintf(fiw,'%s\n', voxsizeline);
    fprintf(fiw,'%s\n','filename_base   *');
    fprintf(fiw,'%s\n',[ 'roi.global.tri  roi aims_Tmtktri.mesh Tmtktri_label']);
    fprintf(fiw,'%s\n','type.global.tri roi.global.tri');
    fprintf(fiw,'\n');
    
    Nstruct = size(ctabr.table,1);
    for j = 1:Nstruct
        ind2remove = find(txtr~=ctabr.table(j,5));
        
        ind2keep = find(ismember(txtr, lindwwall) == 0);
        Surft = Surfwr;
        indfac2remove = find(sum(ismember(Surft.SurfData.faces,ind2remove),2) == 3);
        
        Surft.SurfData.faces(indfac2remove,:) = []; % Medial Wall
        [Surft] = Reorg_Surf(Surft);
        Surfa(j) = Surft; % Medial Wall
        Surfa(j).Name = char(ctabr.struct_names(j));
        
        
        
        fprintf(fiw,'%s\n',['*BEGIN NODE roi ' num2str(j)]);
        fprintf(fiw,'%s\n',['name          ' Surfa(j).Name]);
        fprintf(fiw,'%s\n',['Tmtktri_label ' num2str(j)]);
        fprintf(fiw,'%s\n',['roi_label     ' num2str(j-1)]);
        fprintf(fiw,'%s\n',['surface_area  ' num2str(0)]);
        fprintf(fiw,'%s\n','*END');
        fprintf(fiw,'\n');
    end
    fprintf(fiw,'%s','*END');
    fclose(fiw);
    
    Out_mesh = save_mesh(Surfa,[BaseName '.data' filesep 'aims_Tmtktri.mesh' ]);
    
    %% %%%%%%%%%%%%%%%%%%%% End of  Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%
    
    
end
return;
