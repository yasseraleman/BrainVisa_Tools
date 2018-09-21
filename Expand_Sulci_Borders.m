function Expand_Sulci_Borders(BrainVisaDatabaseDir,IdFile);

FreeSurferDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase';
%IdFile = '/media/MyDisk/PROCESSING_RESULTS/Total_Ids.txt';
IdFile = '/media/MyDisk/23PartIds.txt';

Ids = char(textread(IdFile,'%s'));
Ids = '1Test_HCP_899885-20140807-T1wMPR1';

Files2Delete = '';

Ns = size(Ids, 1);
for i = 1:Ns
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Ns)]);
    %subjId = '1Test_HCP_899885-20140807-T1wMPR1';
    tic;
    
        
    %% %%%%%%%%%%%%%%%%%%%%%%%  Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
       
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
        else
            
            Surft = Sulcal_Face_Labelling(Surf(j));
            Normals = Surft.SurfData.VertexNormals;
            [Trip] = Vert_Neibp(double(Surft.SurfData.faces),size(Surft.SurfData.vertices,1),size(Surft.SurfData.faces,1));
            Temp = sum(Trip);
            Trip(:,Temp==0) = [];
            temp = Trip(:,3:end);
            indz = find(temp == 0);
            
            Coord = [repmat(Trip(:,1),[size(temp,2) 1]) temp(:)];
            ind = find(sum(logical(Coord),2) == 2);
            Coord = Coord(ind,:);
            temp = sort(Coord')';
            Coord = unique(temp,'rows');
            
            T = Surft.Is(Coord);
            ind = find(T(:,1) - T(:,2) ~=0); % Crossing edges between walls
            Coord = Coord(ind,:);
            
            indd = accumarray(Coord(:),Coord(:)*0+1); % Points that only appears in a single edge
            ind = find(indd == 1);
            edge2rem = find(sum(ismember(Coord,ind),2) ~=0);
            Coord(edge2rem,:) = []; % Removing unconnected edges
            Surf(j).SurfData.vertices(Coord(:),:) = Surf(j).SurfData.vertices(Coord(:),:) + 2*Normals(Coord(:),:);
        end
    end
    
    % Removing empty Surfaces
    Surf(indel) = [];
    
    % Saving the new mesh surface
    LSulcTmtktriFile2Save = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'aims_Tmtktri_mod.mesh' ]
    save_mesh(Surf, LSulcTmtktriFile2Save);

    
end
return;



