function varagout = Brainvisa_Sulcal_Basins_Extraction(varargin);


BrainVisaDatabaseDir = '/media/COSAS/8-BrainVISADataBase-HCP';
subjID               = '100206';
HemiChar = 'L';
sulcNames  = {'S.C._';'S.C.sylvian._';'F.C.M.ant._';'F.C.M.post._';'F.C.M.r.AMS.ant._';'F.P.O._';'S.T.s._';'F.Cal.ant.-Sc.Cal._'};


nodeLabels = char(strcat(sulcNames,repmat({'left'},[ length( sulcNames) 1])));

[groupsIds, groupsNames, correctNodes, labelsIds] = Detecting_BrainVisa_Groups(which('Brainvisa_nomenclature_sulci+STS.txt'), 'left');

varargout = Brainvisa_Sulcal_Basins_Extraction_Hemisphere(BrainVisaDatabaseDir, subjID, HemiChar, groupsIds, groupsNames, correctNodes, labelsIds);


function varargout = Brainvisa_Sulcal_Basins_Extraction_Hemisphere(varargin);

BrainVisaDatabaseDir = varargin{1};
subjID               = varargin{2};
HemiChar             = varargin{3};
groupsIds            = varargin{4};
groupsNames          = varargin{5};
correctNodes         = varargin{6}; 
labelsIds            = varargin{7};

sts2look = [5 7 8 9 10];

opts.distthresh  = 5; % mm. Maximum Sulcal Width
opts.epsilon = 10^-5;          % Tolerance to estimate barycentric coordinates
opts.verbose = 1;              % Maximum Sulcal Width

%% %%%%%%%%%%%%%%%%%% For each Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

pialMesh    =     [BrainVisaDatabaseDir filesep 'subjects' filesep subjID filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep subjID '_' HemiChar 'hemi.gii' ];
whiteMesh   =     [BrainVisaDatabaseDir filesep 'subjects' filesep subjID filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep subjID '_' HemiChar 'white.gii' ];
dpfFile     =     [BrainVisaDatabaseDir filesep 'subjects' filesep subjID filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep 'surface_analysis' filesep subjID '_' HemiChar 'white_DPF.gii' ];


SulcTmtktri =     [BrainVisaDatabaseDir filesep 'subjects' filesep subjID filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar subjID '_default_session_auto.data' filesep 'aims_Tmtktri.gii' ];
roiArgFile  =     [BrainVisaDatabaseDir filesep 'subjects' filesep subjID filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep HemiChar subjID '_default_session_auto.arg' ];

% Loading Surfaces
Surfwhite = Read_Surface(whiteMesh);
Surfpial = Read_Surface(pialMesh);

Surfsulci = Read_Surface(SulcTmtktri);


% Read Arg Filoe
[NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(roiArgFile);
NodeNpoint = NodeNpoint(1:length(TmtktriIds));



%% ====================== Processing each sulcus NODE ================== %%

failids = 0;
Is = zeros(size(Surfwhite.SurfData.vertices,1),1);
Nsulc = length(sts2look);


[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(Surfwhite.SurfData.vertices,double(Surfwhite.SurfData.faces));
curvmap =  -1*Cmean;


labSulbasin = Computing_Sulcal_Basins(Surfwhite, curvmap);


HullSurfMat = Compute_Hull_from_Surface(Surfwhite,.4);
opts.verb = 0;opts.nsub = 2;
[vertex,faces] = perform_mesh_subdivision(HullSurfMat.SurfData.vertices',HullSurfMat.SurfData.faces',opts.nsub,opts);
HullSurfMat.SurfData.faces = faces';
HullSurfMat.SurfData.vertices = vertex';




dpfmap = gifti(dpfFile);
dpfmap = double(dpfmap.cdata);

sLines = [0 0 0];
for i = 1:Nsulc
    Is = zeros(size(Surfwhite.SurfData.vertices,1),1);
    tempVar = correctNodes(labelsIds == sts2look(i),:);
     inds = find(ismember(cellstr(SulcLabels),cellstr(tempVar)));
    

    Nn = length(inds); % Number of Nodes in the sulci mesh
    for k = 1:Nn
        nodeSurf = Surfsulci(TmtktriIds(inds(k)));
        
        [Ist] = Compute_Node_Projection_over_Surf(Surfwhite, nodeSurf);
        Is(find(Ist)) = TmtktriIds(inds(k));
        
    end
    
% % % % %     ind = find(Is == 0);
    Surft = Surfwhite;
% % % % %     It = Is*0;
% % % % %     It(ind) = 1;
% % % % %     Surft.Is = It;
% % % % %     [labid] = Recur_Corr(Surft,0,zeros(size(Surft.Is)),1);
% % % % %     labid(ind) = 0;
% % % % %     
% % % % %     temp = accumarray(nonzeros(labid),nonzeros(labid)*0+1);
% % % % %     cont = 0;
% % % % %     opts.iterations = 1;
% % % % %     opts.fill = 1;
% % % % %     opts.reminter = 1;
% % % % %     while length(temp)~=1
% % % % %         cont = cont + 1;
% % % % %         [labid] = Dilate_Surface_Label(Surft, labid);
% % % % %         
% % % % %         
% % % % %         
% % % % %         ind = find(labid == 0);
% % % % %         It = labid*0;
% % % % %         It(ind) = 3000000;
% % % % %         Surft.Is = It;
% % % % %         [labid] = Recur_Corr(Surft,0,zeros(size(Surft.Is)),1);
% % % % %         labid(ind) = 0;
% % % % %         
% % % % %         temp = accumarray(nonzeros(labid),nonzeros(labid)*0+1);
% % % % %         ind2del = find(temp < max(temp)*.1);
% % % % %         labid(ismember(labid,ind2del)) = 0;
% % % % %         temp = accumarray(nonzeros(labid),nonzeros(labid)*0+1);
% % % % %         
% % % % %     end

Is = Is.*logical(labSulbasin);
ind = find(Is);
labid = labSulbasin;
b = accumarray(labid(ind),labid(ind)*0+1);
indpos = find(b <30);
labid(ismember(labid,indpos)) = 0;

sts = nonzeros(unique(labid(ind)));

labid(ismember(labid,sts)==0) = 0;



%     opts.iterations = 1;
%     [labid] = Erode_Surface_Label(Surft, labid, opts);
%     
%     
%     
%     
%     
% 
%     
%     ind2del = find(curvmap <= 0);
    
    
    Surft.Is = logical(labid);
    
    [Surft,~,indCell] = Extract_Sub_Surface(Surft, 1);
    ind = indCell{1};
    if ~isempty(ind)
        Surft.Is = curvmap(ind);
        % Reorganicing surface
        curvmapt = Surft.Is;         % Temporal curvature map
        curvmapt = curvmapt + abs(min(curvmapt)) +eps;
        
        
        locbottom = dsearchn(HullSurfMat.SurfData.vertices,Surft.SurfData.vertices);
        depthmapt = sqrt(sum((HullSurfMat.SurfData.vertices(locbottom,:)- Surft.SurfData.vertices).^2,2));

        opts.depthmap = depthmapt;
        % Extracting Sulcal Line
        [skelEdges] = Extracting_Sulcal_Skeleton(Surft, curvmapt,opts);
        skelEdges = skelEdges(find(skelEdges(:,3) == 1),:);
        

        
% % % % %         addFaces = Surft.SurfData.faces(find(sum(ismember(Surft.SurfData.faces,skelEdges(:,1:2)),2) ==3),:);
% % % % %         addEdges = [addFaces(:,1) addFaces(:,2);addFaces(:,2) addFaces(:,3);addFaces(:,1) addFaces(:,3)];
% % % % %         
% % % % %         Nv = length(Surft.SurfData.vertices); % Graph dimension
% % % % %         Graph = sparse(Nv,Nv); % Creating an empty sparse matrix dimension
% % % % %         allEdges = [skelEdges(:,1:2);addEdges];
% % % % %         X = (Surft.SurfData.vertices(allEdges(:,1),1) - Surft.SurfData.vertices(allEdges(:,2),1)).^2;
% % % % %         Y = (Surft.SurfData.vertices(allEdges(:,1),2) - Surft.SurfData.vertices(allEdges(:,2),2)).^2;
% % % % %         Z = (Surft.SurfData.vertices(allEdges(:,1),3) - Surft.SurfData.vertices(allEdges(:,2),3)).^2;
% % % % %         disEdge = sqrt(X + Y + Z);
% % % % %         ind2graph = sub2ind(size(Graph),allEdges(:,1),allEdges(:,2));
% % % % %         Graph(ind2graph) = disEdge;
% % % % %         start_point = skelEdges(1,1);
% % % % %         end_point = skelEdges(end,2);
% % % % %         [DIST, PATHS]=graphshortestpath(Graph,start_point,end_point);
% % % % %         PATHS = PATHS(:);
% % % % %         
% % % % %         skelEdges = [PATHS(1:end-1) PATHS(2:end)];
% % % % %         mainPathPoints = [skelEdges(:,1);skelEdges(end,1)];
% % % % %         [indp,locpeaks] = findpeaks(curvmapt(mainPathPoints));
% % % % %         
% % % % %         tempVar = remove_outliers(indp);
% % % % %         ind2rem = find(tempVar == Inf);
% % % % %         locpeaks(ind2rem) = [];
% % % % %         range(1) = min(locpeaks);
% % % % %         range(2) = max(locpeaks);
% % % % %         
% % % % %         skelEdges = skelEdges(find(sum(ismember(skelEdges(:,1:2),mainPathPoints(range(1):range(2))),2)==2),:);
        
        
        if ~isempty(skelEdges)
            sLines = [sLines;ind(skelEdges(:,1)) ind(skelEdges(:,2)) ones(size(skelEdges,1),1)*i] ;
        end
        %         if ~isempty(skelEdges)
        %             [skelEdges] = Branch_Labelling(Surft, skelEdges);
        %             MaxCurvPath = Compute_Max_Curvature_Path(Surft, skelEdges);
        %             sLines = [sLines;[ind(skelEdges(:,1:2)) zeros(size(skelEdges,1),1) repmat(i,[size(skelEdges,1) 1])]]; % Saving Branches. Branches are saved
        %             TotalMaxPath  = [TotalMaxPath; [ind(MaxCurvPath(:,1:2)) zeros(size(MaxCurvPath,1),1) repmat(i,[size(MaxCurvPath,1) 1])]]; % Saving Branches. Branches are saved
        %         end
    end
    
    
    
    
    
    
end
sLines(1,:) = [];


[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(Surfwhite.SurfData.vertices,double(Surfwhite.SurfData.faces));
labSulbasin = Computing_Sulcal_Basins(Surfwhite, -1*Cmean);

Is = Is.*logical(labSulbasin);
ind = find(Is);
b = accumarray(labSulbasin(ind),labSulbasin(ind)*0+1);
indpos = find(b <30);
labSulbasin(ismember(labSulbasin,indpos)) = 0;

sts = nonzeros(unique(labSulbasin(ind)));

labSulbasin(ismember(labSulbasin,sts)==0) = 0;


for i = 1:Nsul
end



return;