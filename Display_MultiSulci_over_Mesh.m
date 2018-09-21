function varargout = Display_MultiSulci_over_Mesh(varargin);
% %% ================== End of Checking Input parameters ================== %
BrainVisaDatabaseDir = '/media/COSAS/8-BrainVISADataBase-HCP';
IdFile = '/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';
% IdFile = '810843';

maxNumb = 5;
hemi= 'right';
if exist(IdFile, 'file')
    Ids = char(textread(IdFile,'%s'));
    Ids = Ids(2:end,:);
else
    Ids = IdFile;
end
OutSulcLines = '';
Idf = '';
switch hemi
    case 'left'
        hemiChar = 'L';
    case 'right'
        hemiChar = 'R';
    otherwise
        error('Wrong Hemisphere');
        return
end
label2Process = strvcat('F.P.O._right','S.C._right','F.Cal.ant.-Sc.Cal._right','S.T.s._right','F.C.M.post._right','F.C.M.ant._right');

%% ============================ Main Program ============================ %
% Colors
col = [[1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];[213 221 227;175 206 227;149 196 228; 120 186 230;87 172 231;24 146 232;6 73 119;244 207 154;244 192 117;244 179 86;244 161 43; 212 133 20;158 101 19; 113 71 12]/255];

Ns = size(Ids, 1);
if maxNumb > 0
    indsubj = [1:maxNumb];
else
    indsubj = [1:Ns];
end
Nl = size(label2Process,1);
cont = zeros(Nl,1);
for i = 1:length(indsubj)
    subjId = deblank(Ids(indsubj(i),:));
    
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(length(indsubj))]);
    
    % Reading Arg File
    LArgFile     = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep hemiChar subjId '_default_session_auto.arg' ];
    [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(LArgFile);
    NodeNpoint = NodeNpoint(1:length(TmtktriIds));
    
    % Reading Associeate Surface
    LTmtktriFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep hemiChar subjId '_default_session_auto.data' filesep 'aims_Tmtktri.gii' ];
    Surf = Read_Surface(LTmtktriFile);

    % Creating Spatial Transformation to Talairach Space
    spatTransf = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep  'registration' filesep 'RawT1-' subjId '_default_acquisition_TO_Talairach-MNI.trm'];
    tempMat = load(spatTransf);
    transfMat = [tempMat(2:4,:) tempMat(1,:)';[0 0 0 1] ];
    
    if i == 1;
        refSurfFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'segmentation' filesep 'mesh' filesep subjId '_' hemiChar 'white.gii' ];

        Surfref = Read_Surface(refSurfFile);
        Surfref = Surface_Rotation(Surfref,'rotMat',transfMat);
    end
    
    
    % For each sulci   
    for j =1:Nl
        Slabelo = deblank(label2Process(j,:)); % Sulci Label 
        L = length(Slabelo);
        inds = find(ismember(SulcLabels(:,1:L),Slabelo,'rows')); % Nodes that belong to the same sulci
        
        % Normalization to Talairach
        Surftemp = Surface_Rotation(Surf(TmtktriIds(inds)),'rotMat',transfMat);
        
        % Joinning Surfaces
        if ~isempty(inds)
            cont(j) = cont(j) + 1;
            if cont(j) == 1
                Surfja =  Compound_Surf(Surftemp);
                Surfja.Color = col(j,:);
                Surfj(j) = Surfja;
            else
                Surfja = Compound_Surf({Surfj(j);Surftemp});
                Surfja.Color = col(j,:);
                Surfj(j) = Surfja;
            end
        end
    end
    
    
    
    
    
% % %     RArgFile     = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.arg' ];
% % %     RTmtktriFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'aims_Tmtktri.mesh' ];    if ~exist(LArgFile,'file')|~exist(RArgFile,'file')
        
end
a = 1;

