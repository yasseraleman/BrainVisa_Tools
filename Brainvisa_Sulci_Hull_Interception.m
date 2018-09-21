function Brainvisa_Sulci_Hull_Interception(BrainVisaDatabaseDir,IdFile);

FreeSurferDatabaseDir = '/media/HCPData/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/COSAS/8-BrainVISADataBase-HCP';

%IdFile = '/media/MyDisk/PROCESSING_RESULTS/Total_Ids.txt';
IdFile ='/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';

Ids = char(textread(IdFile,'%s'));
Files2Delete = '';

Ns = size(Ids, 1);
for subj = 2:41
    subjId = deblank(Ids(subj,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(subj) ' of ' num2str(Ns)]);
    
    % ---- Creating Output Directories
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 'surface' ]);
    Outdir = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 'surface' ];
    
    sulcNames  = {'S.C._';'F.C.M.ant._';'F.C.M.post._';'F.C.M.r.AMS.ant._';'F.P.O._';'F.Cal.ant.-Sc.Cal._';'S.T.s._'};
    
    %% %%%%%%%%%%%%%%%%%%%%%%%% Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reading Sulci attributes form Arg file
    LArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.arg' ];
    [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(LArgFile);
    
    % Unifying sulcal labels
    SulcStr = unique(SulcLabels,'rows');
    
    SulcStr = char(strcat(sulcNames,repmat({'left'},[ length( sulcNames) 1])));
    SulcStrL = SulcStr;
    Ns = size(SulcStr,1);
    SulcNames = '';
    for j = 1:Ns
        SulcNames = [SulcNames '''' deblank(SulcStr(j,:)) '''' ' '];
    end
    SulcNames(end) = [];
    
    HullJunctPath   = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'hull_junction.nii'];
    HullJunctTrans  = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'L' subjId '_default_session_auto.data' filesep 'hull_junction_trans.txt'];
    cmd = [ 'siGraph2Label -g ' LArgFile ' -o ' HullJunctPath  ' -ot ' HullJunctTrans  ' -l ' SulcNames ' -a label -b aims_junction -s hull_junction' ];
    system(cmd)
    
    
    [sNames,sLabels] = textread(HullJunctTrans,'%s%u');
    
    
    % %     ind2Keep = find(sLabels);
    % %     sLabels = sLabels(ind2Keep);
    % %     sNames = sNames(ind2Keep);
    % %     [~,b] = ismember(SulcStr,sNames);
    % %     order = nonzeros(b);
    % %     sNames = sNames(order);
    % %     clustIds = sLabels(order);
    
    Vhull = spm_vol(HullJunctPath);
    I = spm_read_vols(Vhull);
    
    voxsize = sum(sqrt(Vhull.mat(1:3,1:3).^2));
    
    Graph = mask2graph(I,voxsize);
    ind = find(I);
    
    [X,Y,Z] = ind2sub(size(I),ind);
    
    tempVar = [Vhull.mat*[[X Y Z] ones(length(X),1)]']';
    mmCoords = tempVar(:,1:3);
    
    
    % % % % %     col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
    % % % % %     figure;
    Nstruct = size(SulcStr,1);
    for i = 1:Nstruct
        inpos = find(ismember(sNames, cellstr(SulcStr(i,:))));
        if sLabels(inpos) ~=0
            ind2 = find(I(ind) == sLabels(inpos));
            newCoords = mmCoords(ind2,:);
            subGraph = Graph(ind2,ind2);
            
            geoDist = graphallshortestpaths(subGraph);
            eucDist = dist(newCoords');
            totDist = geoDist +  eucDist;
            sulcLength = 0;
            
            LabNet = Label_Graph_Components(subGraph);
            Nclust = max(LabNet(:));
            for j = 1:Nclust
                [Xind3,Yind3] = find(LabNet == j);
                ind3 = unique([Xind3(:);Yind3(:)]);
                tempDistMat  = totDist(ind3,ind3);
                subGraphind3 = subGraph(ind3,ind3);
                
                [XextremCoord,YextremCoord] = find(tempDistMat == max(tempDistMat(:)));
                [sLength,temPath] = graphshortestpath(subGraphind3,XextremCoord(1),YextremCoord(1));
                tempCoords = newCoords(ind3,:);
                resLine = [smooth(tempCoords(temPath,1)) smooth(tempCoords(temPath,2)) smooth(tempCoords(temPath,3))];
                
                topCoords = fitCurveTo3DPts(resLine, resLine(1,:), resLine(end,:),floor(sLength/.5), 0 ); % Fitting a 3D curve
                
                sulcLength = sulcLength + sum(sqrt(sum((topCoords(1:end-1,:) - topCoords(2:end,:)).^2,2)));
                % % % % %                         plot3(tempCoords(:,1),tempCoords(:,2),tempCoords(:,3),'.','Color',col(i,:),'Markersize',20);
                % % % % %                         axis image
                % % % % %                         hold on
                % % % % %                         plot3(tempCoords([XextremCoord YextremCoord],1),tempCoords([XextremCoord YextremCoord],2),tempCoords([XextremCoord YextremCoord],3),'.w','Markersize',50);
                % % % % %                         plot3(topCoords(:,1),topCoords(:,2),topCoords(:,3),'-','Color',[1 1 1] -col(i,:),'Linewidth',5);
            end
        else
            sulcLength = 0;
        end
        resMatL(subj-1,i) = sulcLength;
    end
    
    %% %%%%%%%%%%%%%%%%%%%%% End of Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%  Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Reading Sulci attributes form Arg file
    RArgFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.arg' ];
    [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(RArgFile);
    
    % Unifying sulcal labels
    SulcStr = unique(SulcLabels,'rows');
    
    SulcStr = char(strcat(sulcNames,repmat({'right'},[ length( sulcNames) 1])));
    SulcStrR = SulcStr;
    Ns = size(SulcStr,1);
    SulcNames = '';
    for j = 1:Ns
        SulcNames = [SulcNames '''' deblank(SulcStr(j,:)) '''' ' '];
    end
    SulcNames(end) = [];
    
    HullJunctPath   = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'hull_junction.nii'];
    HullJunctTrans  = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'default_analysis' filesep 'folds' filesep '3.1' filesep 'default_session_auto' filesep 'R' subjId '_default_session_auto.data' filesep 'hull_junction_trans.txt'];
    cmd = [ 'siGraph2Label -g ' RArgFile ' -o ' HullJunctPath  ' -ot ' HullJunctTrans  ' -l ' SulcNames ' -a label -b aims_junction -s hull_junction' ];
    system(cmd)
    
    
    [sNames,sLabels] = textread(HullJunctTrans,'%s%u');
    %     ind2Keep = find(sLabels);
    %     sLabels = sLabels(ind2Keep);
    %     sNames = sNames(ind2Keep);
    %     [~,b] = ismember(SulcStr,sNames);
    %     order = nonzeros(b);
    %     sNames = sNames(order);
    %     clustIds = sLabels(order);
    
    Vhull = spm_vol(HullJunctPath);
    I = spm_read_vols(Vhull);
    
    voxsize = sum(sqrt(Vhull.mat(1:3,1:3).^2));
    
    Graph = mask2graph(I,voxsize);
    ind = find(I);
    
    [X,Y,Z] = ind2sub(size(I),ind);
    
    tempVar = [Vhull.mat*[[X Y Z] ones(length(X),1)]']';
    mmCoords = tempVar(:,1:3);
    
    
    col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
    
    figure;
    %     clustIds = nonzeros(unique(I(ind)));
    Nstruct = size(SulcStr,1);
    for i = 1:Nstruct
        inpos = find(ismember(sNames, cellstr(SulcStr(i,:))));
        if sLabels(inpos) ~=0
            ind2 = find(I(ind) == sLabels(inpos));
            
            newCoords = mmCoords(ind2,:);
            subGraph = Graph(ind2,ind2);
            
            geoDist = graphallshortestpaths(subGraph);
            eucDist = dist(newCoords');
            totDist = geoDist +  eucDist;
            sulcLength = 0;
            
            LabNet = Label_Graph_Components(subGraph);
            Nclust = max(LabNet(:));
            for j = 1:Nclust
                [Xind3,Yind3] = find(LabNet == j);
                ind3 = unique([Xind3(:);Yind3(:)]);
                tempDistMat  = totDist(ind3,ind3);
                subGraphind3 = subGraph(ind3,ind3);
                
                [XextremCoord,YextremCoord] = find(tempDistMat == max(tempDistMat(:)));
                [sLength,temPath] = graphshortestpath(subGraphind3,XextremCoord(1),YextremCoord(1));
                tempCoords = newCoords(ind3,:);
                resLine = [smooth(tempCoords(temPath,1)) smooth(tempCoords(temPath,2)) smooth(tempCoords(temPath,3))];
                
                topCoords = fitCurveTo3DPts(resLine, resLine(1,:), resLine(end,:),floor(sLength/.5), 0 ); % Fitting a 3D curve
                
                sulcLength = sulcLength + sum(sqrt(sum((topCoords(1:end-1,:) - topCoords(2:end,:)).^2,2)));
                            plot3(tempCoords(:,1),tempCoords(:,2),tempCoords(:,3),'.','Color',col(i,:),'Markersize',20);
                            axis image
                            hold on
                            plot3(tempCoords([XextremCoord YextremCoord],1),tempCoords([XextremCoord YextremCoord],2),tempCoords([XextremCoord YextremCoord],3),'.w','Markersize',50);
                            plot3(topCoords(:,1),topCoords(:,2),topCoords(:,3),'-','Color',[1 1 1] -col(i,:),'Linewidth',5);
            end
        else
            sulcLength = 0;
        end
        resMatR(subj-1,i) = sulcLength;
    end
    
    %% %%%%%%%%%%%%%%%%%%%%% End of Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%
    
end

   structNames =  {'Left-CentralSulcus(mm)' ,'Left-CingulateSulcus(mm)' ,'Left-ParietoOccipitalFissure(mm)' ,'Left-CalcCarineFissure(mm)' , 'Left-SuperiorTemporalSulcus(mm)',...
                          'Right-CentralSulcus(mm)','Right-CingulateSulcus(mm)','Right-ParietoOccipitalFissure(mm)','Right-CalcCarineFissure(mm)', 'Right-SuperiorTemporalSulcus(mm)'};     

resMatT = [ resMatL(:,1) sum(resMatL(:,2:4),2) resMatL(:,5:7) resMatR(:,1) sum(resMatR(:,2:4),2) resMatR(:,5:7)];
vt = [structNames;num2cell(resMatT)];

cads = repmat(' ',[size(vt,1) 1]);
for i = 1:size(vt,2)
    cadi = [structNames(i);cellstr(num2str(cell2mat(vt(2:end,i))))];
    cads = [cads char(cadi) repmat(';   ',[size(vt,1) 1])];
end

OutFile = [BrainVisaDatabaseDir filesep 'STATS_Results' filesep 'Hull_Major_Sulci_Interception.txt'];
fid = fopen(OutFile,'wt');
for i = 1:size(cads,1)
    line = cads(i,:);
    %     ind = strfind(line,'.');
    %     if ~isempty(ind)
    %         line(ind) = ',';
    %     end
    fprintf(fid, '%s\n', line);
end
fclose all;

HullSulcInt.Results = resMatT;
HullSulcInt.StructNames = char(structNames);
save( [BrainVisaDatabaseDir filesep 'STATS_Results' filesep 'Hull_Major_Sulci_Interception.mat'],'HullSulcInt');




return;


function Graph = mask2graph(I,voxsize);
% I = logical(D);
[X Y Z] = meshgrid(-1:1,-1:1,-1:1);Neib = [X(:) Y(:) Z(:)];
Neib(sum(logical(Neib)')==0,:) = [];
[mx,my,mz] = size(I);
ind = find(I);
[X, Y, Z] = ind2sub(size(I),ind);
Xo = 0;
Yo = 0;
Ndist = 0;
Np = length(ind);
B = sparse(Np,Np);
for j = 1:size(Neib,1)
    indorig = ind;
    CoordN = [X Y Z];
    indpos = [1:Np]';
    indx = find((X+Neib(j,1)<1)|(X+Neib(j,1)>mx));
    indy = find((Y+Neib(j,2)<1)|(Y+Neib(j,2)>my));
    indz = find((Z+Neib(j,3)<1)|(Z+Neib(j,3)>mz));
    inds2del = unique([indx;indy;indz]);
    CoordN(inds2del,:)  = [];
    indorig(inds2del,:) = [];
    
    
    indneigh = sub2ind(size(I),CoordN(:,1)+Neib(j,1), CoordN(:,2)+Neib(j,2), CoordN(:,3)+Neib(j,3));
    inds2del = find(I(indneigh) == 0);
    indneigh(inds2del) = [];
    indorig(inds2del) = [];
    
    
    Xo = [Xo;indpos(ismember(ind,indorig))];
    Yo = [Yo;indpos(ismember(ind,indneigh))];
    
    distad = ones(length(indorig),1).*(sqrt(dot(Neib(j,:).^2,voxsize.^2))); % Creating Graph with Distance Information
    
    Ndist = [Ndist;distad];
end
indtemp =ind;
pos = [1:length(ind)];
Xo(1,:) =[];
Yo(1,:) =[];
Ndist(1,:) = [];
Ndist = double(Ndist);
Graph = sparse([Xo;Yo],[Yo;Xo],[Ndist;Ndist ]);
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
