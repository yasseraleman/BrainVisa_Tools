function varargout = BrainVisa_RII_stats_to_table(varargin);
%
% Syntax :
%  statTable = BrainVisa_RII_stats_to_table(BrainVisaDatabaseDir, IdFile, statTable);
%
% This script creates a stats table using the BrainVisa outputs.
%
% Input Parameters:
%       statFiles               : BrainVisa Stat Files (morpho*)
%       selIds                  : Id List
%
%
%
% Output Parameters:
%       statTable               : Stats table File
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%
% if nargin <2
%     error('Two Inputs are needed');
%     return
% end
% BrainVisaDatabaseDir = varargin{1};
% IdFile = varargin{2};
%
% if nargin == 3
%     statTable = varargin{3};
% else
%     tempFile = deblank(statFiles(1,:));
%     if exist(tempFile,'file')
%         mkdir([BrainVisaDatabaseDir filesep 'RII_Stat_Results']);
%         statTable = [BrainVisaDatabaseDir filesep 'RII_Stat_Results' filesep 'RII_General_Table_Orig.txt'];
%     else
%         error('Incorrect Stat File');
%         return;
%     end
% end
%% ======================= End of Input parameters  ======================%

BrainVisaDatabaseDir =  '/media/COSAS/8-BrainVISADataBase-HCP';
IdFile ='/media/HCPData/5-freesurfer_processing/freeSurferIds.txt';

outDir = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase/Sulci_Morpho_Metrics';
warning off;
mkdir(outDir);


%% ========================== Useful Nodes ============================== %

[groupsIds, groupsNamesL, correctNodesIL, labelsIdsL] = Detecting_BrainVisa_Groups(which('Brainvisa_nomenclature_sulci+STS.txt'), 'left');
colHeadersIL = correctNodesIL;
[groupsIds, groupsNamesR, correctNodesIR, labelsIdsR] = Detecting_BrainVisa_Groups(which('Brainvisa_nomenclature_sulci+STS.txt'), 'right');
colHeadersIR = correctNodesIR;



% ------------------ Lobar Grouping
% -- Left
labels = [ones(19,1);2*ones(11,1);3*ones(12,1);4*ones(7,1);5*ones(2,1);6*ones(7,1);7*ones(4,1);8;9];

% -- Right
% labels = [labels;labels + 9];

%% ========================== End of Useful Nodes ======================= %



repComa = [',';','];
Ids = char(textread(IdFile,'%s'));
Ns = size(Ids, 1);


measNamesL = '';
measNamesR = '';

Ns = 101;

for i = 1:Ns
    subjId = deblank(Ids(i,:));
    disp(['Proccessing Subject: ' subjId '  ==> ' num2str(i) ' of ' num2str(Ns)]);

    % ------------------ Output Folders
    OutdirLENGTH = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII1' filesep 'sulcallength'];
    OutdirDEPTH = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII1' filesep 'sulcaldepth'];
    OutdirSPAM = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII1' filesep 'sulcalspams'];
    
    for j = 1:size(correctNodesIL,1)
        
        %% %%%%%%%%%%%%%%%%%%%%%%%  Left Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%
        SulcNameL = deblank(correctNodesIL(j,:));
        
        % Extracting Length Values
        lengthFile = [OutdirLENGTH filesep subjId '_L_' SulcNameL '.txt'];
        if exist(lengthFile,'file')
            [id,lhbotLength,lhtopLength] = textread(lengthFile,'%s%f%f');
            lhmeanLength = mean([lhtopLength lhbotLength]);
            lhintercLength = lhtopLength;
        else
            lhtopLength = 0;
            lhbotLength = 0;
            lhmeanLength = 0;
            lhintercLength = 0;
        end
        
        
        % Extracting Depth Values
        depthFile = [OutdirDEPTH filesep subjId '_L_' SulcNameL '_SegDepth.txt'];
        if exist(depthFile,'file')
            [temp] = textread(depthFile,'%s');
            if length(temp) > 1
                temp1 = str2num(char(temp{2:end}));
                lhmaxDepth = max(temp1);
                lhmeanDepth = mean(temp1);
                lhminDepth = min(temp1);
                lhnDepth = length(temp1);
            else
                lhmaxDepth = 0;
                lhmeanDepth = 0;
                lhminDepth = 0;
                lhnDepth = 0;
            end
        else
            lhmaxDepth = 0;
            lhmeanDepth = 0;
            lhminDepth = 0;
            lhnDepth = 0;
        end
        
        
        % Extracting Spam Values
        spamFile = [OutdirSPAM filesep subjId '_L_' SulcNameL '.txt'];
        if exist(spamFile,'file')
            [id,lhmaxSpam,lhmeanSpam,lhminSpam, lhnSpam] = textread(spamFile,'%s%f%f%f%u');
        else
            lhmaxSpam = 0;
            lhmeanSpam = 0;
            lhminSpam = 0;
            lhnSpam = 0;
        end
        
        if i == 1
            measNamesL = strvcat(measNamesL,['geodDepthMax-'    SulcNameL '(mm)'],['geodDepthMin-'    SulcNameL '(mm)'],['geodDepthMean-'   SulcNameL '(mm)'],...
                ['topLength-'      SulcNameL '(mm)'],['bottomLength-'   SulcNameL '(mm)'],['intercepLength-' SulcNameL '(mm)'],['lengthMean-'     SulcNameL '(mm)'],...
                ['spamMax-'        SulcNameL '(mm)'],['spamMin-'        SulcNameL '(mm)'],['spamMean-'       SulcNameL '(mm)']);
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%% End of Left Hemisphere %%%%%%%%%%%%%%%%%%
        
        %% %%%%%%%%%%%%%%%%%%%%%%%  Right Hemisphere %%%%%%%%%%%%%%%%%%%%%%
        SulcNameR = deblank(correctNodesIR(j,:));
        
        % Extracting Length Values
        lengthFile = [OutdirLENGTH filesep subjId '_R_' SulcNameR '.txt'];
        if exist(lengthFile,'file')
            [id,rhbotLength,rhtopLength] = textread(lengthFile,'%s%f%f');
            rhmeanLength = mean([rhtopLength rhbotLength]);
            rhintercLength = rhtopLength;
        else
            rhtopLength = 0;
            rhbotLength = 0;
            rhmeanLength = 0;
            rhintercLength = 0;
        end
        
        
        % Extracting Depth Values
        depthFile = [OutdirDEPTH filesep subjId '_R_' SulcNameR '_SegDepth.txt'];
        if exist(depthFile,'file')
            [temp] = textread(depthFile,'%s');
            if length(temp) > 1
                temp1 = str2num(char(temp{2:end}));
                rhmaxDepth = max(temp1);
                rhmeanDepth = mean(temp1);
                rhminDepth = min(temp1);
                rhnDepth = length(temp1);
            else
                rhmaxDepth = 0;
                rhmeanDepth = 0;
                rhminDepth = 0;
                rhnDepth = 0;
            end
        else
            rhmaxDepth = 0;
            rhmeanDepth = 0;
            rhminDepth = 0;
            rhnDepth = 0;
        end
        
        
        % Extracting Spam Values
        spamFile = [OutdirSPAM filesep subjId '_R_' SulcNameR '.txt'];
        if exist(spamFile,'file')
            [id,rhmaxSpam,rhmeanSpam,rhminSpam, rhnSpam] = textread(spamFile,'%s%f%f%f%u');
        else
            rhmaxSpam = 0;
            rhmeanSpam = 0;
            rhminSpam = 0;
            rhnSpam = 0;
        end
        
        if i == 1
            measNamesR = strvcat(measNamesR,['geodDepthMax-'    SulcNameR '(mm)'],['geodDepthMin-'    SulcNameR '(mm)'],['geodDepthMean-'   SulcNameR '(mm)'],...
                ['topLength-'      SulcNameR '(mm)'],['bottomLength-'   SulcNameR '(mm)'],['intercepLength-' SulcNameR '(mm)'],['lengthMean-'     SulcNameR '(mm)'],...
                ['spamMax-'        SulcNameR '(mm)'],['spamMin-'        SulcNameR '(mm)'],['spamMean-'       SulcNameR '(mm)']);
        end
        %% %%%%%%%%%%%%%%%%%%% End of Right Hemisphere %%%%%%%%%%%%%%%%%%%%
        
        if j == 1
            tempVarL = [lhmaxDepth lhminDepth lhmeanDepth lhtopLength lhbotLength lhintercLength lhmeanLength lhmaxSpam lhminSpam lhmeanSpam];
            tempVarR = [rhmaxDepth rhminDepth rhmeanDepth rhtopLength rhbotLength rhintercLength rhmeanLength rhmaxSpam rhminSpam rhmeanSpam];
        else
            tempVarL = [tempVarL lhmaxDepth lhminDepth lhmeanDepth lhtopLength lhbotLength lhintercLength lhmeanLength lhmaxSpam lhminSpam lhmeanSpam];
            tempVarR = [tempVarR rhmaxDepth rhminDepth rhmeanDepth rhtopLength rhbotLength rhintercLength rhmeanLength rhmaxSpam rhminSpam rhmeanSpam];
        end
    end
    metricsMatL(i,:) = tempVarL(:)';
    metricsMatR(i,:) = tempVarR(:)';
end

[metricsMatL,measNamesL]  = Lobar_Stats(metricsMatL,groupsNamesL,labelsIdsL, measNamesL, 'Left-Hemisphere');

[metricsMatR,measNamesR]  = Lobar_Stats(metricsMatR,groupsNamesR,labelsIdsR, measNamesR, 'Right-Hemisphere');


%% %%%%%%%%%%%%%%%%%%%%%%%% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

Nvar = 10;
% Reading
L1 = cellstr(measNamesL);
R1 = cellstr(measNamesR);
Ls1 = reshape(L1,[Nvar length(L1)/Nvar]);
Rs1 = reshape(R1,[Nvar length(R1)/Nvar]);
TotNames = [Ls1;Rs1];
TotNames = TotNames(:);

for j = 1:Ns
    opts.subjid = deblank(Ids(i,:));
    
    opts.subjid = deblank(Ids(j,:));
    L1 = metricsMatL(j,:);
    R1 = metricsMatR(j,:);
    Ls1 = reshape(L1,[Nvar length(L1)/Nvar]);
    Rs1 = reshape(R1,[Nvar length(R1)/Nvar]);
    TotVariab = [Ls1;Rs1];
    TotVariab = cellstr(num2str(TotVariab(:)));
    
    Temp = [TotNames';TotVariab'];
    cad2print = '';
    for vari = 1:size(Temp,2)
        cad2print = [cad2print char(Temp(:,vari)) [',';',']];
    end
    cad2print(:,end) = [];
    cad2print = [strvcat('Subjects',opts.subjid) [',';','] cad2print];
    warning off;
    mkdir([ BrainVisaDatabaseDir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'RII1' filesep 'stats']);
    warning on;
    IstatTable = [ BrainVisaDatabaseDir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'RII1' filesep 'stats' filesep opts.subjid '_sulcal_auto_Stats.txt'];
    fidind = fopen(IstatTable,'wt');
    fprintf(fidind,'%s\n',cad2print(1,:));
    fprintf(fidind,'%s\n',cad2print(2,:));
    fclose(fidind);
end


return;

%%  ======================== Grouping into Lobes ======================= %%

%% ---  Detecting the number of subjects

function [metricsMat,measNames]  = Lobar_Stats(metricsMat,groupsNames,labelsIds, measNames, hemiHeader);

sts = unique(labelsIds);
Nlobes = length(sts);
for j = 1:Nlobes
    ind = find(labelsIds == sts(j));
    T = metricsMat(:,(ind-1)*10+1);
    LmaxDepth        = max(T,[],2);
    T = metricsMat(:,(ind-1)*10+2);
    LminDepth        = min(T,[],2);
    T = metricsMat(:,(ind-1)*10+3);
    LmeanDepth       = sum(T,2)./(sum(logical(T),2) + eps);
    LtopLength       = sum(metricsMat(:,(ind-1)*10+4),2);
    LbottomLength    = sum(metricsMat(:,(ind-1)*10+5),2);
    LinterLength      = sum(metricsMat(:,(ind-1)*10+6),2);
    LmeanLength      = sum(metricsMat(:,(ind-1)*10+7),2);
    T = metricsMat(:,(ind-1)*10+8);
    LmaxSPAM         = max(T,[],2);
    T = metricsMat(:,(ind-1)*10+9);
    LminSPAM         = min(T,[],2);
    T = metricsMat(:,(ind-1)*10+10);
    LmeanSPAM        = sum(T,2)./(sum(logical(T),2) + eps);
    if j == 1
        tempMatL = [LmaxDepth LminDepth LmeanDepth LtopLength LbottomLength LinterLength LmeanLength LmaxSPAM LminSPAM LmeanSPAM];
        measNamesL = strvcat(['geodDepthMax-' deblank(groupsNames(j,:))  '(mm)'],...
            ['geodDepthMin-' deblank(groupsNames(j,:))  '(mm)'],...
            ['geodDepthMean-'  deblank(groupsNames(j,:))  '(mm)'],...
            ['topLength-'  deblank(groupsNames(j,:))  '(mm)'],...
            ['bottomLength-' deblank(groupsNames(j,:))  '(mm)'],...
            ['intercepLength-' deblank(groupsNames(j,:))  '(mm)'],...
            ['lengthMean-' deblank(groupsNames(j,:))  '(mm)'],...
            ['spamMax-' deblank(groupsNames(j,:))  '(mm)'],...
            ['spamMin-' deblank(groupsNames(j,:))  '(mm)'],...
            ['spamMean-' deblank(groupsNames(j,:))  '(mm)']);
    else
        tempMatL = [tempMatL LmaxDepth LminDepth LmeanDepth LtopLength LbottomLength LinterLength LmeanLength LmaxSPAM LminSPAM LmeanSPAM];
        measNamesL = strvcat(measNamesL, ['geodDepthMax-' deblank(groupsNames(j,:))  '(mm)'],...
            ['geodDepthMin-' deblank(groupsNames(j,:))  '(mm)'],...
            ['geodDepthMean-'  deblank(groupsNames(j,:))  '(mm)'],...
            ['topLength-'  deblank(groupsNames(j,:))  '(mm)'],...
            ['bottomLength-' deblank(groupsNames(j,:))  '(mm)'],...
            ['intercepLength-' deblank(groupsNames(j,:))  '(mm)'],...
            ['lengthMean-' deblank(groupsNames(j,:))  '(mm)'],...
            ['spamMax-' deblank(groupsNames(j,:))  '(mm)'],...
            ['spamMin-' deblank(groupsNames(j,:))  '(mm)'],...
            ['spamMean-' deblank(groupsNames(j,:))  '(mm)']);
    end
end
measNames = strvcat(measNamesL,measNames);
metricsMat = [tempMatL metricsMat];

% Creating Hemisphere Values
measNames = strvcat(['geodDepthMax-' hemiHeader  '(mm)'],...
    ['geodDepthMin-' hemiHeader  '(mm)'],...
    ['geodDepthMean-'  hemiHeader  '(mm)'],...
    ['topLength-'  hemiHeader  '(mm)'],...
    ['bottomLength-' hemiHeader  '(mm)'],...
    ['intercepLength-' hemiHeader  '(mm)'],...
    ['lengthMean-' hemiHeader  '(mm)'],...
    ['spamMax-' hemiHeader  '(mm)'],...
    ['spamMin-' hemiHeader  '(mm)'],...
    ['spamMean-' hemiHeader  '(mm)'],...
    measNames);
metricsMat = [max(tempMatL(:,1:10:end),[],2) ...
    min(tempMatL(:,2:10:end),[],2) ...
    mean(tempMatL(:,3:10:end),2) ...
    sum(tempMatL(:,4:10:end),2) ...
    sum(tempMatL(:,5:10:end),2) ...
    sum(tempMatL(:,6:10:end),2) ...
    sum(tempMatL(:,7:10:end),2) ...
    max(tempMatL(:,8:10:end),[],2) ...
    min(tempMatL(:,9:10:end),[],2)...
    mean(tempMatL(:,10:10:end),2) ...
    metricsMat];
return






% % % % % % % % 
% % % % % % % % 
% % % % % % % % % --------------- Grouping into Lobes
% % % % % % % % sts = unique(labels);
% % % % % % % % Nlobes = length(sts);
% % % % % % % % for j = 1:Nlobes
% % % % % % % %     ind = find(labels == sts(j));
% % % % % % % %     T = tempMat(:,(ind-1)*20+1);
% % % % % % % %     lmaxDepth        = max(T,[],2);
% % % % % % % %     T = tempMat(:,(ind-1)*20+2);
% % % % % % % %     lminDepth        = min(T,[],2);
% % % % % % % %     T = tempMat(:,(ind-1)*20+3);
% % % % % % % %     lmeanDepth       = sum(T,2)./(sum(logical(T),2) + eps);
% % % % % % % %     ltopLength       = sum(tempMat(:,(ind-1)*20+4),2);
% % % % % % % %     LbottomLength    = sum(tempMat(:,(ind-1)*20+5),2);
% % % % % % % %     LinterLength      = sum(tempMat(:,(ind-1)*20+6),2);
% % % % % % % %     lmeanLength      = sum(tempMat(:,(ind-1)*20+7),2);
% % % % % % % %     T = tempMat(:,(ind-1)*20+8);
% % % % % % % %     LmaxSPAM         = max(T,[],2);
% % % % % % % %     T = tempMat(:,(ind-1)*20+9);
% % % % % % % %     LminSPAM         = min(T,[],2);
% % % % % % % %     T = tempMat(:,(ind-1)*20+10);
% % % % % % % %     LmeanSPAM        = sum(T,2)./(sum(logical(T),2) + eps);
% % % % % % % %     
% % % % % % % %     T = tempMat(:,(ind-1)*20+11);
% % % % % % % %     RmaxDepth        = max(T,[],2);;
% % % % % % % %     T = tempMat(:,(ind-1)*20+12);
% % % % % % % %     RminDepth        = sum(T,2)./(sum(logical(T),2) + eps);
% % % % % % % %     T = tempMat(:,(ind-1)*20+13);
% % % % % % % %     RmeanDepth       = sum(T,2)./(sum(logical(T),2) + eps);
% % % % % % % %     RtopLength       = sum(tempMat(:,(ind-1)*20+14),2);
% % % % % % % %     RbottomLength    = sum(tempMat(:,(ind-1)*20+15),2);
% % % % % % % %     RinterLength      = sum(tempMat(:,(ind-1)*20+16),2);
% % % % % % % %     RmeanLength      = sum(tempMat(:,(ind-1)*20+17),2);
% % % % % % % %     T = tempMat(:,(ind-1)*20+18);
% % % % % % % %     RmaxSPAM         = max(T,[],2);
% % % % % % % %     T = tempMat(:,(ind-1)*20+19);
% % % % % % % %     RminSPAM         = min(T,[],2);
% % % % % % % %     T = tempMat(:,(ind-1)*20+20);
% % % % % % % %     RmeanSPAM        = sum(T,2)./(sum(logical(T),2) + eps);
% % % % % % % %     
% % % % % % % %     if j == 1
% % % % % % % %         tempMatL = [lmaxDepth lminDepth lmeanDepth ltopLength LbottomLength LinterLength lmeanLength LmaxSPAM LminSPAM LmeanSPAM ...
% % % % % % % %                     RmaxDepth RminDepth RmeanDepth RtopLength RbottomLength RinterLength RmeanLength RmaxSPAM RminSPAM RmeanSPAM];
% % % % % % % %     else
% % % % % % % %         tempMatL = [tempMatL lmaxDepth lminDepth lmeanDepth ltopLength LbottomLength LinterLength lmeanLength LmaxSPAM LminSPAM LmeanSPAM ...
% % % % % % % %                              RmaxDepth RminDepth RmeanDepth RtopLength RbottomLength RinterLength RmeanLength RmaxSPAM RminSPAM RmeanSPAM];
% % % % % % % %     end
% % % % % % % % end
% % % % % % % % 
% % % % % % % % % --------------- Grouping into Hemispheres
% % % % % % % % 
% % % % % % % % T = tempMat(:,1:20:end);
% % % % % % % % lmaxDepth        = max(T,[],2);
% % % % % % % % T = tempMat(:,2:20:end);
% % % % % % % % lminDepth        = min(T,[],2);
% % % % % % % % T = tempMat(:,3:20:end);
% % % % % % % % lmeanDepth       = sum(T,2)./(sum(logical(T),2) + eps);
% % % % % % % % 
% % % % % % % % ltopLength       = sum(tempMat(:,4:20:end),2);
% % % % % % % % LbottomLength    = sum(tempMat(:,5:20:end),2);
% % % % % % % % LinteLength      = sum(tempMat(:,6:20:end),2);
% % % % % % % % lmeanLength      = sum(tempMat(:,7:20:end),2);
% % % % % % % % 
% % % % % % % % T = tempMat(:,8:20:end);
% % % % % % % % LmaxSPAM         = max(T,[],2);;
% % % % % % % % T = tempMat(:,9:20:end);
% % % % % % % % LminSPAM         = min(T,[],2);;
% % % % % % % % T = tempMat(:,10:20:end);
% % % % % % % % LmeanSPAM        = sum(T,2)./(sum(logical(T),2) + eps);
% % % % % % % % 
% % % % % % % % T = tempMat(:,11:20:end);
% % % % % % % % RmaxDepth        = max(T,[],2);;
% % % % % % % % T = tempMat(:,12:20:end);
% % % % % % % % RminDepth        = min(T,[],2);;
% % % % % % % % T = tempMat(:,13:20:end);
% % % % % % % % RmeanDepth       = sum(T,2)./(sum(logical(T),2) + eps);
% % % % % % % % 
% % % % % % % % RtopLength       = sum(tempMat(:,14:20:end),2);
% % % % % % % % RbottomLength    = sum(tempMat(:,15:20:end),2);
% % % % % % % % RinteLength      = sum(tempMat(:,16:20:end),2);
% % % % % % % % RmeanLength      = sum(tempMat(:,17:20:end),2);
% % % % % % % % 
% % % % % % % % T = tempMat(:,18:20:end);
% % % % % % % % RmaxSPAM         = max(T,[],2);;
% % % % % % % % T = tempMat(:,19:20:end);
% % % % % % % % RminSPAM         = min(T,[],2);
% % % % % % % % T = tempMat(:,20:20:end);
% % % % % % % % RmeanSPAM        = sum(T,2)./(sum(logical(T),2) + eps);
% % % % % % % % 
% % % % % % % % TotalMat = [lmaxDepth lminDepth lmeanDepth ltopLength LbottomLength LinteLength lmeanLength LmaxSPAM LminSPAM LmeanSPAM ...
% % % % % % % %             RmaxDepth RminDepth RmeanDepth RtopLength RbottomLength RinteLength RmeanLength RmaxSPAM RminSPAM RmeanSPAM...
% % % % % % % %             tempMatL  tempMat];
% % % % % % % % 
% % % % % % % % % Creating Column Headers
% % % % % % % % temp = [correctNodes repmat('(mm),*',[size(correctNodes,1) 1])];
% % % % % % % % temp = temp';
% % % % % % % % temp = temp(:)';
% % % % % % % % temp(isspace(temp)) = [];
% % % % % % % % a = strread(temp,'%s','delimiter','*');
% % % % % % % % b = reshape(a,[Nstruct/2 2])';
% % % % % % % % b = b(:);
% % % % % % % % c = {'Left-Hemisphere(mm),';'Right-Hemisphere(mm),';'Left-FrontalLobe(mm),';'Right-FrontalLobe(mm),';'Left-ParietalLobe(mm),';'Right-ParietalLobe(mm),';'Left-TemporalLobe(mm),';'Right-TemporalLobe(mm),';'Left-OccipitalLobe(mm),';'Right-OccipitalLobe(mm),';
% % % % % % % %     'Left-CentralSulcus(mm),';'Right-CentralSulcus(mm),';'Left-LateralFissure(mm),';'Right-LateralFissure(mm),';'Left-CingulateSulcus(mm),';'Right-CingulateSulcus(mm),'; 'Left-ParietoOccipitalFissure(mm),';'Right-ParietoOccipitalFissure(mm),' ;'Left--SubcallosalSulcus(mm),';'Right-SubcallosalSulcus(mm),'};
% % % % % % % % regNames = [c;b];
% % % % % % % % 
% % % % % % % % LHs = strvcat('geodDepthMax-','geodDepthMin-','geodDepthMean-','topLength-','bottomLength-','intercepLength-','lengthMean-','spamMax-','spamMin-','spamMean-');
% % % % % % % % temp = repmat(regNames,[1 size(LHs,1)])';
% % % % % % % % aaaa = [repmat(LHs,[size(temp(:),1)/size(LHs,1) 1]) char(temp(:))];
% % % % % % % % aba = aaaa';aba = aba(:)';aba(isspace(aba))=[];
% % % % % % % % colHeaders = strread(aba,'%s','delimiter',',');
% % % % % % % % 
% % % % % % % % % Creating Table
% % % % % % % % repComa = repmat(',',[Ns+1 1]);
% % % % % % % % for j = 1:length(colHeaders)
% % % % % % % %     if j == 1
% % % % % % % %         cad2print = [strvcat('SubjID',Ids) repComa strvcat(colHeaders{j},num2str(TotalMat(:,j)))];
% % % % % % % %     else
% % % % % % % %         cad2print = [cad2print repComa strvcat(colHeaders{j},num2str(TotalMat(:,j))) ];
% % % % % % % %     end
% % % % % % % % end
% % % % % % % % 
% % % % % % % % % -------- Saving General Stat File
% % % % % % % % fid = fopen(statTable,'wt')
% % % % % % % %     
% % % % % % % % for i = 1:size(cad2print,1)
% % % % % % % %     fprintf(fid,'%s\n',cad2print(i,:));
% % % % % % % %     
% % % % % % % %     if i > 1
% % % % % % % %         subjId = deblank(Ids(i-1,:));
% % % % % % % %         warning off;
% % % % % % % %         mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'stats']);
% % % % % % % %         IstatTable = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'stats' filesep subjId '_RII_sulcalStats.txt'];
% % % % % % % %         fidind = fopen(IstatTable,'wt');
% % % % % % % %         fprintf(fidind,'%s\n',cad2print(1,:));
% % % % % % % %         fprintf(fidind,'%s\n',cad2print(i-1,:));
% % % % % % % %         fclose(fidind);
% % % % % % % %     end
% % % % % % % % end
% % % % % % % % fclose(fid);
% % % % % % % % 
% % % % % % % % varargout{1} = statTable;
% % % % % % % % return;