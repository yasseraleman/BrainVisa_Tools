function varargout = BrainVisa_stats_to_table_2(varargin);
%
% Syntax :
%  statTable = BrainVisa_stats_to(BrainVisaDatabaseDir, IdFile, statTable);
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

IndivDir = '/media/COSAS/8-BrainVISADataBase-HCP/BrainVISA_Morpho_Metrics';
BrainVisaDatabaseDir =  '/media/COSAS/8-BrainVISADataBase-HCP';

outDir = '/media/COSAS/8-BrainVISADataBase-HCP/BrainVISA_Morpho_Metrics';
warning off;
mkdir(outDir);
statTable = [outDir filesep 'BioMarkers_Table-BrainVISA.txt'];
Nbvstatcols = 35;
%% ---  Detecting the number of subjects
temp = dir([IndivDir filesep '*morpho*']);
if ~isempty(temp)
    testFile = [IndivDir filesep deblank(temp(1).name)];
    % -------- Creating String Cad to read the stat file
    a = repmat('%s',[1 Nbvstatcols]);
    temps = repmat('Temp',[Nbvstatcols 1]);
    b = [temps num2str([1:Nbvstatcols]') repmat(',',[Nbvstatcols 1])]';
    b = b(:)';
    b(isspace(b)) = [];
    b(end) = [];
    cad = ['[' b '] = textread(''' testFile '''' ',' '''' a '''' ',' '''' 'delimiter'  '''' ',' '''' ' ' '''' ');'];
    % -------- E    nd of Creating String Cad to read the stat file
    eval(cad); % Evaluating the string Cad
    Ids = char(Temp1(2:end));
    Ns = size(Ids,1);
else
    return;
end


%% ===================== Stats for Individual Nodes ==================== %%
% ========================== Individual Useful Nodes ==================== %


[groupsIds, groupsNamesL, correctNodesIL, labelsIdsL] = Detecting_BrainVisa_Groups('/media/COSAS/scripts/BrainVisa_myTools/Brainvisa_nomenclature_sulci+STS.txt', 'left');
colHeadersIL = correctNodesIL;
[groupsIds, groupsNamesR, correctNodesIR, labelsIdsR] = Detecting_BrainVisa_Groups('/media/COSAS/scripts/BrainVisa_myTools/Brainvisa_nomenclature_sulci+STS.txt', 'right');
colHeadersIR = correctNodesIR;

% % [correctNodesL, colHeadersL] = Detecting_BrainVisa_Groups('/media/COSAS/scripts/BrainVisa_myTools/Brainvisa_nomenclature_sulci.txt', 'left');
% % [correctNodesR, colHeadersR] = Detecting_BrainVisa_Groups('/media/COSAS/scripts/BrainVisa_myTools/Brainvisa_nomenclature_sulci.txt', 'right');

% ---- Columns separator
repComa = repmat(',',[Ns+1 1]);


Nstruct = size(correctNodesIL,1);
cont = 0;

measNamesL = '';
measNamesR = '';

for j = 1:Nstruct
    
    % ----------------------  Left Hemisphere --------------------------- %
    SulcNameL = deblank(correctNodesIL(j,:));
    
    Temp23 = num2str(zeros(Ns+1,1));
    Temp24 = num2str(zeros(Ns+1,1));
    Temp25 = num2str(zeros(Ns+1,1));
    Temp31 = num2str(zeros(Ns+1,1));
    Temp33 = num2str(zeros(Ns+1,1));
    
    filename = [IndivDir filesep  'morpho_' SulcNameL '.dat']; % Single Stat FIle
    if exist(filename, 'file')
        cont = cont + 1;
        % -------- Creating String Cad to read the stat file
        a = repmat('%s',[1 Nbvstatcols]);
        temps = repmat('Temp',[Nbvstatcols 1]);
        b = [temps num2str([1:Nbvstatcols]') repmat(',',[Nbvstatcols 1])]';
        b = b(:)';
        b(isspace(b)) = [];
        b(end) = [];
        cad = ['[' b '] = textread(''' filename '''' ',' '''' a '''' ',' '''' 'delimiter'  '''' ',' '''' ' ' '''' ');'];
        % -------- End of Creating String Cad to read the stat file
        eval(cad); % Evaluating the string Cad
    end
        
    lhmaxDepth      = strvcat(['geodDepthMax-'    SulcNameL '(mm)'],char(Temp23(2:end)));
    lhminDepth      = strvcat(['geodDepthMin-'    SulcNameL '(mm)'],char(Temp24(2:end)));
    lhmeanDepth     = strvcat(['geodDepthMean-'   SulcNameL '(mm)'],char(Temp25(2:end)));
    
    
    lhtopLength      = strvcat(['topLength-'      SulcNameL '(mm)'],repmat('0',[Ns 1]));
    lhbotLength      = strvcat(['bottomLength-'   SulcNameL '(mm)'],repmat('0',[Ns 1]));
    lhintercLength   = strvcat(['intercepLength-' SulcNameL '(mm)'],char(Temp31(2:end)));
    lhmeanLength     = strvcat(['lengthMean-'     SulcNameL '(mm)'],repmat('0',[Ns 1]));
    
    lhmaxSpam        = strvcat(['spamMax-'        SulcNameL '(mm)'],repmat('0',[Ns 1]));
    lhminSpam        = strvcat(['spamMin-'        SulcNameL '(mm)'],repmat('0',[Ns 1]));
    lhmeanSpam       = strvcat(['spamMean-'       SulcNameL '(mm)'],char(Temp33(2:end)));
    
    
    measNamesL = strvcat(measNamesL,['geodDepthMax-'    SulcNameL '(mm)'],['geodDepthMin-'    SulcNameL '(mm)'],['geodDepthMean-'   SulcNameL '(mm)'],...
        ['topLength-'      SulcNameL '(mm)'],['bottomLength-'   SulcNameL '(mm)'],['intercepLength-' SulcNameL '(mm)'],['lengthMean-'     SulcNameL '(mm)'],...
        ['spamMax-'        SulcNameL '(mm)'],['spamMin-'        SulcNameL '(mm)'],['spamMean-'       SulcNameL '(mm)']);
    
    % --------------------  End of Left Hemisphere ---------------------- %
    
    
    % ---------------------  Right Hemisphere --------------------------- %
    SulcNameR = deblank(correctNodesIR(j,:));
    Temp23 = num2str(zeros(Ns+1,1));
    Temp24 = num2str(zeros(Ns+1,1));
    Temp25 = num2str(zeros(Ns+1,1));
    Temp31 = num2str(zeros(Ns+1,1));
    Temp33 = num2str(zeros(Ns+1,1));
    
    filename = [IndivDir filesep  'morpho_' SulcNameR '.dat']; % Single Stat FIle
    if exist(filename, 'file')
        cont = cont + 1;
        % -------- Creating String Cad to read the stat file
        a = repmat('%s',[1 Nbvstatcols]);
        temps = repmat('Temp',[Nbvstatcols 1]);
        b = [temps num2str([1:Nbvstatcols]') repmat(',',[Nbvstatcols 1])]';
        b = b(:)';
        b(isspace(b)) = [];
        b(end) = [];
        cad = ['[' b '] = textread(''' filename '''' ',' '''' a '''' ',' '''' 'delimiter'  '''' ',' '''' ' ' '''' ');'];
        % -------- End of Creating String Cad to read the stat file
        eval(cad); % Evaluating the string Cad
    end
    
    
    % -------- Creacting Columns
    rhmaxDepth      = strvcat(['geodDepthMax-'      SulcNameR '(mm)'],char(Temp23(2:end)));
    rhminDepth      = strvcat(['geodDepthMin-'      SulcNameR '(mm)'],char(Temp24(2:end)));
    rhmeanDepth     = strvcat(['geodDepthMean-'     SulcNameR '(mm)'],char(Temp25(2:end)));
    
    
    rhtopLength      = strvcat(['topLength-'        SulcNameR '(mm)'],repmat('0',[Ns 1]));
    rhbotLength      = strvcat(['bottomLength-'     SulcNameR '(mm)'],repmat('0',[Ns 1]));
    rhintercLength   = strvcat(['intercepLength-'   SulcNameR '(mm)'],char(Temp31(2:end)));
    rhmeanLength     = strvcat(['lengthMean-'       SulcNameR '(mm)'],repmat('0',[Ns 1]));
    
    rhmaxSpam        = strvcat(['spamMax-'          SulcNameR '(mm)'],repmat('0',[Ns 1]));
    rhminSpam        = strvcat(['spamMin-'          SulcNameR '(mm)'],repmat('0',[Ns 1]));
    rhmeanSpam       = strvcat(['spamMean-'         SulcNameR '(mm)'],char(Temp33(2:end)));
    
    
    measNamesR = strvcat(measNamesR,['geodDepthMax-'    SulcNameR '(mm)'],['geodDepthMin-'    SulcNameR '(mm)'],['geodDepthMean-'   SulcNameR '(mm)'],...
        ['topLength-'      SulcNameR '(mm)'],['bottomLength-'   SulcNameR '(mm)'],['intercepLength-' SulcNameR '(mm)'],['lengthMean-'     SulcNameR '(mm)'],...
        ['spamMax-'        SulcNameR '(mm)'],['spamMin-'        SulcNameR '(mm)'],['spamMean-'       SulcNameR '(mm)']);
    
    % --------------------  End of Right Hemisphere --------------------- %
    
    % Joinning Values
    if j == 1
        Indiv2print = [lhmaxDepth repComa lhminDepth repComa lhmeanDepth repComa lhtopLength repComa lhbotLength repComa lhintercLength repComa lhmeanLength repComa lhmaxSpam repComa lhminSpam repComa lhmeanSpam repComa...
            rhmaxDepth repComa rhminDepth repComa rhmeanDepth repComa rhtopLength repComa rhbotLength repComa rhintercLength repComa rhmeanLength repComa rhmaxSpam repComa rhminSpam repComa rhmeanSpam];
        
        metricsMatL = [str2num(lhmaxDepth(2:end,:)) str2num(lhminDepth(2:end,:)) str2num(lhmeanDepth(2:end,:)) str2num(lhtopLength(2:end,:)) str2num(lhbotLength(2:end,:)) str2num(lhintercLength(2:end,:)) str2num(lhmeanLength(2:end,:)) str2num(lhmaxSpam(2:end,:)) str2num(lhminSpam(2:end,:)) str2num(lhmeanSpam(2:end,:))];
        metricsMatR = [str2num(rhmaxDepth(2:end,:)) str2num(rhminDepth(2:end,:)) str2num(rhmeanDepth(2:end,:)) str2num(rhtopLength(2:end,:)) str2num(rhbotLength(2:end,:)) str2num(rhintercLength(2:end,:)) str2num(rhmeanLength(2:end,:)) str2num(rhmaxSpam(2:end,:)) str2num(rhminSpam(2:end,:)) str2num(rhmeanSpam(2:end,:))];
    else
        Indiv2print = [Indiv2print repComa lhmaxDepth repComa lhminDepth repComa lhmeanDepth repComa lhtopLength repComa lhbotLength repComa lhintercLength repComa lhmeanLength repComa lhmaxSpam repComa lhminSpam repComa lhmeanSpam ...
            repComa rhmaxDepth repComa rhminDepth repComa rhmeanDepth repComa rhtopLength repComa rhbotLength repComa rhintercLength repComa rhmeanLength repComa rhmaxSpam repComa rhminSpam repComa rhmeanSpam];
        
        metricsMatL = [metricsMatL str2num(lhmaxDepth(2:end,:)) str2num(lhminDepth(2:end,:)) str2num(lhmeanDepth(2:end,:)) str2num(lhtopLength(2:end,:)) str2num(lhbotLength(2:end,:)) str2num(lhintercLength(2:end,:)) str2num(lhmeanLength(2:end,:)) str2num(lhmaxSpam(2:end,:)) str2num(lhminSpam(2:end,:)) str2num(lhmeanSpam(2:end,:))];
        metricsMatR = [metricsMatR str2num(rhmaxDepth(2:end,:)) str2num(rhminDepth(2:end,:)) str2num(rhmeanDepth(2:end,:)) str2num(rhtopLength(2:end,:)) str2num(rhbotLength(2:end,:)) str2num(rhintercLength(2:end,:)) str2num(rhmeanLength(2:end,:)) str2num(rhmaxSpam(2:end,:)) str2num(rhminSpam(2:end,:)) str2num(rhmeanSpam(2:end,:))];
    end
end

%% ================ End of Stats for Individual Nodes ================== %%

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

for j = 1:length(Ids)
    
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
    mkdir([ BrainVisaDatabaseDir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'BrainVISA' filesep 'stats']);
    warning on;
    IstatTable = [ BrainVisaDatabaseDir filesep 'subjects' filesep opts.subjid filesep 't1mri' filesep 'default_acquisition' filesep 'BrainVISA' filesep 'stats' filesep opts.subjid '_sulcal_auto_Stats.txt'];
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















% % % % 
% % % % 
% % % % 
% % % % 
% % % % temp = dir([LobarDir filesep '*morpho*']);
% % % % if ~isempty(temp)
% % % %     testFile = [LobarDir filesep deblank(temp(1).name)];
% % % %     % -------- Creating String Cad to read the stat file
% % % %     a = repmat('%s',[1 34]);
% % % %     temps = repmat('Temp',[34 1]);
% % % %     b = [temps num2str([1:34]') repmat(',',[34 1])]';
% % % %     b = b(:)';
% % % %     b(isspace(b)) = [];
% % % %     b(end) = [];
% % % %     cad = ['[' b '] = textread(''' testFile '''' ',' '''' a '''' ',' '''' 'delimiter'  '''' ',' '''' ' ' '''' ');'];
% % % %     % -------- E    nd of Creating String Cad to read the stat file
% % % %     eval(cad); % Evaluating the string Cad
% % % %     Ids = char(Temp1(2:end));
% % % %     Ns = size(Ids,1);
% % % % else
% % % %     return;
% % % % end
% % % % 
% % % % repComa = repmat(',',[Ns+1 1]);
% % % % 
% % % % 
% % % % Nstruct = size(correctNodesL,1);
% % % % cont = 0;
% % % % for j = 1:Nstruct
% % % %     
% % % %     % ----------------------  Left Hemisphere --------------------------- %
% % % %     SulcNameL = deblank(correctNodesL(j,:));
% % % %     
% % % %     Temp23 = num2str(zeros(Ns+1,1));
% % % %     Temp24 = num2str(zeros(Ns+1,1));
% % % %     Temp25 = num2str(zeros(Ns+1,1));
% % % %     Temp31 = num2str(zeros(Ns+1,1));
% % % %     Temp33 = num2str(zeros(Ns+1,1));
% % % %     
% % % %     filename = [LobarDir filesep  'morpho_' SulcNameL '.dat']; % Single Stat FIle
% % % %     if exist(filename, 'file')
% % % %         cont = cont + 1;
% % % %         % -------- Creating String Cad to read the stat file
% % % %         a = repmat('%s',[1 34]);
% % % %         temps = repmat('Temp',[34 1]);
% % % %         b = [temps num2str([1:34]') repmat(',',[34 1])]';
% % % %         b = b(:)';
% % % %         b(isspace(b)) = [];
% % % %         b(end) = [];
% % % %         cad = ['[' b '] = textread(''' filename '''' ',' '''' a '''' ',' '''' 'delimiter'  '''' ',' '''' ' ' '''' ');'];
% % % %         % -------- End of Creating String Cad to read the stat file
% % % %         eval(cad); % Evaluating the string Cad
% % % %     end
% % % %     % -------- Creacting Columns
% % % %     if j == 1
% % % %         lhmaxdepth      =  [str2num(char(Temp23(2:end)))];
% % % %         lhmindepth      =  [str2num(char(Temp24(2:end)))];
% % % %         lhmeandepth     =  [str2num(char(Temp25(2:end)))];
% % % %         lhinterclength  =  [str2num(char(Temp31(2:end)))];
% % % %         lhmeanspam      =  [str2num(char(Temp33(2:end)))];
% % % %         
% % % %     else
% % % %         lhmaxdepth      =  [lhmaxdepth      str2num(char(Temp23(2:end)))];
% % % %         lhmindepth      =  [lhmindepth      str2num(char(Temp24(2:end)))];
% % % %         lhmeandepth     =  [lhmeandepth     str2num(char(Temp25(2:end)))];
% % % %         lhinterclength  =  [lhinterclength  str2num(char(Temp31(2:end)))];
% % % %         lhmeanspam      =  [lhmeanspam      str2num(char(Temp33(2:end)))];
% % % %     end
% % % %     
% % % %     
% % % %     % -------- Creacting Columns
% % % %     
% % % %     NameLabel = deblank(colHeadersL(j,:));
% % % %     lhmaxDepth      = strvcat(['geodDepthMax-'   NameLabel '(mm)'],char(Temp23(2:end)));
% % % %     lhminDepth      = strvcat(['geodDepthMin-'   NameLabel '(mm)'],char(Temp24(2:end)));
% % % %     lhmeanDepth     = strvcat(['geodDepthMean-'  NameLabel '(mm)'],char(Temp25(2:end)));
% % % %     
% % % %     
% % % %     lhtopLength      = strvcat(['topLength-'      NameLabel '(mm)'],repmat('0',[Ns 1]));
% % % %     lhbotLength      = strvcat(['bottomLength-'   NameLabel '(mm)'],repmat('0',[Ns 1]));
% % % %     lhintercLength   = strvcat(['intercepLength-' NameLabel '(mm)'],char(Temp31(2:end)));
% % % %     lhmeanLength     = strvcat(['lengthMean-'     NameLabel '(mm)'],repmat('0',[Ns 1]));
% % % %     
% % % %     lhmaxSpam        = strvcat(['spamMax-'        NameLabel '(mm)'],repmat('0',[Ns 1]));
% % % %     lhminSpam        = strvcat(['spamMin-'        NameLabel '(mm)'],repmat('0',[Ns 1]));
% % % %     lhmeanSpam       = strvcat(['spamMean-'       NameLabel '(mm)'],char(Temp33(2:end)));
% % % %     % --------------------  End of Left Hemisphere ---------------------- %
% % % %     
% % % %     
% % % %     % ---------------------  Right Hemisphere --------------------------- %
% % % %     SulcNameR = deblank(correctNodesR(j,:));
% % % %     Temp23 = num2str(zeros(Ns+1,1));
% % % %     Temp24 = num2str(zeros(Ns+1,1));
% % % %     Temp25 = num2str(zeros(Ns+1,1));
% % % %     Temp31 = num2str(zeros(Ns+1,1));
% % % %     Temp33 = num2str(zeros(Ns+1,1));
% % % %     
% % % %     filename = [LobarDir filesep  'morpho_' SulcNameR '.dat']; % Single Stat FIle
% % % %     if exist(filename, 'file')
% % % %         cont = cont + 1;
% % % %         % -------- Creating String Cad to read the stat file
% % % %         a = repmat('%s',[1 34]);
% % % %         temps = repmat('Temp',[34 1]);
% % % %         b = [temps num2str([1:34]') repmat(',',[34 1])]';
% % % %         b = b(:)';
% % % %         b(isspace(b)) = [];
% % % %         b(end) = [];
% % % %         cad = ['[' b '] = textread(''' filename '''' ',' '''' a '''' ',' '''' 'delimiter'  '''' ',' '''' ' ' '''' ');'];
% % % %         % -------- End of Creating String Cad to read the stat file
% % % %         eval(cad); % Evaluating the string Cad
% % % %     end
% % % %     
% % % %     % -------- Creacting Columns
% % % %     if j == 1
% % % %         rhmaxdepth      =  [str2num(char(Temp23(2:end)))];
% % % %         rhmindepth      =  [str2num(char(Temp24(2:end)))];
% % % %         rhmeandepth     =  [str2num(char(Temp25(2:end)))];
% % % %         rhinterclength  =  [str2num(char(Temp31(2:end)))];
% % % %         rhmeanspam      =  [str2num(char(Temp33(2:end)))];
% % % %         
% % % %     else
% % % %         rhmaxdepth      =  [rhmaxdepth      str2num(char(Temp23(2:end)))];
% % % %         rhmindepth      =  [rhmindepth      str2num(char(Temp24(2:end)))];
% % % %         rhmeandepth     =  [rhmeandepth     str2num(char(Temp25(2:end)))];
% % % %         rhinterclength  =  [rhinterclength  str2num(char(Temp31(2:end)))];
% % % %         rhmeanspam      =  [rhmeanspam      str2num(char(Temp33(2:end)))];
% % % %     end
% % % %     
% % % %     
% % % %     NameLabel = deblank(colHeadersR(j,:));
% % % %        
% % % %     % -------- Creacting Columns
% % % %     rhmaxDepth      = strvcat(['geodDepthMax-'     NameLabel '(mm)'],char(Temp23(2:end)));
% % % %     rhminDepth      = strvcat(['geodDepthMin-'     NameLabel '(mm)'],char(Temp24(2:end)));
% % % %     rhmeanDepth     = strvcat(['geodDepthMean-'    NameLabel '(mm)'],char(Temp25(2:end)));
% % % %     
% % % %         
% % % %     rhtopLength      = strvcat(['topLength-'        NameLabel '(mm)'],repmat('0',[Ns 1]));
% % % %     rhbotLength      = strvcat(['bottomLength-'     NameLabel '(mm)'],repmat('0',[Ns 1]));
% % % %     rhintercLength   = strvcat(['intercepLength-'   NameLabel '(mm)'],char(Temp31(2:end)));
% % % %     rhmeanLength     = strvcat(['lengthMean-'       NameLabel '(mm)'],repmat('0',[Ns 1]));
% % % %        
% % % %     rhmaxSpam        = strvcat(['spamMax-'          NameLabel '(mm)'],repmat('0',[Ns 1]));
% % % %     rhminSpam        = strvcat(['spamMin-'          NameLabel '(mm)'],repmat('0',[Ns 1]));
% % % %     rhmeanSpam       = strvcat(['spamMean-'         NameLabel '(mm)'],char(Temp33(2:end)));
% % % %     
% % % %     % --------------------  End of Right Hemisphere --------------------- %
% % % %     
% % % %     % --- Joinning Lobar Values
% % % %     if j == 1
% % % %         Lobe2print = [lhmaxDepth repComa lhminDepth repComa lhmeanDepth repComa lhtopLength repComa lhbotLength repComa lhintercLength repComa lhmeanLength repComa lhmaxSpam repComa lhminSpam repComa lhmeanSpam repComa...
% % % %                       rhmaxDepth repComa rhminDepth repComa rhmeanDepth repComa rhtopLength repComa rhbotLength repComa rhintercLength repComa rhmeanLength repComa rhmaxSpam repComa rhminSpam repComa rhmeanSpam];
% % % %     else
% % % %         Lobe2print = [Lobe2print repComa lhmaxDepth repComa lhminDepth repComa lhmeanDepth repComa lhtopLength repComa lhbotLength repComa lhintercLength repComa lhmeanLength repComa lhmaxSpam repComa lhminSpam repComa lhmeanSpam ...
% % % %                                  repComa rhmaxDepth repComa rhminDepth repComa rhmeanDepth repComa rhtopLength repComa rhbotLength repComa rhintercLength repComa rhmeanLength repComa rhmaxSpam repComa rhminSpam repComa rhmeanSpam];
% % % %     end
% % % % end
% % % % %%  ===================== End Grouping into Lobes ====================== %%
% % % % 
% % % % 
% % % % %%  =================== Grouping into Hemispheres ====================== %%
% % % % 
% % % % % -------------------------  Left Hemisphere ---------------------------- %
% % % % lhemimaxdepth    = max(lhmaxdepth,[],2);
% % % % lhmaxDepth       = strvcat(['geodDepthMax-Left-Hemisphere(mm)'],num2str(lhemimaxdepth));
% % % % lhemimindepth    = min(lhmindepth,[],2);
% % % % lhminDepth       = strvcat(['geodDepthMin-Left-Hemisphere(mm)'],num2str(lhemimindepth));
% % % % lhemimeandepth   = mean(lhmeandepth,2);
% % % % lhmeanDepth      = strvcat(['geodDepthMean-Left-Hemisphere(mm)'],num2str(lhemimeandepth));
% % % % 
% % % % lhtopLength      = strvcat(['topLength-Left-Hemisphere(mm)'],repmat('0',[Ns 1]));
% % % % lhbotLength      = strvcat(['bottomLength-Left-Hemisphere(mm)'],repmat('0',[Ns 1]));
% % % % lhemiintercdepth = sum(lhinterclength,2);
% % % % lhintercLength   = strvcat(['intercepLength-Left-Hemisphere(mm)'],num2str(lhemiintercdepth));
% % % % lhmeanLength     = strvcat(['lengthMean-Left-Hemisphere(mm)'],repmat('0',[Ns 1]));
% % % % 
% % % % lhmaxSpam        = strvcat(['spamMax-Left-Hemisphere(mm)'],repmat('0',[Ns 1]));
% % % % lhminSpam        = strvcat(['spamMin-Left-Hemisphere(mm)'],repmat('0',[Ns 1]));
% % % % lhemimeanspam    = mean(lhmeanspam,2);
% % % % lhmeanSpam       = strvcat(['spamMean-Left-Hemisphere(mm)'],num2str(lhemimeanspam));
% % % % % ------------------------  End of Left Hemisphere ---------------------- %
% % % %     
% % % %     
% % % % % -------------------------  Right Hemisphere --------------------------- %
% % % % rhemimaxdepth    = max(rhmaxdepth,[],2);
% % % % rhmaxDepth      = strvcat(['geodDepthMax-Right-Hemisphere(mm)'],num2str(rhemimaxdepth));
% % % % rhemimindepth    = min(rhmindepth,[],2);
% % % % rhminDepth      = strvcat(['geodDepthMin-Right-Hemisphere(mm)'],num2str(rhemimindepth));
% % % % rhemimeandepth   = mean(rhmeandepth,2);
% % % % rhmeanDepth     = strvcat(['geodDepthMean-Right-Hemisphere(mm)'],num2str(rhemimeandepth));
% % % % 
% % % % 
% % % % rhtopLength      = strvcat(['topLength-Right-Hemisphere(mm)'],repmat('0',[Ns 1]));
% % % % rhbotLength      = strvcat(['bottomLength-Right-Hemisphere(mm)'],repmat('0',[Ns 1]));
% % % % rhemiintercdepth = sum(rhinterclength,2);
% % % % rhintercLength   = strvcat(['intercepLength-Right-Hemisphere(mm)'],num2str(rhemiintercdepth));
% % % % rhmeanLength     = strvcat(['lengthMean-Right-Hemisphere(mm)'],repmat('0',[Ns 1]));
% % % % 
% % % % rhmaxSpam        = strvcat(['spamMax-Right-Hemisphere(mm)'],repmat('0',[Ns 1]));
% % % % rhminSpam        = strvcat(['spamMin-Right-Hemisphere(mm)'],repmat('0',[Ns 1]));
% % % % rhemimeanspam    = mean(rhmeanspam,2);
% % % % rhmeanSpam       = strvcat(['spamMean-Right-Hemisphere(mm)'],num2str(rhemimeanspam));
% % % % % --------------------  End of Left Hemisphere -------------------------- %
% % % %     
% % % %     
% % % % % -------------- Joinning Hemispheric Values Hemisphere ----------------- %
% % % % Hemi2print = [lhmaxDepth repComa lhminDepth repComa lhmeanDepth repComa lhtopLength repComa lhbotLength repComa lhintercLength repComa lhmeanLength repComa lhmaxSpam repComa lhminSpam repComa lhmeanSpam repComa...
% % % %               rhmaxDepth repComa rhminDepth repComa rhmeanDepth repComa rhtopLength repComa rhbotLength repComa rhintercLength repComa rhmeanLength repComa rhmaxSpam repComa rhminSpam repComa rhmeanSpam];
% % % % %%  =================== End of Grouping into Hemispheres =============== %%
% % % % 
% % % % % --------------- Global Stats Matrix 
% % % % cad2print = [strvcat('SubjID',Ids) repComa Hemi2print repComa Lobe2print repComa Indiv2print];


% % % % % -------- Saving General Stat File
% % % % if exist('statTable','var')
% % % %     [pth, nm, ext ] = fileparts(statTable);
% % % %     warning off;
% % % %     mkdir(pth);
% % % %     fid = fopen(statTable,'wt')
% % % %     
% % % %     for i = 1:size(cad2print,1)
% % % %         fprintf(fid,'%s\n',cad2print(i,:));
% % % %         
% % % %         if i > 1
% % % %             subjId = deblank(Ids(i-1,:));
% % % %             mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'BrainVISA' filesep 'stats']);
% % % %             IstatTable = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'BrainVISA' filesep 'stats' filesep subjId '_BrainVISA_sulcalStats.txt'];
% % % %             fidind = fopen(IstatTable,'wt');
% % % %             fprintf(fidind,'%s\n',cad2print(1,:));
% % % %             fprintf(fidind,'%s\n',cad2print(i-1,:));
% % % %             fclose(fidind);
% % % %         end
% % % %     end
% % % %     fclose(fid);
% % % % else
% % % %     for i = 2:size(cad2print,1)
% % % %             subjId = deblank(Ids(i-1,:));
% % % %             mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'BrainVISA' filesep 'stats']);
% % % %             IstatTable = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'BrainVISA' filesep 'stats' filesep subjId '_BrainVISA_sulcalStats.txt'];
% % % %             fidind = fopen(IstatTable,'wt');
% % % %             fprintf(fidind,'%s\n',cad2print(1,:));
% % % %             fprintf(fidind,'%s\n',cad2print(i,:));
% % % %             fclose(fidind);
% % % %     end
% % % %     statTable = '';
% % % % end
% % % % varargout{1} = statTable;
return;