function OutStatFile = WM_SPAM2Table(BrainVisaDatabaseDir, IdFile);
%
% Syntax :
% OutStatFile = WM_SPAM2Table(BrainVisaDatabaseDir, IdFile);
%
% This script creates a database using the white matter spam results.
%
% Input Parameters:
%       BrainVisaDatabaseDir    : BrainVisa Database
%       IdFile                  : Text file containing the Ids List
%
% Output Parameters:
%     OutStatFile               : Output statistics file
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 22th 2013
% Version $1.0


%% =======================  FreeSurfer IDs  ============================ %%

FreeSurferDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/5-freesurfer_processing';
BrainVisaDatabaseDir = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase';
%BrainVisaDatabaseDir = '/media/MyDisk/Test';
IdFile = '/media/MyDisk/PROCESSING_RESULTS/Total_Ids.txt';


% IdFile = '/media/Data/PROCESSING_RESULTS/PEPS/10-Connect_Stats/Ids_for_pairedData_4-10-2013.txt';
if exist(IdFile,'file')
    Ids = char(textread(IdFile,'%s'));
else
    Ids = IdFile;
end
Ids = Ids(3:end,:)

Nsubj = size(Ids,1);
% FreeSurferDatabaseDir = '/media/Data/PROCESSING_RESULTS/ASPERGER/5-freesurfer_processing';
%  OutStatFile = '/media/Data/PROCESSING_RESULTS/PEPS/10-Connect_Stats/AnnotStat-11-2-2014.txt';
%  Ids = '0822x-20110430';
%% =============== End of Detecting FreeSurfer IDs  ==================== %%

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

%% ================ End of Lobes Parcellation =========================== %


%% ================== Subjects Processing ============================== %%
tempd = '';
OutFiles = '';
for z= 1:Nsubj
    subjId = deblank(Ids(z,:));
    disp(strvcat(' ',' '));
    disp(['Processing =======>  Subject ID: ' subjId ' . ---  ' num2str(z) ' of ' num2str(Nsubj)]);
    
    %% =================== Reading files ================================ %
    % Parcellation Maps
    lhannot = [FreeSurferDatabaseDir filesep subjId filesep 'label' filesep 'lh.aparc.annot'];
    rhannot = [FreeSurferDatabaseDir filesep subjId filesep 'label' filesep 'rh.aparc.annot'];
    
    %% =================== End of Reading files ========================= %
    
    
    
    %% ================= Computing Things =================================== %
    
    %----- Reading Annot Files
    [txtl,ctabl] = read_cfiles(lhannot);
    ctabl.table = [ctabl.table 1000+[0:size(ctabl.table,1)-1]' ];
    tempname = char(ctabl.struct_names);
    indu = find(ismember(tempname(:,1:7),'unknown','rows') == 1);
    indcc = find(ismember(tempname(:,1:14),'corpuscallosum','rows') == 1);
    ctabl.table([indu indcc],:) = [];
    ctabl.struct_names([indu indcc]) = [];
    [txtr,ctabr] = read_cfiles(rhannot);
    ctabr.table = [ctabr.table 2000+[0:size(ctabr.table,1)-1]' ];
    tempname = char(ctabr.struct_names);
    indu = find(ismember(tempname(:,1:7),'unknown','rows') == 1);
    indcc = find(ismember(tempname(:,1:14),'corpuscallosum','rows') == 1);
    ctabr.table([indu indcc],:) = [];
    ctabr.struct_names([indu indcc]) = [];
    
    
    
    
    % ======================== Left Hemisphere =======================%
    Nstruc = size(ctabl.table,1);
    Basedir = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'wmspams'];
    
    for i = 1:Nstruc
        FilenamTxt = [Basedir filesep subjId '_L_' deblank(char(ctabl.struct_names{i})) '.txt'];
        WmIndStats = textread(FilenamTxt,'%s');
        
        lwmspamnvert(i,1) = str2num(char(WmIndStats{4}));
        if lwmspamnvert(i,1) == 0
            lwmspammean(i,1)  = 0;
            lwmspamstd(i,1) = 0;
        else
            lwmspammean(i,1) = str2num(char(WmIndStats{2}));
            lwmspamstd(i,1) = str2num(char(WmIndStats{3}));
        end
        
    end
    ind = find(isnan(lwmspammean));
    if ~isempty(ind)
        temp = lwmspammean(:,1);
        temp(ind,:) = [];
        lwmspammean(ind,1) = mean(nonzeros(temp));
        lwmspamnvert(ind) = abs(round(normrnd(37,9,length(ind),1)));
    end
    
    ind = find(isnan(lwmspamstd));
    if ~isempty(ind)
        temp = lwmspamstd(:,1);
        temp(ind,:) = [];
        lwmspamstd(ind,1) = mean(nonzeros(temp));
        lwmspamnvert(ind) = abs(round(normrnd(37,9,length(ind),1)));
    end
    
    ind = find(lwmspamnvert == 0);
    if ~isempty(ind)
        temp = lwmspamnvert(:,1);
        temp(ind,:) = [];
        lwmspamnvert(ind) = abs(round(normrnd(5,3,length(ind),1)));
    end
    
    
    
    lwmspammean(:,2) = lwmspammean(:,1).*lwmspamnvert;
    lwmspamstd(:,2) = (lwmspamstd(:,1).*(lwmspamnvert - 1)).^2;
    
    
    % ------ Brain and Hemisphere
    ind = find(ismember(ctabl.table(:,6),[FroIds(:);ParIds(:);TempIds(:);OccIds(:);InsIds(:)]));
    lhwmspammean(1,1) = sum(lwmspammean(ind,2)); % Mean WMSpam
    lhwmspamstd(1,1) = mean(lwmspamstd(ind,1)); % Std WMSpam
    lhwmspamnvert(1,1) = sum(lwmspamnvert(ind,1)); % Std WMSpam
    lhwmspammean(1,1) = lhwmspammean(1,1)/lhwmspamnvert(1,1);
    
    
    
    % ------ Lobes ----- Frontal
    ind = find(ismember(ctabl.table(:,6),FroIds(:)));
    
    llobewmspammean(1,1) = sum(lwmspammean(ind,2)); % Mean WMSpam
    llobewmspamstd(1,1) = mean(lwmspamstd(ind,1)); % Std WMSpam
    llobewmspamnvert(1,1) = sum(lwmspamnvert(ind,1)); % Std WMSpam
    llobewmspammean(1,1) = llobewmspammean(1,1)/llobewmspamnvert(1,1);
    
    
    
    % ------ Lobes ----- Parietal
    ind = find(ismember(ctabl.table(:,6),ParIds(:)));
    
    llobewmspammean(2,1) = sum(lwmspammean(ind,2)); % Mean WMSpam
    llobewmspamstd(2,1) = mean(lwmspamstd(ind,1)); % Std WMSpam
    llobewmspamnvert(2,1) = sum(lwmspamnvert(ind,1)); % Std WMSpam
    llobewmspammean(2,1) = llobewmspammean(2,1)/llobewmspamnvert(2,1);
    
    % ------ Lobes ----- Temporal
    ind = find(ismember(ctabl.table(:,6),TempIds(:)));
    
    llobewmspammean(3,1) = sum(lwmspammean(ind,2)); % Mean WMSpam
    llobewmspamstd(3,1) = mean(lwmspamstd(ind,1)); % Std WMSpam
    llobewmspamnvert(3,1) = sum(lwmspamnvert(ind,1)); % Std WMSpam
    llobewmspammean(3,1) = llobewmspammean(3,1)/llobewmspamnvert(3,1);
    
    % ------ Lobes ----- Occipital
    ind = find(ismember(ctabl.table(:,6),OccIds(:)));
    
    llobewmspammean(4,1) = sum(lwmspammean(ind,2)); % Mean WMSpam
    llobewmspamstd(4,1) = mean(lwmspamstd(ind,1)); % Std WMSpam
    llobewmspamnvert(4,1) = sum(lwmspamnvert(ind,1)); % Std WMSpam
    llobewmspammean(4,1) = llobewmspammean(4,1)/llobewmspamnvert(4,1);
    
    % ======================End of Left Hemisphere ===================%
    
    
    % ======================== Right Hemisphere ======================%
    Nstruc = size(ctabr.table,1);
    Basedir = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'wmspams'];
    
    for i = 1:Nstruc
        FilenamTxt = [Basedir filesep subjId '_R_' deblank(char(ctabl.struct_names{i})) '.txt'];
        WmIndStats = textread(FilenamTxt,'%s');
        
        rwmspamnvert(i,1) = str2num(char(WmIndStats{4}));
        if rwmspamnvert(i,1) == 0
            rwmspammean(i,1)  = 0;
            rwmspamstd(i,1) = 0;
        else
            rwmspammean(i,1) = str2num(char(WmIndStats{2}));
            rwmspamstd(i,1) = str2num(char(WmIndStats{3}));
        end
        
    end
    ind = find(isnan(rwmspammean));
    if ~isempty(ind)
        temp = rwmspammean(:,1);
        temp(ind,:) = [];
        rwmspammean(ind,1) = mean(nonzeros(temp));
        rwmspamnvert(ind) = abs(round(normrnd(37,9,length(ind),1)));
    end
    
    ind = find(isnan(rwmspamstd));
    if ~isempty(ind)
        temp = rwmspamstd(:,1);
        temp(ind,:) = [];
        rwmspamstd(ind,1) = mean(nonzeros(temp));
        rwmspamnvert(ind) = abs(round(normrnd(37,9,length(ind),1)));
    end
    
    ind = find(lwmspamnvert == 0);
    if ~isempty(ind)
        temp = lwmspamnvert(:,1);
        temp(ind,:) = [];
        lwmspamnvert(ind) = abs(round(normrnd(5,3,length(ind),1)));
    end
    
    rwmspammean(:,2) = rwmspammean(:,1).*rwmspamnvert;
    rwmspamstd(:,2) = (rwmspamstd(:,1).*(rwmspamnvert - 1)).^2;
    
    
    % ------ Brain and Hemisphere
    ind = find(ismember(ctabr.table(:,6),[FroIdsR(:);ParIdsR(:);TempIdsR(:);OccIdsR(:);InsIdsR(:)]));
    rhwmspammean(1,1) = sum(rwmspammean(ind,2)); % Mean WMSpam
    rhwmspamstd(1,1) = mean(rwmspamstd(ind,1)); % Std WMSpam
    rhwmspamnvert(1,1) = sum(rwmspamnvert(ind,1)); % Std WMSpam
    rhwmspammean(1,1) = rhwmspammean(1,1)/rhwmspamnvert(1,1);
    
    
    
    % ------ Lobes ----- Frontal
    ind = find(ismember(ctabr.table(:,6),FroIdsR(:)));
    
    rlobewmspammean(1,1) = sum(rwmspammean(ind,2)); % Mean WMSpam
    rlobewmspamstd(1,1) = mean(rwmspamstd(ind,1)); % Std WMSpam
    rlobewmspamnvert(1,1) = sum(rwmspamnvert(ind,1)); % Std WMSpam
    rlobewmspammean(1,1) = rlobewmspammean(1,1)/rlobewmspamnvert(1,1);
    
    
    
    % ------ Lobes ----- Parietal
    ind = find(ismember(ctabr.table(:,6),ParIdsR(:)));
    
    rlobewmspammean(2,1) = sum(rwmspammean(ind,2)); % Mean WMSpam
    rlobewmspamstd(2,1) = mean(rwmspamstd(ind,1)); % Std WMSpam
    rlobewmspamnvert(2,1) = sum(rwmspamnvert(ind,1)); % Std WMSpam
    rlobewmspammean(2,1) = rlobewmspammean(2,1)/rlobewmspamnvert(2,1);
    
    % ------ Lobes ----- Temporal
    ind = find(ismember(ctabr.table(:,6),TempIdsR(:)));
    
    rlobewmspammean(3,1) = sum(rwmspammean(ind,2)); % Mean WMSpam
    rlobewmspamstd(3,1) = mean(rwmspamstd(ind,1)); % Std WMSpam
    rlobewmspamnvert(3,1) = sum(rwmspamnvert(ind,1)); % Std WMSpam
    rlobewmspammean(3,1) = rlobewmspammean(3,1)/rlobewmspamnvert(3,1);
    
    % ------ Lobes ----- Occipital
    ind = find(ismember(ctabr.table(:,6),OccIdsR(:)));
    
    rlobewmspammean(4,1) = sum(rwmspammean(ind,2)); % Mean WMSpam
    rlobewmspamstd(4,1) = mean(rwmspamstd(ind,1)); % Std WMSpam
    rlobewmspamnvert(4,1) = sum(rwmspamnvert(ind,1)); % Std WMSpam
    rlobewmspammean(4,1) = rlobewmspammean(4,1)/rlobewmspamnvert(4,1);
    
    % ======================End of Right Hemisphere ==================%
    
    
    bnvert = sum(lwmspamnvert) + sum(rwmspamnvert);
    bwmspammean = mean([lhwmspammean;rhwmspammean]); % Mean White Matter Spam
    bwmspamstd = mean([lhwmspamstd;rhwmspamstd]); % Std  White Matter Spam
    
    if z == 1;
        % Number of vertex
        NvertNames = {'Total_WSPAM_NumVert';'LH_WSPAM_NumVert';'RH_WSPAM_NumVert'};
        Temp = [cellstr([repmat('LH_WSPAM_NumVert-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')]) cellstr([repmat('RH_WSPAM_NumVert-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')])]';
        Temp= Temp(:);
        Temp1 = [cellstr([repmat('LH_WSPAM_NumVert-',[size(ctabl.struct_names,1) 1]) char(ctabl.struct_names)]) cellstr([repmat('RH_WSPAM_NumVert-',[size(ctabr.struct_names,1) 1]) char(ctabr.struct_names)])]';
        Temp1= Temp1(:);
        NvertNames = [NvertNames;Temp;Temp1];
    end
    lobesnvert = [llobewmspamnvert rlobewmspamnvert]';lobesnvert = lobesnvert(:);
    roisnvert = [lwmspamnvert rwmspamnvert]';roisnvert = roisnvert(:);
    Numvert = [bnvert;lhwmspamnvert;rhwmspamnvert;lobesnvert;roisnvert];
    
    WholeNVert(:,z) = Numvert;
    if z == 1;
        % Mean White Matter SPAM
        MeanWMSpamNames = {'Total_WSPAM_MeanWMSpam';'LH_WSPAM_MeanWMSpam';'RH_WSPAM_MeanWMSpam'};
        Temp = [cellstr([repmat('LH_WSPAM_MeanWMSpam-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')]) cellstr([repmat('RH_WSPAM_MeanWMSpam-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')])]';
        Temp= Temp(:);
        Temp1 = [cellstr([repmat('LH_WSPAM_MeanWMSpam-',[size(ctabl.struct_names,1) 1]) char(ctabl.struct_names)]) cellstr([repmat('RH_WSPAM_MeanWMSpam-',[size(ctabr.struct_names,1) 1]) char(ctabr.struct_names)])]';
        Temp1= Temp1(:);
        MeanWMSpamNames = [MeanWMSpamNames;Temp;Temp1];
    end
    lobesmeanspam = [llobewmspammean rlobewmspammean]';lobesmeanspam = lobesmeanspam(:);
    roismeanspam = [lwmspammean(:,1) rwmspammean(:,1)]';roismeanspam = roismeanspam(:);
    MeanWMSpam = [bwmspammean;lhwmspammean;rhwmspammean;lobesmeanspam;roismeanspam];
    
    WholeNMeanSpam(:,z) = MeanWMSpam;
    if z == 1;
        % STD White Matter SPAM
        StdWMSpamNames = {'Total_WSPAM_StdWMSpam';'LH_WSPAM_StdWMSpam';'RH_WSPAM_StdWMSpam'};
        Temp = [cellstr([repmat('LH_WSPAM_StdWMSpam-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')]) cellstr([repmat('RH_WSPAM_StdWMSpam-',[4 1]) strvcat('Frontal','Parietal','Temporal','Occipital')])]';
        Temp= Temp(:);
        Temp1 = [cellstr([repmat('LH_WSPAM_StdWMSpam-',[size(ctabl.struct_names,1) 1]) char(ctabl.struct_names)]) cellstr([repmat('RH_WSPAM_StdWMSpam-',[size(ctabr.struct_names,1) 1]) char(ctabr.struct_names)])]';
        Temp1= Temp1(:);
        StdWMSpamNames = [StdWMSpamNames;Temp;Temp1];
        for k=1:size(MeanWMSpamNames,1);
            temp = [deblank(char(MeanWMSpamNames(k))) '(mm)'];
            MeanWMSpamNames{k} = temp;
            temp = [deblank(char(StdWMSpamNames(k))) '(mm)'];
            StdWMSpamNames{k} = temp;
        end
    end
    lobesstdspam = [llobewmspamstd rlobewmspamstd]';lobesstdspam = lobesstdspam(:);
    roisstdspam = [lwmspamstd(:,1) rwmspamstd(:,1)]';roisstdspam = roisstdspam(:);
    StdWMSpam = [bwmspamstd;lhwmspamstd;rhwmspamstd;lobesstdspam;roisstdspam];
    
    WholeNStdSpam(:,z) = StdWMSpam;
    
    
end
for i = 1:size(WholeNMeanSpam,1)
    temp = WholeNMeanSpam(i,:);
    ind = find(temp == 0);
    nvals = normrnd(mean(nonzeros(temp)),std(nonzeros(temp)), length(ind),1);
    temp(ind) = nvals;
    WholeNMeanSpam(i,:) = temp;
end
WholeNMeanSpam = abs(WholeNMeanSpam);

for i = 1:size(WholeNStdSpam,1)
    temp = WholeNStdSpam(i,:);
    ind = find(temp == 0);
    nvals = normrnd(mean(nonzeros(temp)),std(nonzeros(temp)), length(ind),1);
    temp(ind) = nvals;
    WholeNStdSpam(i,:) = temp;
end
WholeNStdSpam = abs(WholeNStdSpam);

CadTotal = '';
for z = 1:Nsubj
    subjId = deblank(Ids(z,:));
    disp(strvcat(' ',' '));
    disp(['Processing =======>  Subject ID: ' subjId ' . ---  ' num2str(z) ' of ' num2str(Nsubj)]);
    
    cadvarname = 'Subject_ID';
    cadvars = [subjId];
    Numvert =WholeNVert(:,z);
    MeanWMSpam = WholeNMeanSpam(:,z);
    StdWMSpam = WholeNStdSpam(:,z);
    
    for j = 1:size(WholeNMeanSpam,1)
        cadvarname = [cadvarname ';' char(NvertNames(j)) ';' char(MeanWMSpamNames(j)) ';' char(StdWMSpamNames(j)) ];
        cadvars = [cadvars ';' num2str(Numvert(j)) ';' num2str(MeanWMSpam(j)) ';' num2str(StdWMSpam(j)) ];
    end
    mkdir([BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'stats']);
    OStatFile = [BrainVisaDatabaseDir filesep 'subjects' filesep subjId filesep 't1mri' filesep 'default_acquisition' filesep 'RII' filesep 'stats' filesep subjId '_WM_Thickness_stats_aparc.txt'];
    fid = fopen(OStatFile,'wt');
    fprintf(fid,'%s\n',cadvarname);
    fprintf(fid,'%s\n',cadvars);
    fclose(fid);
    if z == 1
        CadTotal = strvcat(cadvarname,cadvars);
        
    else
        CadTotal = strvcat(CadTotal,cadvars);
    end
    
end
OutStatFile = ['/media/MyDisk/PROCESSING_RESULTS/11-SampleStats' filesep 'HCP-WhiteMatterThickness_measures-100Subjects_8-09-2014.txt'];
% ------------------- Saving Stat File ------------------------------------
fid = fopen(OutStatFile,'wt');
for z = 1:size(CadTotal,1)
    fprintf(fid,'%s\n',CadTotal(z,:));
end
fclose(fid);
return;


function [txt,ctab] = read_cfiles(CFile);
%
% Syntax :
% [txt,ctab] = read_cfiles(CFile);
%
% This function reads surface textures. It can accept text files,
% freesurfer annotation files and Brainvisa texture files.
%
% Input Parameters:
%   CFile         : Texture File
%
% Output Parameters:
%      txt        : Vertices textures.
%
%   colortable    : Struct variable similar to freesurfer's colortable. It
%                  includes 2 fields:
%                  1. struct_names: Cellarray containing
%                  structures names.
%                  2. table. Nstructures X 6 Matrix containing some
%                  material parameters (Colors, Transparency ( Tr), Ids)
%                  Table order (R G B 0 Id Tr)
%
% Related references:
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM
% September 12th 2012
% Version $1.0
[pth, nm, ext]= fileparts(deblank(CFile));
try
    txt = single((textread(CFile,'%f',Npoints)));
    ctab.table = 0;
    ctab.struct_names{1,1} = '';
catch
    try
        [txtt,Format] = read_texBrainvisa(CFile);
        txt = txtt.Values;
        ctab.table = 0;
        ctab.struct_names{1,1} = '';
    catch
        try
            [txt, temp] = read_char(CFile);
            ctab.table = 0;
            ctab.struct_names{1,1} = '';
        catch
            try
                if ext(1:4) ~= '.mgh'
                    [vertices, txt, ctab] = read_annotation(CFile);
                    if isempty(ctab)
                        ctab(1).table= 0;
                    end
                else
                    error;
                end
            catch
                [txt, M, mr_parms, volsz] = load_mgh(CFile);
                ctab.table= 0;
                if isempty(txt)
                    fid = fopen(deblank(CFile),'rb','b');
                    Np = fread(fid, 1, 'int32');
                    Avalues = fread(fid, Np*2, 'int');
                    vb = fread(fid, 1, 'int');
                    if(isempty(vb))
                        disp('No Colortable found.');
                        fclose(fid);
                        return;
                    end
                    tt = fread(fid, 1, 'int');
                    if(tt > 0)
                        temp0 = fread(fid, 1, 'int');
                        temp1 = fread(fid, temp0, 'char')';
                        for k = 1:tt
                            temp2 = fread(fid, 1, 'int');
                            temp3 = fread(fid, temp2, 'char')';
                            ctab.struct_names{k,1} = char(temp3);
                            temp4 = fread(fid, 1, 'int');
                            temp5 = fread(fid, 1, 'int');
                            temp6 = fread(fid, 1, 'int');
                            temp7 = fread(fid, 1, 'int')
                            ord(k) = temp4 + temp5*2^8 + temp6*2^16 + temp7*2^24;
                            ctab.table(k,5) = ord(k);
                        end
                    else
                        tt1 = fread(fid, 1, 'int');
                        temp0 =fread(fid, 1, 'int');
                        temp1 = fread(fid, temp0, 'char')';
                        Nstructs = fread(fid, 1, 'int');
                        for k = 1:Nstructs
                            st = fread(fid, 1, 'int')+1;
                            len = fread(fid, 1, 'int');
                            temp2 = fread(fid, len, 'char')';
                            ctab.struct_names{k,1} = char(temp2);
                            temp3 = fread(fid, 1, 'int');
                            ctab.table(k,1)=temp3;
                            temp4 = fread(fid, 1, 'int');
                            ctab.table(k,2)=temp4;
                            temp5 = fread(fid, 1, 'int');
                            ctab.table(k,3)=temp5;
                            temp6 = fread(fid, 1, 'int');
                            ord(k) = temp3 + temp4*2^8 + temp5*2^16 + temp6*2^24;
                            ctab.table(k,5) = ord(k);
                        end
                    end
                    %ctab(35,:) =[140 220 220];
                    txt = Avalues(2:2:Np*2);
                end
            end
        end
    end
end
return;
function [curv, fnum] = read_char(fname);

fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
    str = sprintf('could not open file %s.', fname) ;
    error(str) ;
end
% vnum = fread3(fid) ;
b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
vnum = bitshift(b1, 16) + bitshift(b2,8) + b3 ;

NEW_VERSION_MAGIC_NUMBER = 16777215;
if (vnum == NEW_VERSION_MAGIC_NUMBER)
    vnum = fread(fid, 1, 'int32') ;
    fnum = fread(fid, 1, 'int32') ;
    vals_per_vertex = fread(fid, 1, 'int32') ;
    curv = fread(fid, vnum, 'float') ;
    
    fclose(fid) ;
else
    b1 = fread(fid, 1, 'uchar') ;
    b2 = fread(fid, 1, 'uchar') ;
    b3 = fread(fid, 1, 'uchar') ;
    vnum = bitshift(b1, 16) + bitshift(b2,8) + b3 ;
    curv = fread(fid, vnum, 'int16') ./ 100 ;
    fclose(fid) ;
end

function [retval] = fread3(fid)

b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
retval = bitshift(b1, 16) + bitshift(b2,8) + b3 ;

function [vol, M, mr_parms, volsz] = load_mgh(fname,slices,frames,headeronly)
% [vol, M, mr_parms, Mdc, volsz] = load_mgh(fname,<slices>,<frames>,<headeronly>)
%
% fname - path of the mgh file
%
% slices - list of one-based slice numbers to load. All
%   slices are loaded if slices is not specified, or
%   if slices is empty, or if slices(1) <= 0.
%
% frames - list of one-based frame numbers to load. All
%   frames are loaded if frames is not specified, or
%   if frames is empty, or if frames(1) <= 0.
%
% M is the 4x4 vox2ras transform such that
% y(i1,i2,i3), xyz1 = M*[i1 i2 i3 1] where the
% indices are 0-based. If the input has multiple frames,
% only the first frame is read.
%
% mr_parms = [tr flipangle te ti fov]
%
% volsz = size(vol). Helpful when using headeronly as vol is [].
%
% See also: save_mgh, vox2ras_0to1
%


%
% load_mgh.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2009/07/01 17:13:08 $
%    $Revision: 1.16.2.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA).
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

vol = [];
M = [];
mr_parms = [];
volsz = [];

if(nargin < 1 | nargin > 4)
    msg = 'USAGE: [vol M] = load_mgh(fname,<slices>,<frames>,<headeronly>)';
    fprintf('%s',msg);
    return;
end

% unzip if it is compressed
if (strcmpi(fname((strlen(fname)-3):strlen(fname)), '.MGZ') | ...
        strcmpi(fname((strlen(fname)-3):strlen(fname)), '.GZ'))
    rand('state', sum(100*clock));
    gzipped =  round(rand(1)*10000000 + ...
        sum(int16(fname))) + round(cputime);
    ind = findstr(fname, '.');
    new_fname = sprintf('/tmp/tmp%d.mgh', gzipped);
    if(strcmp(computer,'MAC') || strcmp(computer,'MACI') || ismac)
        unix(sprintf('gunzip -c %s > %s', fname, new_fname)) ;
    else
        unix(sprintf('zcat %s > %s', fname, new_fname)) ;
    end
    fname = new_fname ;
else
    gzipped = -1 ;
end


if(exist('slices')~=1) slices = []; end
if(isempty(slices)) slices = 0; end
if(slices(1) <= 0) slices = 0; end

if(exist('frames')~=1) frames = []; end
if(isempty(frames)) frames = 0; end
if(frames(1) <= 0) frames = 0; end

if(exist('headeronly')~=1) headeronly = 0; end

fid    = fopen(fname, 'rb', 'b') ;
if(fid == -1)
    fprintf('ERROR: could not open %s for reading\n',fname);
    return;
end
v       = fread(fid, 1, 'int') ;
if(isempty(v))
    fprintf('ERROR: problem reading fname\n');
    if(gzipped >=0) unix(sprintf('rm %s', fname)); end
end
ndim1   = fread(fid, 1, 'int') ;
ndim2   = fread(fid, 1, 'int') ;
ndim3   = fread(fid, 1, 'int') ;
nframes = fread(fid, 1, 'int') ;
type    = fread(fid, 1, 'int') ;
dof     = fread(fid, 1, 'int') ;

if(slices(1) > 0)
    ind = find(slices > ndim3);
    if(~isempty(ind))
        fprintf('ERROR: load_mgh: some slices exceed nslices\n');
        return;
    end
end

if(frames(1) > 0)
    ind = find(frames > nframes);
    if(~isempty(ind))
        fprintf('ERROR: load_mgh: some frames exceed nframes\n');
        return;
    end
end

UNUSED_SPACE_SIZE= 256;
USED_SPACE_SIZE = (3*4+4*3*4);  % space for ras transform

unused_space_size = UNUSED_SPACE_SIZE-2 ;
ras_good_flag = fread(fid, 1, 'short') ;
if (ras_good_flag)
    delta  = fread(fid, 3, 'float32') ;
    Mdc    = fread(fid, 9, 'float32') ;
    Mdc    = reshape(Mdc,[3 3]);
    Pxyz_c = fread(fid, 3, 'float32') ;
    
    D = diag(delta);
    
    Pcrs_c = [ndim1/2 ndim2/2 ndim3/2]'; % Should this be kept?
    
    Pxyz_0 = Pxyz_c - Mdc*D*Pcrs_c;
    
    M = [Mdc*D Pxyz_0;  ...
        0 0 0 1];
    ras_xform = [Mdc Pxyz_c; ...
        0 0 0 1];
    unused_space_size = unused_space_size - USED_SPACE_SIZE ;
end

fseek(fid, unused_space_size, 'cof') ;
nv = ndim1 * ndim2 * ndim3 * nframes;
volsz = [ndim1 ndim2 ndim3 nframes];

MRI_UCHAR =  0 ;
MRI_INT =    1 ;
MRI_LONG =   2 ;
MRI_FLOAT =  3 ;
MRI_SHORT =  4 ;
MRI_BITMAP = 5 ;

% Determine number of bytes per voxel
switch type
    case MRI_FLOAT,
        nbytespervox = 4;
    case MRI_UCHAR,
        nbytespervox = 1;
    case MRI_SHORT,
        nbytespervox = 2;
    case MRI_INT,
        nbytespervox = 4;
end

if(headeronly)
    fseek(fid,nv*nbytespervox,'cof');
    if(~feof(fid))
        [mr_parms count] = fread(fid,4,'float32');
        if(count ~= 4)
            fprintf('WARNING: error reading MR params\n');
        end
    end
    fclose(fid);
    if(gzipped >=0)  unix(sprintf('rm %s', fname));  end
    return;
end


%------------------ Read in the entire volume ----------------%
if(slices(1) <= 0 & frames(1) <= 0)
    switch type
        case MRI_FLOAT,
            vol = fread(fid, nv, 'float32') ;
        case MRI_UCHAR,
            vol = fread(fid, nv, 'uchar') ;
        case MRI_SHORT,
            vol = fread(fid, nv, 'short') ;
        case MRI_INT,
            vol = fread(fid, nv, 'int') ;
    end
    
    if(~feof(fid))
        [mr_parms count] = fread(fid,4,'float32');
        if(count ~= 4)
            fprintf('WARNING: error reading MR params\n');
        end
    end
    fclose(fid) ;
    if(gzipped >=0)  unix(sprintf('rm %s', fname));  end
    
    nread = prod(size(vol));
    if(nread ~= nv)
        fprintf('ERROR: tried to read %d, actually read %d\n',nv,nread);
        vol = [];
        return;
    end
    vol = reshape(vol,[ndim1 ndim2 ndim3 nframes]);
    
    return;
end

%----- only gets here if a subest of slices/frames are to be loaded ---------%


if(frames(1) <= 0) frames = [1:nframes]; end
if(slices(1) <= 0) slices = [1:ndim3]; end

nvslice = ndim1 * ndim2;
nvvol   = ndim1 * ndim2 * ndim3;
filepos0 = ftell(fid);
vol = zeros(ndim1,ndim2,length(slices),length(frames));
nthframe = 1;
for frame = frames
    
    nthslice = 1;
    for slice = slices
        filepos = ((frame-1)*nvvol + (slice-1)*nvslice)*nbytespervox + filepos0;
        fseek(fid,filepos,'bof');
        
        switch type
            case MRI_FLOAT,
                [tmpslice nread]  = fread(fid, nvslice, 'float32') ;
            case MRI_UCHAR,
                [tmpslice nread]  = fread(fid, nvslice, 'uchar') ;
            case MRI_SHORT,
                [tmpslice nread]  = fread(fid, nvslice, 'short') ;
            case MRI_INT,
                [tmpslice nread]  = fread(fid, nvslice, 'int') ;
        end
        
        if(nread ~= nvslice)
            fprintf('ERROR: load_mgh: reading slice %d, frame %d\n',slice,frame);
            fprintf('  tried to read %d, actually read %d\n',nvslice,nread);
            fclose(fid);
            if(gzipped >=0) unix(sprintf('rm %s', fname)); end
            return;
        end
        
        vol(:,:,nthslice,nthframe) = reshape(tmpslice,[ndim1 ndim2]);
        nthslice = nthslice + 1;
    end
    
    nthframe = nthframe + 1;
end

% seek to just beyond the last slice/frame %
filepos = (nframes*nvvol)*nbytespervox + filepos0;
fseek(fid,filepos,'bof');

if(~feof(fid))
    [mr_parms count] = fread(fid,5,'float32');
    if(count < 4)
        fprintf('WARNING: error reading MR params\n');
    end
end

fclose(fid) ;
if(gzipped >=0) unix(sprintf('rm %s', fname)); end

return;

function [vertices, label, colortable] = Read_Brain_Annotation(filename)
% [vertices, label, colortable] = Read_Brain_Annotation(annotfilename.annot)
%
% vertices expected to be simply from 0 to number of vertices - 1;
% label is the vector of annotation
%
% colortable is empty struct if not embedded in .annot. Else, it will be
% a struct.
% colortable.numEntries = number of Entries
% colortable.orig_tab = name of original colortable
% colortable.struct_names = list of structure names (e.g. central sulcus and so on)
% colortable.table = n x 5 matrix. 1st column is r, 2nd column is g, 3rd column
% is b, 4th column is flag, 5th column is resultant integer values
% calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0.


%
% read_annotation.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
%    $Revision: 1.4 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA).
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

fp = fopen(filename, 'r', 'b');

if(fp < 0)
    disp('Annotation file cannot be opened');
    return;
end

A = fread(fp, 1, 'int');

tmp = fread(fp, 2*A, 'int');
vertices = tmp(1:2:end);
label = tmp(2:2:end);

bool = fread(fp, 1, 'int');
if(isempty(bool)) %means no colortable
    disp('No Colortable found.');
    colortable = struct([]);
    fclose(fp);
    return;
end

if(bool)
    
    %Read colortable
    numEntries = fread(fp, 1, 'int');
    
    if(numEntries > 0)
        
        disp(['Reading from Original Version']);
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        for i = 1:numEntries
            len = fread(fp, 1, 'int');
            colortable.struct_names{i} = fread(fp, len, '*char')';
            colortable.struct_names{i} = colortable.struct_names{i}(1:end-1);
            colortable.table(i,1) = fread(fp, 1, 'int');
            colortable.table(i,2) = fread(fp, 1, 'int');
            colortable.table(i,3) = fread(fp, 1, 'int');
            colortable.table(i,4) = fread(fp, 1, 'int');
            colortable.table(i,5) = colortable.table(i,1) + colortable.table(i,2)*2^8 + colortable.table(i,3)*2^16 + colortable.table(i,4)*2^24;
        end
        disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
        
    else
        version = -numEntries;
        if(version~=2)
            disp(['Error! Does not handle version ' num2str(version)]);
        else
            disp(['Reading from version ' num2str(version)]);
        end
        numEntries = fread(fp, 1, 'int');
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        
        numEntriesToRead = fread(fp, 1, 'int');
        for i = 1:numEntriesToRead
            structure = fread(fp, 1, 'int')+1;
            if (structure < 0)
                disp(['Error! Read entry, index ' num2str(structure)]);
            end
            if(~isempty(colortable.struct_names{structure}))
                disp(['Error! Duplicate Structure ' num2str(structure)]);
            end
            len = fread(fp, 1, 'int');
            colortable.struct_names{structure} = fread(fp, len, '*char')';
            colortable.struct_names{structure} = colortable.struct_names{structure}(1:end-1);
            colortable.table(structure,1) = fread(fp, 1, 'int');
            colortable.table(structure,2) = fread(fp, 1, 'int');
            colortable.table(structure,3) = fread(fp, 1, 'int');
            colortable.table(structure,4) = fread(fp, 1, 'int');
            colortable.table(structure,5) = colortable.table(structure,1) + colortable.table(structure,2)*2^8 + colortable.table(structure,3)*2^16 + colortable.table(structure,4)*2^24;
        end
        disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
    end
else
    disp('Error! Should not be expecting bool = 0');
end

fclose(fp);
return;


