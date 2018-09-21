function varargout = BrainVisa_stats_to_table(varargin);
%
% Syntax :
%  statTable = BrainVisa_stats_to_table(statFiles, selIds);
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
if nargin == 0
    error('Three Inputs are needed');
    return
end
statFiles = varargin{1};

if nargin == 2
    selIds = varargin{2};
    if exist(selIds,'file')
        selIds = char(textread(selIds,'%s'));
    end
end

if nargin == 3
    statTable = varargin{3};
else
    tempFile = deblank(statFiles(1,:));
    if exist(tempFile,'file')
        [pth,nm,ext] = fileparts(tempFile);
        statTable = [pth filesep 'BrainVISA_General_Table_Orig.txt'];
    else
        error('Incorrect Stat File');
        return;
    end
end
%% ======================= End of Input parameters  ======================%
%% =========================== Main Program ============================= %
Ns = size(statFiles,1);
for i = 1:Ns
    filename = deblank(statFiles(i,:)); % Single Stat FIle
    
    % -------- Creating String Cad to read the stat file
    a = repmat('%s',[1 34]); 
    temps = repmat('Temp',[34 1]);
    b = [temps num2str([1:34]') repmat(',',[34 1])]';
    b = b(:)';
    b(isspace(b)) = [];
    b(end) = [];
    cad = ['[' b '] = textread(''' filename '''' ',' '''' a '''' ',' '''' 'delimiter'  '''' ',' '''' ' ' '''' ');'];
    % -------- End of Creating String Cad to read the stat file
    eval(cad); % Evaluating the string Cad
    
    % -------- Selecting Ids
    IdList = char(Temp1(2:end));
    if nargin < 2
        selIds = IdList;
    end
    Ncol = size(selIds,2);
    ind = find(ismember(IdList(:,1:Ncol),selIds,'rows'));
    IdList = IdList(ind,:);
    ind = [1;ind + 1];
    
    % -------- Sulcus Label
    sulcLabel    = [char(Temp2(2)) '_' char(Temp3(2))];
    
    % -------- Creacting Columns
    maxDepth    = strvcat(['geodDepthMax-'  sulcLabel '(mm)'],char(Temp23(2:end)));
    maxDepth = maxDepth(ind,:);
    minDepth    = strvcat(['geodDepthMin-'  sulcLabel '(mm)'],char(Temp24(2:end)));
    minDepth = minDepth(ind,:);
    meanDepth   = strvcat(['geodDepthMean-' sulcLabel '(mm)'],char(Temp25(2:end)));
    meanDepth = meanDepth(ind,:);
    topLength   = strvcat(['topLengh-'      sulcLabel '(mm)'],char(Temp31(2:end)));
    topLength = topLength(ind,:);
    foldOpening = strvcat(['foldOpening-'   sulcLabel '(mm)'],char(Temp33(2:end)));
    foldOpening = foldOpening(ind,:);
    
    % -------- Creacting the table
    if i == 1
        Nl = size(maxDepth,1);
        repComa = repmat(',',[Nl 1]);
        sIds = strvcat(['SubjectsID'],IdList);
        cad2table = [sIds repComa maxDepth repComa minDepth repComa meanDepth repComa topLength repComa foldOpening];
    else
        cad2table = [cad2table repComa maxDepth repComa minDepth repComa meanDepth repComa topLength repComa foldOpening];
    end
end

% -------- Saving Stat File
fid = fopen(statTable,'wt')
for i = 1:size(cad2table,1)
    fprintf(fid,'%s\n',cad2table(i,:));
end
fclose(fid);

varargout{1} = statTable;
return;