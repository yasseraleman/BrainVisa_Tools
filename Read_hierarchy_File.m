function varargout = Read_hierarchy_File(varargin);
%
% Syntax :
%  [allStructLabels,allColors] = Read_hierarchy_File(hieFilename);
%
% This script reads brainvisa hierarchy file. 
% (i.e ~/brainvisa/share/brainvisa-share-4.5/nomenclature/hierarchy/sulcal_root_colors.hie).
%
% Input Parameters:
%       hieFilename             : Brainvisa hierarchy filename.
%
% Output Parameters:
%      allStructLabels          : Structures Labels.
%      allColors                : Structures Colors.
%
%
% See also: Sulci_Nodal_Processing  Compute_Node_Metrics
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2016
% Version $1.0

%% ===================== Checking Input parameters =======================%
if nargin == 0
    
    [~,B] = system('locate sulcal_root_colors.hie');
    ind = find(isspace(B));
    hieFilename = deblank(B(1:ind(1)-1));
elseif nargin == 1
    hieFilename = varargin{1};
else
    error('To many inputs');
    return;
end
%% ================== End of Checking Input parameters ===================%

%% ==================== Reading hierearchy Filename ======================%
fid = fopen(hieFilename);
contcolor = 0;
cont = 0;
alllines = '';
while 1
    
    line = fgetl(fid);
    if ~ischar(line),   break,   end
    
    if ~isempty(line)
        cont = cont + 1;
        alllines = strvcat(alllines,line);
        
        indc = strfind(line,'color');
        %         indl = strfind(line,'label');
        indn = strfind(line,'name');
        
        if isempty(indc)
            indc = 0;
        end
        %         if isempty(indl)
        %             indl = 0;
        %         end
        if isempty(indn)
            indn = 0;
        else
            if indn(1)~=1
                indn = 0;
            end
        end
        %         if indc(1)==1||indl(1)==1||indn(1)==1
        if indc(1)==1||indn(1)==1
            
            contcolor = contcolor + 1;
            clines(contcolor) = cont;
        end
    end
end
fclose(fid);
clear cont;
%% ================== End of Reading hierearchy Filename ==================%

%% ============ Detecting colors and its corresponding labels =============%

a1 = (clines(1:end-1)-clines(2:end));
% a2 = (clines(2:end-1)-clines(3:end));
indexes = clines(find(a1==-1|a1==-2));

ind2del = find(sum(ismember([cellstr(alllines(indexes,1:5)) cellstr(alllines(indexes+1,1:5)) cellstr(alllines(indexes+2,1:5))],{'name'}),2)>=2);
indexes(ind2del) = [];

ind2del = find(ismember([cellstr(alllines(indexes+1,1:5))],{'*BEGI'}));
indexes(ind2del) = [];

allStructLabels = {};
for i  = 1:length(indexes)
    
    switch alllines(indexes(i),1:5)
        case 'color'
            [~,r,g,b] = strread(alllines(indexes(i),:),'%s%u%u%u','delimiter',' ');
            allColors(i,:) = [r g b];
         case 'name '
            [~,label] = strread(alllines(indexes(i),:),'%s%s','delimiter',' ');
            allStructLabels(i,1) = label;
    end
    
    switch alllines(indexes(i)+1,1:5)
        case 'color'
            [~,r,g,b] = strread(alllines(indexes(i)+1,:),'%s%u%u%u','delimiter',' ');
            allColors(i,:) = [r g b];
         case 'name '
             [~,label] = strread(alllines(indexes(i)+1,:),'%s%s','delimiter',' ');
             allStructLabels(i,1) = label;
    end
    
    switch alllines(indexes(i)+2,1:5)
        case 'color'
            [~,r,g,b] = strread(alllines(indexes(i)+2,:),'%s%u%u%u','delimiter',' ');
            allColors(i,:) = [r g b];
        case 'name '
            [~,label] = strread(alllines(indexes(i)+2,:),'%s%s','delimiter',' ');
            allStructLabels(i,1) = label;
    end
end
%% ========= End of Detecting colors and its corresponding labels =========%

% Ouputs
varargout{1} = allStructLabels;
varargout{2} = allColors;
return;




