function   varargout = Correct_BrainVisa_Median_Mesh(varargin);
%
% Syntax :
%     Surfsulci = Correct_BrainVisa_Median_Mesh(Surfsulci);
%
% This function corrects possible topological errors in brainvisa sulcal median meshes.
%
% Input Parameters:
%        Surfsulci              : Surface variable (file, struct or cellarray).
%
% Output Parameters:
%        Surfout                : Output Surface variable.
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% November 13th 2014
% Version $1.0

%% ============================= Checking Inputs ======================= %%
if nargin == 0
   error('One Input is mandatory');
   return; 
end
if nargin > 1
    error('Too Many Input Parameters');
    return;
end
if nargout > 1
    error('Too Many Output Parameters');
    return;
end
Surfsulci = varargin{1};
%% ========================= End of Checking Inputs ==================== %%

%% ============================= Main Program ========================== %%
Surfsulci = Surface_Checking(Surfsulci);
Nsulc = length(Surfsulci);
for sulc = 1:Nsulc
    tempSurf = Surfsulci(sulc);
    [V,F] = meshfixmex(tempSurf.SurfData.vertices, tempSurf.SurfData.faces);
    tempSurf.SurfData.vertices = V;
    tempSurf.SurfData.faces = F;
    
    S = statistics(tempSurf.SurfData.vertices,tempSurf.SurfData.faces);
    factor = 1;
    
    while S.num_handles ~=0 && factor > 0
        tempSurf = Compute_Hull_from_Surface(tempSurf,factor);
        S = statistics(tempSurf.SurfData.vertices,tempSurf.SurfData.faces);
        factor = factor - 0.1;
    end
    
    [V,F] = meshfixmex(tempSurf.SurfData.vertices, tempSurf.SurfData.faces);
    tempSurf.SurfData.vertices = V;
    tempSurf.SurfData.faces = F;
    Surfsulci(sulc) = tempSurf;
end
%% ====================== End of Main Program ========================== %%

% Outputs
varargout{1} = Surfsulci;