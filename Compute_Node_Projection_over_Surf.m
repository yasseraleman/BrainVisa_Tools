function varargout = Compute_Node_Projection_over_Surf(varargin)
%
% Syntax :
%       [Is] = Compute_Sulcal_Spam(refSurf, nodeSurf, opts);
%
% This function computes sulcal spam intercepting the median mesh normals
% in the Pial surface
%
%
% Input Parameters:
%       refSurf             : Pial Surface (Matlab Format)
%       nodeSurf                   : Median Mesh (Matlab Format).
%
% Output Parameters:
%      Surfwidth                : Sulcal Spam lines. They can be plotted
%                                 using Plot_Surf.
%      sulcalWidth              : Sulcal Width Vector.
%
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% February 13th 2015
% Version $1.0

%% =========================== Input parameters  =========================%
if nargin<2 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return
else
    refSurf  = varargin{1};
    nodeSurf = varargin{2};
    refSurf = Surface_Checking(refSurf);
    nodeSurf = Surface_Checking(nodeSurf);
end
%% =========================== Input parameters  =========================%
if nargin < 3
    opts.distthresh  = 5; % mm. Maximum Sulcal Width
    opts.epsilon = 10^-5;          % Tolerance to estimate barycentric coordinates
    opts.verbose = 1;              % Maximum Sulcal Width
else
    opts = varargin{3};
    if ~isfield(opts,'maxsulcalwidth')
        opts.maxsulcalwidth = 12;      % Maximum Sulcal Width
    end
    if ~isfield(opts,'epsilon')
        opts.epsilon = 1;              % Tolerance to estimate barycentric coordinates
    end
    if ~isfield(opts,'verbose')
        opts.verbose = 1;              % Maximum Sulcal Width
    end
end
if nargin > 3
    error('To many inputs');
    return;
end
if nargout > 1
    error('To many outputs');
    return;
end
%% ==================== End of Input parameters  =========================%

% % % % load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/matlab_test1.mat');

%% ========================== Main Program ============================== %
if opts.verbose
    disp(' ');
    disp('Computing Projection lines ... ');
    tic;
end

if ~isfield(nodeSurf.SurfData,'VertexNormals');
    Normals = Compute_Surface_Normals(nodeSurf);
    
else
    Normals = nodeSurf.SurfData.VertexNormals;
end

% Detecting Sulcal border
Surft = Sulcal_Face_Labelling(nodeSurf);
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





% Surface vertices and faces
VertP = refSurf.SurfData.vertices;FacesP = refSurf.SurfData.faces;
Nvrefsurf = size(VertP,1);

%  Surface planes
np = cross(VertP(FacesP(:,1),:)-VertP(FacesP(:,2),:),VertP(FacesP(:,3),:)-VertP(FacesP(:,2),:),2); Dp = -1*dot(np,VertP(FacesP(:,1),:),2);

% Maximum edge length
maxdist = sqrt(max(sum((VertP(FacesP(:,1),:)-VertP(FacesP(:,2),:)).^2,2)));
maxdist = max(maxdist,sqrt(max(sum((VertP(FacesP(:,1),:)-VertP(FacesP(:,3),:)).^2,2))));
maxdist = max(maxdist,sqrt(max(sum((VertP(FacesP(:,2),:)-VertP(FacesP(:,3),:)).^2,2))));

% Creating perpendicular lines to sulcus walls
P1 = nodeSurf.SurfData.vertices; %P2 = nodeSurf.SurfData.vertices+Normals; 

P1 =[ P1;[nodeSurf.SurfData.vertices(Coord(:,1),:)+nodeSurf.SurfData.vertices(Coord(:,2),:)]/2];
Normals = [Normals;[Normals(Coord(:,1),:)+Normals(Coord(:,2),:)]/2];
P2 = P1+ Normals;

Npoints = length(P1);
tempValue = ones(Npoints,1);
tempValue(size(nodeSurf.SurfData.vertices,1)+1:Npoints) = 2;

Is = zeros(Nvrefsurf,1);
for po =1:Npoints
    
    % Distance from all Surface points to the line defined by the normal
    % vector
    A = (P1(po,1) - VertP(:,1)).^2 + (P1(po,2) - VertP(:,2)).^2 +(P1(po,3) - VertP(:,3)).^2;
    B = (P2(po,1) - P1(po,1)).*(P1(po,1) - VertP(:,1)) + (P2(po,2) - P1(po,2)).*(P1(po,2) - VertP(:,2)) +( P2(po,3) - P1(po,3)).*(P1(po,3) - VertP(:,3));
    C = (P2(po,1) - P1(po,1)).^2 + (P2(po,2) - P1(po,2)).^2 +(P2(po,3) - P1(po,3)).^2;
    t = -1*dot((repmat(P1(po,:),[Nvrefsurf 1]) - VertP),repmat(Normals(po,:),[Nvrefsurf 1]) ,2);
    d = sqrt(A+ 2*t.*B+ t.^2*C);
    
    % Points with distance to the line higher than the maximum edge length
    % will be removed
    vertx = find(abs(d)<maxdist);
    indf = ismember(FacesP,vertx);
    [rindf,cindf] = find(indf); %#ok
    rindf = unique(rindf);
    NintFaces = length(rindf); % Number of remaining faces

    % Interception between face planes and the line defined by the normal
    % vector
    Num = dot(np(rindf,:),repmat(P1(po,:),[length(rindf) 1]),2)+Dp(rindf); Den = dot(np(rindf,:),repmat(P2(po,:)-P1(po,:),[length(rindf) 1]),2); % Creating Lines
    t = -1*Num./Den;clear Num Den; % Lines parameters
    xint = single(repmat(P1(po,1)',[NintFaces 1])+t.*(repmat((P2(po,1)-P1(po,1))',[NintFaces 1]))); % Line parametrization
    yint = single(repmat(P1(po,2)',[NintFaces 1])+t.*(repmat((P2(po,2)-P1(po,2))',[NintFaces 1])));
    zint = single(repmat(P1(po,3)',[NintFaces 1])+t.*(repmat((P2(po,3)-P1(po,3))',[NintFaces 1])));clear t;
    
    % Computing barycentric coordinates
    v2 = [xint yint zint]- VertP(FacesP(rindf,1),:);
    v0 = VertP(FacesP(rindf,3),:) - VertP(FacesP(rindf,1),:); % dot products
    v1 = VertP(FacesP(rindf,2),:) - VertP(FacesP(rindf,1),:);
    dot00 = dot(v0, v0, 2);
    dot01 = dot(v0, v1, 2);
    dot11 = dot(v1, v1, 2);
    dot02 = dot(v0, v2, 2);
    dot12 = dot(v1, v2, 2);
    invDenom = 1 ./ (dot00 .* dot11 - dot01 .* dot01);
    u = (dot11 .* dot02 - dot01 .* dot12) .* invDenom;
    v = (dot00 .* dot12 - dot01 .* dot02) .* invDenom;
    
    % Detecting faces with an interior interception
    out = u >= -opts.epsilon & v >= -opts.epsilon & (u+v-opts.epsilon) <= 1;
    ind = find(out);
    
    [inte,indord] = unique([xint(ind) yint(ind) zint(ind)],'rows');ind = ind(indord); % Line interception with pial surface
    
    if ~isempty(inte)
        % Apply distance threshold. Distance between the interception point
        % and the starting point will be lower than opts.distthresh.
        distances = sqrt(sum((inte - repmat(P1(po,:),[size(inte,1) 1])).^2,2)); % Computing distance to the pial surface
        ind2del = distances > opts.distthresh;
        distances(ind2del) = [];
        inte(ind2del,:) = [];
        ind(ind2del) = [];
        
        if ~isempty(inte)
            % Keep only the points in the same direction than the normal vector.
            orientsign = sign((dot(inte -repmat(P1(po,:),[size(inte,1) 1]),repmat(P2(po,:)-P1(po,:),[size(inte,1) 1]),2))); % Detecting Orientation
            ind_same_dir = (orientsign>0);
            if ~isempty(ind_same_dir)
                % Select the closest point and assing the label 1 to the vertices in the interception face. 
                selectFaces = ind(ind_same_dir);
                [~,loc] = min(distances(ind_same_dir));
                interFace = selectFaces(loc,:);
                Is(FacesP(rindf(interFace),:)) = tempValue(po);
            end
        end
    end 
end
%% ===================== End of Main Program ============================ %

% Outputs
varargout{1} = Is;
if opts.verbose
    toc;
end
return;