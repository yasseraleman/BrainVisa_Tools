function Surf = Collapse_Surf(Surf);
%
% Syntax :
% Surf = Join(Surf);
%
% This function Collapse all the surface information staraged in Surf into a 
% new Surface
%
% Input Parameters:
%   Surf        : Matlab variable containing all the surfaces.
%
% Output Parameters:
%   OutputFiles  : Atlased Surfaces Files.
%   Surf         : Matlab variable containing the resulting surface.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Plot_oversurf
% Exp_Surf
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1st 2007
% Version $1.0

%============== Creating new structure and joining Surfaces ==============%
cont = 0;
Ns = length(Surf);
vert = [0 0 0];faces = [0 0 0]; normals = [0 0 0];Is = 0;
for i = 1:Ns
    faces = [faces;Surf(i).SurfData.faces+size(vert,1)];
    vert = [vert;Surf(i).SurfData.vertices];
    normals = [normals;Surf(i).SurfData.VertexNormals];
    if ~isfield(Surf(i),'Is')|isempty(Surf(i).Is)
        Surf(i).Is = repmat(i,[size(Surf(i).SurfData.vertices,1) 1]);
    end
    Is = [Is;Surf(i).Is];
end
faces = faces-1;faces(1,:) = [];vert(1,:) = [];normals(1,:) = []; Is(1) = [];
Surft.Name = 'Joined_Surface';
Surft.SurfData.vertices = vert;
Surft.SurfData.faces = faces;
Surft.SurfData.VertexNormals = normals;
Surft.Is = Is;
[Tri] = Vert_Neib(double(Surft.SurfData.faces),size(Surft.SurfData.vertices,1),size(Surft.SurfData.faces,1));
Temp = sum(Tri);
Tri(:,Temp==0) = [];
Surft.Tri = Tri;
Surft.Code = '';
Surft.Type = 'Joint';
%=========================================================================%
Surf = Surft;
return