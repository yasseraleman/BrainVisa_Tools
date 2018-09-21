function Comp_MeanGyriSpam
FilenameSpam = '/media/COSAS/Test/Joost/spanmesh/ASPER_00001__101-20060510_L_gyrus8_left_gyvec.mesh';
Surf = load_mesh(FilenameSpam);
GyriParc = '/media/COSAS/Test/Joost/spanmesh/ASPER_00001__101-20060510_Lwhite_gyri_default_session_auto.tex';
[Textth,Format] = read_texBrainvisa(GyriParc);
FilenameMesh = '/media/COSAS/Test/Joost/spanmesh/ASPER_00001__101-20060510_Lwhite.mesh';
SurfM = load_mesh(FilenameMesh);

a = (Surf.SurfData.vertices(Surf.SurfData.faces(:,1),1)-Surf.SurfData.vertices(Surf.SurfData.faces(:,2),1)).^2;
b = (Surf.SurfData.vertices(Surf.SurfData.faces(:,1),2)-Surf.SurfData.vertices(Surf.SurfData.faces(:,2),2)).^2;
c = (Surf.SurfData.vertices(Surf.SurfData.faces(:,1),3)-Surf.SurfData.vertices(Surf.SurfData.faces(:,2),3)).^2;
d = sqrt(a+b+c);

X = [Surf.SurfData.vertices(Surf.SurfData.faces(:,1),1) Surf.SurfData.vertices(Surf.SurfData.faces(:,2),1)];
Y = [Surf.SurfData.vertices(Surf.SurfData.faces(:,1),2) Surf.SurfData.vertices(Surf.SurfData.faces(:,2),2)];
Z = [Surf.SurfData.vertices(Surf.SurfData.faces(:,1),3) Surf.SurfData.vertices(Surf.SurfData.faces(:,2),3)];
P1 = [X(:,1) Y(:,1) Z(:,1)];
P2 = [X(:,2) Y(:,2) Z(:,2)];
Mat = [1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1];
P1 = Mat*[P1 ones(size(P1,1),1)]';
P1 = P1';P1(:,4) = [];
P2 = Mat*[P2 ones(size(P2,1),1)]';
P2 = P2';P2(:,4) = [];
X = [P1(:,1) P2(:,1)];
Y = [P1(:,2) P2(:,2)];
Z = [P1(:,3) P2(:,3)];
SurfM.SurfData.vertices = Mat*[SurfM.SurfData.vertices ones(size(SurfM.SurfData.vertices,1),1)]';
SurfM.SurfData.vertices = SurfM.SurfData.vertices';SurfM.SurfData.vertices(:,4) = [];
SurfM.Is = Textth.Values(:,3);
Plot_Surf(SurfM,[0.8],'y');hold on;
line(X',Y',Z','Color',[0 1 0])/media/COSAS/scripts/BrainVisaTools