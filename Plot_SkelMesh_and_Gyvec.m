function Plot_SkelMesh_and_Gyvec(SkelFile,GyvecFiles );

SkelFile = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase/subjects/1Test_HCP_899885-20140807-T1wMPR1/t1mri/default_acquisition/RII/mesh/L1Test_HCP_899885-20140807-T1wMPR1_skel.mesh';
GyvecFiles = '/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase/subjects/1Test_HCP_899885-20140807-T1wMPR1/t1mri/default_acquisition/RII/sulcalspams/1Test_HCP_899885-20140807-T1wMPR1_L_F.C.L.a._left_FS.mesh';

Ns = size(GyvecFiles,1);
col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
re = floor(Ns/Ncolor); col = repmat(col,[re+1 1]);

Surfskel = load_mesh(SkelFile);
Plot_Surf(Surfskel);
hold on

for i = 1:Ns
    %GyVec = deblank(GyvecFiles(i,:));
    %Surfgyvec = load_mesh(GyVec);
    X = [Surfgyvec.SurfData.vertices(Surfgyvec.SurfData.faces(:,1),1) Surfgyvec.SurfData.vertices(Surfgyvec.SurfData.faces(:,2),1)];
    Y = [Surfgyvec.SurfData.vertices(Surfgyvec.SurfData.faces(:,1),2) Surfgyvec.SurfData.vertices(Surfgyvec.SurfData.faces(:,2),2)];
    Z = [Surfgyvec.SurfData.vertices(Surfgyvec.SurfData.faces(:,1),3) Surfgyvec.SurfData.vertices(Surfgyvec.SurfData.faces(:,2),3)];
    plot3(X',Y',Z','Color',col(2,:));
end