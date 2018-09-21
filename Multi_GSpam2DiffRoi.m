function Multi_GSpam2DiffRoi;

% ArgFiles =  strvcat('/media/COSAS/Test/Joost/aquivan/aquivienen/ASPER_00001__101-20060510_FS2BV.lh.white+tal.aims_Tmtkmtri.arg','/media/COSAS/Test/Joost/aquivan/aquivienen/ASPER_00001__101-20060510_FS2BV.rh.white+tal.aims_Tmtkmtri.arg'); %Mandatory 
% Imfile = '/media/COSAS/Test/Joost/aquivan/aquivienen/ASPER_00001__101-20060510_T1.nii'; %Mandatory 
% Famap = '/media/Data/PROCESSING_RESULTS/ASPERGER/7-connectome/ASPER_00001__101-20060510/dtreconst/fsl/ASPER_00001__101-20060510_fa.nii'; %Mandatory
% Transf = '/media/Data/PROCESSING_RESULTS/ASPERGER/7-connectome/ASPER_00001__101-20060510/transforms/ASPER_00001__101-20060510_T2Bet_2_diffAffine.txt';
InputFolder = '/media/LA-PUBLIC/TEA_ABIDE_brainvisa43/FS_gyral_span/lobar';
% ResultsFolder = '/media/Data/PROCESSING_RESULTS/ASPERGER/7-connectome';
% T1Files = sel_files(InputFolder,'*T1.nii');
% ind = strfind(T1Files(1,:),filesep);
% Ids = T1Files(:,ind(end)+1:ind(end)+25);
% Ids = unique(Ids,'rows');
% Ns = size(Ids,1);
% cadtotvol = 'Subjects          ';
% cadtotfa = 'Subjects          ';
% for i = 1:Ns
%     disp(['Processing Subject ' deblank(Ids(i,:)) ' : Subject ' num2str(i) ' of ' num2str(Ns) ]);
%     Lspamvecs = sel_files(InputFolder,[deblank(Ids(i,:)) '*L*.mesh']);
%     Rspamvecs = sel_files(InputFolder,[deblank(Ids(i,:)) '*R*.mesh']);
%     Larg = sel_files(InputFolder,[deblank(Ids(i,:)) '*.lh.*.arg']);
%     Rarg = sel_files(InputFolder,[deblank(Ids(i,:)) '*.rh.*.arg']);
%     ArgFiles = strvcat(Larg,Rarg);
%     Famap = [ResultsFolder filesep deblank(Ids(i,:)) filesep 'dtreconst' filesep 'fsl' filesep deblank(Ids(i,:)) '_fa.nii'];
%     Transf = [ResultsFolder filesep deblank(Ids(i,:)) filesep 'transforms' filesep deblank(Ids(i,:)) '_T2Bet_2_diffAffine.txt'];
%     OutAtlasFile = [ResultsFolder filesep deblank(Ids(i,:)) filesep 'preproc' filesep deblank(Ids(i,:)) '_GyralSpam_Atlas.nii']
%     [OutAtlasFile,Oufiles] = GSpam2DiffRoi(Lspamvecs,Rspamvecs,ArgFiles,deblank(T1Files(i,:)),Famap,Transf,'Fa2Gy',OutAtlasFile);
%     Vseg = spm_vol(deblank(Oufiles(2,:)));
%     Vfa = spm_vol(deblank(Oufiles(1,:)));
%     Is = spm_read_vols(Vseg);
%     If = spm_read_vols(Vfa);
%     voxsize = sqrt(sum(Vseg.mat(1:3,1:3).^2)); voxvol = prod(voxsize);
%     %ind = unique(Is(:));ind(ind==0) = [];
%     ind = [1:3 5:34 101:103 105:134];
%     if i == 1
%         for j = 1:length(ind)
%             if ind(j)<100
%                 cadtotfa = [cadtotfa '     ' ';Left_Gyrus_' num2str(ind(j))];
%                 cadtotvol = [cadtotvol '     ' ';Left_Gyrus_' num2str(ind(j))];
%             else
%                 cadtotfa = [cadtotfa '     ' ';Right_Gyrus_' num2str(ind(j)-100)];
%                 cadtotvol = [cadtotvol '     ' ';Right_Gyrus_' num2str(ind(j)-100)];
%             end
%         end
%     end
%     cadfa = deblank(Ids(i,:));
%     cadvol = deblank(Ids(i,:));
%     for j=1:length(ind)
%         indv = find(Is ==ind(j));
%         if ~isempty(indv)
%             meanfa = mean(If(indv));
%             vol = length(indv)*voxvol/1000;
%         else
%             meanfa = 0;
%             vol = 0;
%         end
%         cadfa = [cadfa '     ;' num2str(meanfa)];
%         cadvol = [cadvol '     ;' num2str(vol)];
%     end
%     cadtotfa = strvcat(cadtotfa,cadfa);
%     cadtotvol = strvcat(cadtotvol,cadvol);
% end

%% NUEVOS Reclutamientos
%ResultsFolder = '/media/Data/PROCESSING_RESULTS/NUEVRECLU/7-connectome';
T1Files = sel_files(InputFolder,'*T1.nii');
ind = strfind(T1Files(1,:),filesep);
%Ids = T1Files(:,ind(end)+1:ind(end)+25);
%Ids = unique(Ids,'rows');
Ns = size(T1Files,1);
cad = ['IDs                     LH_Gyrus_1        LH_Gyrus_3        LH_Gyrus_5        LH_Gyrus_6        RH_Gyrus_1        RH_Gyrus_3        RH_Gyrus_5        RH_Gyrus_6'];
disp(['No      ' cad]);
FailedIds = '';
for i = 1:Ns
    [pth,nm,ext] = fileparts(deblank(T1Files(i,:)));
    indn = strfind(nm,'_');
    Id = nm(1:indn(end)-1);
   % disp(['Processing Subject ' Id ' : Subject ' num2str(i) ' of ' num2str(Ns) ]);
    Lspamvecs = sel_files(InputFolder,[Id '*L*_gyvec.mesh']);
    Rspamvecs = sel_files(InputFolder,[Id '*R*_gyvec.mesh']);
%     for j = 1:
%     end
    Larg = sel_files(InputFolder,[Id '*_lh.white+tal.aims_Tmtkmtri.arg']);
    Rarg = sel_files(InputFolder,[Id '*_rh.white+tal.aims_Tmtkmtri.arg']);
    ArgFiles = strvcat(Larg,Rarg);
    try
        %Famap = [ResultsFolder filesep Id filesep 'dtreconst' filesep 'fsl' filesep Id '_fa.nii'];
        %Transf = [ResultsFolder filesep Id filesep 'transforms' filesep Id '_T2Bet_2_diffAffine.txt'];
        OutAtlasFile = [InputFolder filesep Id '_GyralSpam_Atlas.nii'];
        %[OutAtlasFile,Oufiles] = GSpam2DiffRoi(Lspamvecs,Rspamvecs,ArgFiles,deblank(T1Files(i,:)),Famap,Transf,'Fa2Gy',OutAtlasFile);
        [OutAtlasFile, Vols] = GSpam2Roi(Lspamvecs,Rspamvecs,ArgFiles,deblank(T1Files(i,:)),OutAtlasFile);
         cadt = [Id '        ' num2str(Vols(1)) '             ' num2str(Vols(2)) '              ' num2str(Vols(3)) '             ' num2str(Vols(4))...
                    '             ' num2str(Vols(5)) '             ' num2str(Vols(6)) '              ' num2str(Vols(7)) '             ' num2str(Vols(8))];
                disp([sprintf('%.3d',i) '     ' cadt]);
          cad = strvcat(cad,cadt);
%         Vseg = spm_vol(deblank(Oufiles(2,:)));
%         Vfa = spm_vol(deblank(Oufiles(1,:)));
%         Is = spm_read_vols(Vseg);
%         If = spm_read_vols(Vfa);
%         voxsize = sqrt(sum(Vseg.mat(1:3,1:3).^2)); voxvol = prod(voxsize);
%         %ind = unique(Is(:));ind(ind==0) = [];
%         ind = [1:3 5:34 101:103 105:134];
%         if i == 1
%             for j = 1:length(ind)
%                 if ind(j)<100
%                     cadtotfa = [cadtotfa '     ' ';Left_Gyrus_' num2str(ind(j))];
%                     cadtotvol = [cadtotvol '     ' ';Left_Gyrus_' num2str(ind(j))];
%                 else
%                     cadtotfa = [cadtotfa '     ' ';Right_Gyrus_' num2str(ind(j)-100)];
%                     cadtotvol = [cadtotvol '     ' ';Right_Gyrus_' num2str(ind(j)-100)];
%                 end
%             end
%         end
%         cadfa = deblank(Ids(i,:));
%         cadvol = deblank(Ids(i,:));
%         for j=1:length(ind)
%             indv = find(Is ==ind(j));
%             if ~isempty(indv)
%                 meanfa = mean(If(indv));
%                 vol = length(indv)*voxvol/1000;
%             else
%                 meanfa = 0;
%                 vol = 0;
%             end
%             cadfa = [cadfa '     ' num2str(meanfa)];
%             cadvol = [cadvol '     ' num2str(vol)];
%         end
%         cadtotfa = strvcat(cadtotfa,cadfa);
%         cadtotvol = strvcat(cadtotvol,cadvol);
    catch
        FailedIds = strvcat(FailedIds,Id);
    end
end

%% Saving Results

 fid = fopen([InputFolder filesep 'Gyral_Volume_BasedOnSpam.txt'],'wt');
 for i = 1:size(cad,1)
     fprintf(fid, '%s\n', cad(i,:));
 end
 fclose all;
 disp(['Output file: ==>>  ' InputFolder filesep 'Gyral_Volume_BasedOnSpam.txt']);
% 
% fid = fopen(['/media/Data/PROCESSING_RESULTS/ASPERGER/Volume_Gyral_SPam.txt'],'wt');
% for i = 1:size(cadtotvol,1)
%     fprintf(fid, '%s\n', cadtotvol(i,:));
% end
% fclose all;
% disp(['Output file: ==>>  ' InputFolder filesep 'Volume_Gyral_SPam.txt']);


