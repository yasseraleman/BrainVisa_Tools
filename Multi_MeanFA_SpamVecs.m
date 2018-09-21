function Multi_MeanFA_SpamVecs;
Idfile = '/media/Data/Joost/gspan/gyvecs/gyri/Ids_ASP.txt';
InputFolder = '/media/Data/Joost/gspan/gyvecs/gyri';
FreeSDir = '/media/Data/PROCESSING_RESULTS/ASPERGER/5-freesurfer_processing';
ConnectDir = '/media/Data/PROCESSING_RESULTS/ASPERGER/7-connectome';
OutputDir = '/media/Data/PROCESSING_RESULTS/ASPERGER';
Ids = char(textread(Idfile,'%s'));
cadtot = '';
% Ids = Ids(82,:);
for i = 1:size(Ids,1)
    tic
    disp(i);
    Id = deblank(Ids(i,:));
    Atlasfile = ['/media/Data/Joost/gspan/volumes/gyri/' Id '.wmparc.gyri.nii.gz'];
    try
        cad = ['gunzip -df ' Atlasfile];
        system(cad);
    end
     Atlasfile = ['/media/Data/Joost/gspan/volumes/gyri/' Id '.wmparc.gyri.nii'];
    %try
        [OutFile, cads] = Atlasing_Skeletonized_FA(OutputDir,Id, Atlasfile);
        %[cads] = MeanFA_SpamVecs(InputFolder, Id,FreeSDir,ConnectDir);
        if i  == 1
            cadnames = cads(1,:);
        end
        cadtot = strvcat(cadtot,cads(2,:));
   % end
    toc
end
cadtot = strvcat(cadnames,cadtot);
fid = fopen('/home/yaleman/MeanMD_WM_Skel_ASP.txt','wt');
for i = 1:size(cadtot,1)
    line = cadtot(i,:);
    %     ind = strfind(line,'.');
    %     if ~isempty(ind)
    %         line(ind) = ',';
    %     end
    fprintf(fid, '%s\n', line);
end
fclose all;
disp(['Output file: ==>>  /home/yaleman/MeanFA_GSPams_PEP.txt']);