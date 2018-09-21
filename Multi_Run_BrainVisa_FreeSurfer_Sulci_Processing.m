function Multi_Run_BrainVisa_FreeSurfer_Sulci_Processing(Idfile);

% Idfile = '/media/Data/PROCESSING_RESULTS/PEPS/10-Connect_Stats/PhD_Thesis_Stats-5-6-2014_yaleman/PEPS-Ids_for_unpairedData_5-06-2014.txt';
%Idfile = '/media/Data/PROCESSING_RESULTS/5-freesurfer_processing/Connectome_Ids.txt';


% Idfile = strvcat('HCP_201111-20140807-T1wMPR1');

if exist(Idfile,'file')
    Ids = char(textread(Idfile,'%s'));
else
    Ids = Idfile;
end
Ns = size(Ids,1);
failed = '';
FreeSDir = '/media/Data/PROCESSING_RESULTS/HCP/8-BrainVisaDataBase';
% Ids = char(textread(Idfile,'%s'));
cad = '';
Ido = '';
%%matlabpool('open',4);
for  i = 1:Ns
    
    Id = deblank(Ids(i,:));
    disp(['Processing Subject ' Id '. Subject ' num2str(i) ' of ' num2str(Ns)]);
    delete([ FreeSDir filesep Id '_SULCMETRICS_FAILED.txt']);
    try
        [OutdirMesh] = BrainVisa_FreeSurfer_Sulci_Processing(Id);
    catch
        fid = fopen([ FreeSDir filesep Id '_SULCMETRICS_FAILED.txt'],'wt');
        fclose(fid);
    end
    disp(' ');
end
%%matlabpool('close');
return;
