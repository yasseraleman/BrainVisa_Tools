function [OutImage,Out_mesh] = Gyral_WM_Skeleton_Computation(InputImage, OutputDirectory);


Files2delete = '';
cad = ['VipSingleThreshold -t 150 -m gt -i ' InputImage ' -o ' OutputDirectory filesep 'tmp_thr.nii'];
system(cad);
Files2delete = strvcat(Files2delete,[OutputDirectory filesep 'tmp_thr.nii']);
Files2delete = strvcat(Files2delete,[OutputDirectory filesep 'tmp_thr.nii.minf']);

cad = ['AimsMedianSmoothing --dz 3 -i ' OutputDirectory filesep 'tmp_thr.nii -o ' OutputDirectory filesep 'tmp_smooth.nii'];
system(cad);
Files2delete = strvcat(Files2delete,[OutputDirectory filesep 'tmp_smooth.nii']);
Files2delete = strvcat(Files2delete,[OutputDirectory filesep 'tmp_smooth.nii.minf']);

cad = ['VipSkeleton -sk s -im a -gcs 2 -e 1.5 -i ' OutputDirectory filesep 'tmp_smooth.nii -so ' OutputDirectory filesep 'tmp_skel.nii -g ' OutputDirectory filesep 'tmp_smooth.nii'];
system(cad);
Files2delete = strvcat(Files2delete,[OutputDirectory filesep 'tmp_skel.nii']);
Files2delete = strvcat(Files2delete,[OutputDirectory filesep 'tmp_skel.nii.minf']);

cad = ['VipSingleThreshold -m gt -t 1 -c b -i ' OutputDirectory filesep 'tmp_skel.nii -o ' OutputDirectory filesep 'tmp_bin.nii'];
system(cad);
Files2delete = strvcat(Files2delete,[OutputDirectory filesep 'tmp_bin.nii']);
Files2delete = strvcat(Files2delete,[OutputDirectory filesep 'tmp_bin.nii.minf']);

cad = ['AimsMorphoMath -r 3.5 -m clo -i ' OutputDirectory filesep 'tmp_bin.nii  -o ' OutputDirectory filesep 'final_skel.nii'];
system(cad);
OutImage = [OutputDirectory filesep 'final_skel.nii'];
Files2delete = strvcat(Files2delete,[OutputDirectory filesep 'final_skel.nii.minf']);

cad = ['AimsMesh --smooth --smoothIt 20 -i ' OutputDirectory filesep 'final_skel.nii  -o ' OutputDirectory filesep 'final_skel.mesh'];
system(cad);

a = dir([OutputDirectory filesep 'final_skel_255_*.mesh']);
for i = 1:length(a);
    Files2delete = strvcat(Files2delete,[OutputDirectory filesep a(i).name],[OutputDirectory filesep a(i).name '.minf']);
    file = [OutputDirectory filesep a(i).name];
    Surf(i) = Read_Surface(file);
end
Surfj = Compound_Surf(Surf);
Out_mesh = save_mesh(Surfj,[OutputDirectory filesep 'final_skel.mesh']);


for i = 1:size(Files2delete,1)
    cad = ['rm -r ' deblank(Files2delete(i,:))];
    system(cad);
end
return;

