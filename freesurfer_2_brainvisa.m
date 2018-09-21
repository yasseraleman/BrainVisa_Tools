opts.freesurferdir = '/media/yaleman/MyThings4T/HCPData/5-freesurfer_processing';
opts.subjid = '1SUBJECT_100206_TestSubject';
outDir = '';

cad = ['mri_convert -i ' opts.freesurferdir  filesep opts.subjid filesep 'mri' filesep 'ribbon.mgz -o ' opts.freesurferdir  filesep opts.subjid filesep 'tmp' filesep 'ribbon.nii'] ;
system(cad);
Vfs = spm_vol([opts.freesurferdir  filesep opts.subjid filesep 'tmp' filesep 'ribbon.nii']);
opts.freevol = Vfs;
delete(Vfs.fname);

surfFiles = {[opts.freesurferdir filesep opts.subjid filesep 'surf' filesep 'lh.pial'],...
    [opts.freesurferdir filesep opts.subjid filesep 'surf' filesep 'rh.pial'],...
    [opts.freesurferdir filesep opts.subjid filesep 'surf' filesep 'lh.white'],...
    [opts.freesurferdir filesep opts.subjid filesep 'surf' filesep 'rh.white']};
for i = 1:Nsurf
    Surf = Read_Surface( surfFiles{i});
    
    % Reading Talairach Transformation
    Taltransf = [opts.freesurferdir filesep opts.subjid filesep 'mri' filesep 'transforms' filesep 'talairach.lta'];
    cras1 = textread(Taltransf,'%s',5,'headerlines',20);cras = char(cras1);cras = [str2num(cras(3,:))  str2num(cras(4,:)) str2num(cras(5,:))];
    
    % Adding Thalairach center
    Surf.SurfData.vertices =Surf.SurfData.vertices+repmat(cras,[size(Surf.SurfData.vertices,1) 1]); % adding RAS center
    
    % Converting to FreeSurfer Coordinate system
    Surf = freesCS2brainvisaCS(Surf,opts.freevol,'f2b');
    
    [pthDir, nam, ext] = fileparts(surfFiles{i});
    if isempty(outDir)
        outDir = pthDir;
    end
    outFile = [outDir filesep nam ext '.mesh'];
    save_mesh(Surf,outFile);
end