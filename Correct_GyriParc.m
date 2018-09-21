function Correct_GyriParc
gtexfiles = sel_files('/media/COSAS/Test/Joost/parc/temp','*_gyri_default_session_auto.tex');
stexfiles = sel_files('/media/COSAS/Test/Joost/parc/temp','*_sulci_default_session_auto.tex');
% Sfiles = sel_files('/media/COSAS/Test/Joost/parc/temp','*.mesh');
Ns = size(gtexfiles,1);
for j =1:Ns
    %disp(['Computing Subject ']);
    [pth,nam,ext] = fileparts(deblank(gtexfiles(j,:)));
    ind = strfind(nam,'-');
    ID = nam(1:ind+15);
    disp(['Computing Subject:  ==>  ' ID]);
    Sfile = [pth filesep ID '.mesh'];
    [Text,Format] = read_texBrainvisa(deblank(gtexfiles(j,:)));
    [OutFiles, SurfF] = Exp_Surf(Sfile, '0', '','', 'imp','n');
    Surf = SurfF{1};
    Npoints = size(Text.Values,1);
    Surf.Is = Text.Values(:,2);
    Surft = Surf;
    ind = find(Surf.Is == 0);
    str = unique(Surf.Is);Surfs = Surf;
    Surfs.Is = Text.Values(:,1);
    for i =1:size(str,1);
        t = zeros(Npoints,1);
        ind = find(Surf.Is ==str(i));
        t(ind) = 1;
        indz = find(t==0);
        t(indz) = 2;
        Surft.Is = t;
        %Plot_Surf(Surft);
        [Surft] = Surf_Corr(Surft);
        indr = find(Surft.Is ==2);
        a = ismember(ind,indr);
        Surfs.Is(ind(a)) = 0;
        %[Surf] = Surf_Corr(Surf); Para arreglar el Gyri
    end
    ind = find(Text.Values(:,1)-Surfs.Is);
    disp([num2str(length(unique(Text.Values(:,1)-Surfs.Is))-1) ' sulci were corrected:' num2str(length(ind)) ' Points']);
    [pth,nam,ext] = fileparts(deblank(stexfiles(j,:)));
    [Text,Format] = read_texBrainvisa(deblank(stexfiles(j,:)));
    Text.Values = Surfs.Is;
    NsTxtfile = [pth filesep nam '_mod.tex'];
    [Txtfile] = save_texBrainvisa(Text,NsTxtfile,Format);
    fclose all;
end