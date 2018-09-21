function Mean_CortThickness;
%pth = '/media/COSAS/Test/Joost'
pth = '/media/LA-PUBLIC/curvature/';
gtexfiles = sel_files(pth,'*_Lwhite_gyri_default_session_auto.tex');
ctexfiles = sel_files(pth,'*_Lwhite_curv.tex');
%ctexfiles = sel_files('/media/COSAS/Test/Joost/','*_curv_bary.tex');
thtexfiles = sel_files(pth,'*_L_wm_thick.tex');
% thtexfiles = ctexfiles;
Ns = size(gtexfiles,1);Ids = '';
for j =1:Ns
    %disp(['Computing Subject ']);
    [pth,nam,ext] = fileparts(deblank(ctexfiles(j,:)));
    ind = strfind(nam,'-');
    ID = nam(1:ind+15);
    Ids = strvcat(Ids,ID);
    disp(['Computing Subject:  ==>  ' ID]);
    Sfile = ['/media/COSAS/Test/Joost/parc/temp' filesep ID '.mesh'];
    [Textg,Format] = read_texBrainvisa(deblank(gtexfiles(j,:)));
    [OutFiles, SurfF] = Exp_Surf(Sfile, '0', '','', 'imp','n');
    Surf = SurfF{1};
    Npoints = size(Textg.Values,1);
    Surf.Is = Textg.Values(:,3);
    ind = find(Surf.Is ==0);
    str = unique(Surf.Is);
    Ns = max(str);Surf.Is(ind) = Ns+1;
    [Surf] = Surf_Corr(Surf);
    ind = find(Surf.Is ==Ns+1); Surf.Is(ind) = 0;
    [Textc,Format] = read_texBrainvisa(deblank(ctexfiles(j,:)));
    ThSup = 0.1*pi;
    ThInf = -0.1*pi;
    indvc = find((Textc.Values(:,1)<=ThSup)&(Textc.Values(:,1)>=ThInf));
    
    [Textth,Format] = read_texBrainvisa(deblank(thtexfiles(j,:)));
    
    for i =1:size(str,1);
        if (str(i)~=0)|(str(i)~=0)|(str(i)~=0)
            ind = find(Surf.Is==str(i));
            indi = ismember(ind,indvc);
            if sum(indi)==0
                disp(['Subject:  ==>  ' ID ' BLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA']);
            end
            cth(j,i) = mean(Textth.Values(ind(indi),1));
            stdcth(j,i) = std(Textth.Values(ind(indi),1));
            fv.vertices = Surf.SurfData.vertices;
            indf = ismember(Surf.SurfData.faces,ind);temp = sum(indf');
            fv.faces = Surf.SurfData.faces(logical(temp),:);
            [At] = Area_Comp(fv);
            area(j,i) = At;
        else
            cth(j,i) = 0;
            stdcth(j,i) = 0;
            area(j,i) = 0;
        end
    end
end
charac = ';';
Names = ['SubjectsIDs' charac 'background' charac 'gyrus14' charac	'Frontal-Superior' charac 'Lingual' charac 'gyrus12' charac 'gyrus19'...
    charac 'gyrus8' charac 'Frontal-Inferior' charac 'Temporal-Superior' charac	'gyrus20' charac 'gyrus10' charac 'Pre-Central'...
    charac '12	Post-Central' charac 'gyrus15' charac	'Frontal-Middle' charac 'Orbital' charac 'Temporal-Inferior' charac 'Temporal-Middle' charac 'Cuneus'];
Mcads = createcads(cth,charac);Mcads = [Ids Mcads];
Name = '/media/LA-PUBLIC/curvature/Results_MeanCTh_Left.txt';
save_resfiles(Name,Names,Mcads);
Name = '/media/LA-PUBLIC/curvature/Results_STDCTh_Left.txt';
Mcads = createcads(stdcth,charac);Mcads = [Ids Mcads];
save_resfiles(Name,Names,Mcads);
Name = '/media/LA-PUBLIC/curvature/Results_Area_Left.txt';
Mcads = createcads(area,charac);Mcads = [Ids Mcads];
save_resfiles(Name,Names,Mcads);






%pth = '/media/COSAS/Test/Joost'
pth = '/media/LA-PUBLIC/curvature/';
gtexfiles = sel_files(pth,'*_Rwhite_gyri_default_session_auto.tex');
ctexfiles = sel_files(pth,'*_Rwhite_curv.tex');
%ctexfiles = sel_files('/media/COSAS/Test/Joost/','*_curv_bary.tex');
thtexfiles = sel_files(pth,'*_R_wm_thick.tex');
% thtexfiles = ctexfiles;
Ns = size(gtexfiles,1);Ids = '';
for j =1:Ns
    %disp(['Computing Subject ']);
    [pth,nam,ext] = fileparts(deblank(ctexfiles(j,:)));
    ind = strfind(nam,'-');
    ID = nam(1:ind+15);
    Ids = strvcat(Ids,ID);
    disp(['Computing Subject:  ==>  ' ID]);
    Sfile = ['/media/COSAS/Test/Joost/parc/temp' filesep ID '.mesh'];
    [Textg,Format] = read_texBrainvisa(deblank(gtexfiles(j,:)));
    [OutFiles, SurfF] = Exp_Surf(Sfile, '0', '','', 'imp','n');
    Surf = SurfF{1};
    Npoints = size(Textg.Values,1);
    Surf.Is = Textg.Values(:,3);
    ind = find(Surf.Is ==0);
    str = unique(Surf.Is);
    Ns = max(str);Surf.Is(ind) = Ns+1;
    [Surf] = Surf_Corr(Surf);
    ind = find(Surf.Is ==Ns+1); Surf.Is(ind) = 0;
    [Textc,Format] = read_texBrainvisa(deblank(ctexfiles(j,:)));
    ThSup = 0.1*pi;
    ThInf = -0.1*pi;
    indvc = find((Textc.Values(:,1)<=ThSup)&(Textc.Values(:,1)>=ThInf));
    
    [Textth,Format] = read_texBrainvisa(deblank(thtexfiles(j,:)));
    
    for i =1:size(str,1);
        if (str(i)~=0)|(str(i)~=0)|(str(i)~=0)
            ind = find(Surf.Is==str(i));
            indi = ismember(ind,indvc);
            if sum(indi)==0
                disp(['Subject:  ==>  ' ID ' BLAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA']);
            end
            cth(j,i) = mean(Textth.Values(ind(indi),1));
            stdcth(j,i) = std(Textth.Values(ind(indi),1));
            fv.vertices = Surf.SurfData.vertices;
            indf = ismember(Surf.SurfData.faces,ind);temp = sum(indf');
            fv.faces = Surf.SurfData.faces(logical(temp),:);
            [At] = Area_Comp(fv);
            area(j,i) = At;
        else
            cth(j,i) = 0;
            stdcth(j,i) = 0;
            area(j,i) = 0;
        end
    end
end
charac = ';';
Names = ['SubjectsIDs' charac 'background' charac 'gyrus14' charac	'Frontal-Superior' charac 'Lingual' charac 'gyrus12' charac 'gyrus19'...
    charac 'gyrus8' charac 'Frontal-Inferior' charac 'Temporal-Superior' charac	'gyrus20' charac 'gyrus10' charac 'Pre-Central'...
    charac '12	Post-Central' charac 'gyrus15' charac	'Frontal-Middle' charac 'Orbital' charac 'Temporal-Inferior' charac 'Temporal-Middle' charac 'Cuneus'];
Mcads = createcads(cth,charac);Mcads = [Ids Mcads];
Name = '/media/LA-PUBLIC/curvature/Results_MeanCTh_Right.txt';
save_resfiles(Name,Names,Mcads);
Name = '/media/LA-PUBLIC/curvature/Results_STDCTh_Right.txt';
Mcads = createcads(stdcth,charac);Mcads = [Ids Mcads];
save_resfiles(Name,Names,Mcads);
Name = '/media/LA-PUBLIC/curvature/Results_Area_Right.txt';
Mcads = createcads(area,charac);Mcads = [Ids Mcads];
save_resfiles(Name,Names,Mcads);

return;


function save_resfiles(Name,Names,Val);
fid = fopen(Name,'wt');
Nf = size(Val,1);
fprintf(fid,'%s\n',Names);
for i = 1:Nf
    fprintf(fid,'%s\n',deblank(Val(i,:)));
end
return

function Mcads = createcads(Mat,charact);
[m,n] = size(Mat);
Mcads = '';
for i = 1:m
    ncad = '';
    for j = 1:n
        ncad= [ncad charact num2str(Mat(i,j))];
    end
    Mcads = strvcat(Mcads,ncad);
end
return