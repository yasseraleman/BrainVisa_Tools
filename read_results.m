Filename = '/media/COSAS/Test/Joost/allvars2.csv'
Ns = 102; %Number of subjects;
deli = ',';
a = textread(Filename,'%s','delimiter',deli);
for i = 1:size(a,1)
    if strcmp(char(a{i}),'" "')
        a{i} = 'NaN';
    end
end
b = reshape(a,[size(a,1)/(Ns+1) Ns+1]);b = b';
% bz = zeros(size(b,1),size(b,2));
% for i =2:size(b,1)
%     for j = 4:size(b,2)
%         bz(i-1,j-1) = str2num(b{i,j});
%     end
% end
fid = fopen(Filename,'rt');
line = fgetl(fid);
names = strread(line,'%s','delimiter',deli);
nam = char(names);post =0;it = 0;news = ''; newst = '';
test = char(b);
for i = 1:size(names,1)
    namt = char(names{i});
    indp = strfind(namt,'_L_');
    if ~isempty(indp)
        newname = [namt(1:indp) 'R_' namt(indp+3:end)];
        ind = ismember(nam,newname,'rows');
        pos = find(ind);
        if pos
            post = [post;i;pos];
            left = test([(i-1)*(Ns+1)+1:i*(Ns+1)],:);
            right = test([(pos-1)*(Ns+1)+1:pos*(Ns+1)],:);
            suma = str2num(left(2:end,:))+str2num(right(2:end,:));
            sumaname = [namt(1:indp) 'SUM_' namt(indp+3:end)];
            means = (str2num(left(2:end,:))+str2num(right(2:end,:)))/2;
            meansname = [namt(1:indp) 'MEAN_' namt(indp+3:end)];
            news = strvcat(news,sumaname,num2str(suma),meansname,num2str(means));
        else
            it = [it;i];
            newst = strvcat(newst,test([(i-1)*(Ns+1)+1:i*(Ns+1)],:));
        end
    else
        indt = strfind(namt,'_R_');
        if isempty(indt)
            it = [it;i];
            newst = strvcat(newst,test([(i-1)*(Ns+1)+1:i*(Ns+1)],:));
        else
            a = 1;
        end
        
    end
end
newst = strvcat(newst,news);
newstotal = newst;
for j = 1:size(newst,1)
    indp = strfind(newst(j,:),'NaN');
    if ~isempty(indp)
        if indp~=1
            newst(j,:) = [newst(j,1:indp-1) '" "' newst(j,indp+3:end)];
        else
            newst(j,:) = ['" "' newst(j,indp+3:end)];
        end
    end
end
nvar = size(newst,1)/(Ns+1);
[pth,nm,ext] = fileparts(Filename);
fid = fopen([pth filesep nm '_Mean-Sum.csv' ],'wt');
 for i = 1:Ns+1
    cad = '';
    for j = 1:nvar
        cad = [cad deli deblank(newst((j-1)*(Ns+1)+i,:))];
    end
    cad(1) = [];
    fprintf(fid,'%s\n',cad);
 end
fclose(fid);
