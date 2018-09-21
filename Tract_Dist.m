function Tract_Dist(Surf)
j = 1;
P1 = [ 0 0 0]; P2 = [ 0 0 0]; P3 = [ 0 0 0]; indexs =0;
for j = 1:length(Surf.SurfData.lines)
    Vert = Surf.SurfData.vertices(Surf.SurfData.lines{j},:);
    P1 = [P1;Vert(1:size(Vert,1)-2,:)];
    P2 = [P2;Vert(2:size(Vert,1)-1,:)];
    P3 = [P3;Vert(3:size(Vert,1),:)];
    indexs(j+1) = size(P1,1)-1;
end
P1(1,:) = [];P2(1,:) = [];P3(1,:) = [];

n = cross(P1-P2,P3-P2,2); D = -1*dot(n,P1,2);
t = P2-P1; u = P3-P1; v = P3-P2;
t2 = sum((t.^2)')'; u2 = sum((u.^2)')'; w2 = sum((n.^2)')';
c = P1+(repmat(t2.*dot(u,v,2),[1 size(u,2)]).*u-repmat(u2.*dot(t,v,2),[1 size(t,2)]).*t)./(2*repmat(w2,[1 size(P1,2)]));
r = 1/2*sqrt(t2.*u2.*dot(v,v,2)./w2);
n = cross(P1-c,P2-c,2); D = -1*dot(n,c,2);
n1 = cross(P2-c,n,2); D1 = -1*dot(n1,c,2);

for i = 1:size(n1,1)
    ind = min(find(indexs-i)>0));
    ind = find(indexs-i == min(indexs-i)&(indexs-i)>0);
    P1temp = P1;
    P2temp = P2;
    
    N = size(P1temp,1);
     Num = dot(repmat(n1(i,:),[N 1]),P1temp,2)+repmat(D1(i,:),[N 1]); 
     Den = dot(repmat(n1(i,:),[N 1]),P2temp-P1temp,2);
     t = -1*Num./Den;clear Num Den;
     intersc = single(P1temp)+repmat(t,[1 3]).*(P2temp-P1temp);
     dots = dot(P1temp-intersc,intersc-P2temp,2)./(sqrt(dot(P1temp-intersc,P1temp-intersc,2)).*sqrt(dot(intersc-P2temp,intersc-P2temp,2))); ind = find(dots ==1);
     points = intersc(ind,:);
     clear t;
end

