function Perp_planes_Comp(Surf)

P1 = [ 0 0 0]; P2 = [ 0 0 0]; P3 = [ 0 0 0]; indexs = 0;
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
    ind = min(find(indexs-i>0));
    P1temp = P1;
    %P1temp(indexs(ind-1)+1:indexs(ind),:)=[];
    P2temp = P2;
    %P2temp(indexs(ind-1)+1:indexs(ind),:)=[];
    N = size(P1temp,1);
     Num = dot(repmat(n1(i,:),[N 1]),P1temp,2)+repmat(D1(i,:),[N 1]); 
     Den = dot(repmat(n1(i,:),[N 1]),P2temp-P1temp,2);
     t = -1*Num./Den;clear Num Den;
     intersc = single(P1temp)+repmat(t,[1 3]).*(P2temp-P1temp);
     dots = dot(P1temp-intersc,intersc-P2temp,2)./(sqrt(dot(P1temp-intersc,P1temp-intersc,2)).*sqrt(dot(intersc-P2temp,intersc-P2temp,2))); ind = find(round(dots) ==1);
     
     
     temp = zeros(size(n1),1);temp(indexs(2:end)) = 2;temp(ind) =1;ind2 = find(temp);temp2 = temp(ind2);a = diff(temp2);  ind3 = find((a ==0)&temp2(1:end-1)~=2); indexes = [ind2(ind3) ind2(ind3+1)];
     a = repmat(P2temp(i,:),[size(indexes(:),1) 1])-intersc(indexes(:),:);
     distance = sqrt(sum((a.^2)'))';distemp = reshape(distance,[size(indexes(:),1)/2 2]);
     [a1,b1] = sort(distemp');
     intee = 1:size(b1,2); indtemps = sub2ind(size(indexes), intee',b1(2,:)');indsdef = indexes(indtemps);indr = ismember(ind,indsdef);ind(indr) = [];
     intersc = intersc(ind,:);
     a = repmat(P2temp(i,:),[size(ind(:),1) 1])-intersc;distance = sqrt(sum((a.^2)'))';temp1 = diff(sort(distance)); indr = find(temp1 == max(temp1));
     if indr ~=length(temp1)
         th = (temp1(indr)+temp1(indr+1))/2;
     else
         th = (temp1(indr))/2;
     end
     
     ind = find(distance>=th); intersc(ind,:) = [];
     hold on;plot3(intersc(:,1), intersc(:,2),intersc(:,3),'.b','Markersize',20);
end









%Plot_Normals(Vert(2:end-1,:), -1*n1);
% for i = 1:size(n1,1)
%     if ~isnan(sum(n1(i,:)))
%         N = repmat([n1(i,:) D1(i)],[size(X,1) 1]);
%         Z = ((-N(:,4)-N(:,1).*X-N(:,2).*Y))./N(:,3);
%     else
%         Z = (i+scont)*ones(size(X,1),1);
%     end
end