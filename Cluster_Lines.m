function Cluster_Lines(Surf)
N = 100;
P1 = [ 0 0 0]; P2 = [ 0 0 0]; P3 = [ 0 0 0]; indexs = 0;
Mat = zeros(length(Surf.SurfData.lines),length(Surf.SurfData.lines));figure;
for j = 1:length(Surf.SurfData.lines)
    Vertj = Surf.SurfData.vertices(Surf.SurfData.lines{j},:);
    for i = j+1:length(Surf.SurfData.lines)
        Verti = Surf.SurfData.vertices(Surf.SurfData.lines{i},:);
        index = repmat([1:size(Verti,1)],[size(Vertj,1) 1]);
        index = index(:);
        a = repmat(Vertj,[size(Verti,1) 1])-Verti(index,:);
        distance = sqrt(sum((a.^2)'))';
        [a,v] = sort(distance);val = mean(distance(v(1:N)));
         %aa = reshape(distance,[size(Vertj,1) size(Verti,1)]);
        %imagesc(aa);drawnow;
%         line(Surf.SurfData.vertices(Surf.SurfData.lines{j},1),Surf.SurfData.vertices(Surf.SurfData.lines{j},2),  Surf.SurfData.vertices(Surf.SurfData.lines{j},3),'Color',[0 1 0],'Linewidth',5);
%         line(Surf.SurfData.vertices(Surf.SurfData.lines{i},1),Surf.SurfData.vertices(Surf.SurfData.lines{i},2),  Surf.SurfData.vertices(Surf.SurfData.lines{i},3),'Color',[1 0 0],'Linewidth',5);
%         bb = repmat(Vertj,[size(Verti,1) 1]);
%         bb = bb(v(1:N),:); hold on;plot3(bb(:,1),bb(:,2),bb(:,3),'.b','Markersize',40)
%         bb = Verti(index,:);
%         bb = bb(v(1:N),:); hold on;plot3(bb(:,1),bb(:,2),bb(:,3),'.y','Markersize',40);drawnow
%         pause(0.3)
%         cla(gca);
        Mat(j,i) = val;
        Mat(i,j) = val;
    end
    
end
a = 1;
