function Plot3_tracts(Surf)
figure;
for j = 1:100%length(Surf.SurfData.lines)
    Vert = Surf.SurfData.vertices(Surf.SurfData.lines{j},:);
    Faces = [(1:size(Vert,1)-2)' (2:size(Vert,1)-1)'];
    p1 = Vert(Faces(:,1),:);
    p2 = Vert(Faces(:,2),:);
    t = p2-p1;
    t2 = sum((t.^2)')';
    distance = 0;
    for i = 2:size(t2,1)+1
        distance(i) = distance(i-1)+t2(i-1);
    end
    ind = Surf.SurfData.lines{j};%ind = ind(2:end);
    hold on;
    x = ones(size(ind,1),1)*j;
    distance(end+1) = distance(end)+min(t2); y = distance'; z = Surf.Is(ind);z(end) = 0; z(1) = 0;
    plot3(x,y,z,'-r')
    
    

end
return;