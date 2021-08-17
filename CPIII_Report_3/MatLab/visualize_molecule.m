function  d = visualize_molecule(coord,thr,figure_num,molecule_name)
%VISUALIZE_MOLECULE generates a 3D representation of a system of points.
%
%Invoked as d=visualize_molecule(coord,thr,molecule_name), coord is the
%matrix of atomic coordinates, thr is the largest atomic distance to draw a
%bond and molecule_name is a string containing the title of the plot.
%
%Coord has to be a Nx3 matrix with the x,y,z coordinates of each point
%stored as rows. The function generates a figure, save it as an eps file and gives
%the distance map of the points d

%Search for connectivities
n=size(coord,1);

d=zeros(n,n);

nconn=0;

for i1=1:n
    for i2=i1+1:n
        d(i1,i2)=sum((coord(i1,:)-coord(i2,:)).^2)^(0.5);
        
        if d(i1,i2) < thr
            nconn=nconn+1;
            conn(nconn,1)=i1; conn(nconn,2)=i2;
        end
    end
end

d=d+d';

figure(figure_num);

hold on;

for i1=1:n
     plot3(coord(i1,1),coord(i1,2),coord(i1,3),'-r.','MarkerSize',30);
end

for i1=1:nconn
    line([coord(conn(i1,1),1) coord(conn(i1,2),1)],[coord(conn(i1,1),2) coord(conn(i1,2),2)],[coord(conn(i1,1),3) coord(conn(i1,2),3)],'LineWidth',2);
end

axis equal;
grid on;

xlim([min(coord(:,1))-1 max(coord(:,1))+1]);
ylim([min(coord(:,2))-1 max(coord(:,2))+1]);
zlim([min(coord(:,3))-1 max(coord(:,3))+1]);

xlabel('X (sigma)');
ylabel('Y (sigma)');
zlabel('Z (sigma)');

title(molecule_name);

view([-60 26]);

filename=strcat(molecule_name,'.eps');
print(gcf,'-depsc',filename)

hold off;


end
