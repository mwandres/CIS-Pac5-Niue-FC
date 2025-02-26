%% Determine len in ccw order for UNSWAN
clearvars; clc;
%addpath(genpath('D:\Tools\OceanMesh2D\utilities\'))
%addpath('D:\Projects_SPC\CK_GCF\CK_hindcast\Soft\Build_UNSWAN_boundaries\codes')

%% load the positions of the hidcast points to be used as forcing
d=shaperead('NOAA_GFS_points_Rarotonga.shp');
xbnd=[];ybnd=[];
for i=1:length(d);
    xbnd=[xbnd;d(i).X];
    ybnd=[ybnd;d(i).Y];
end
% clear d;
nc_tmp_file = 'C:\Users\moritzw\OneDrive - SPC\Documents\Projects\ECIKS\Niue\Forecast\CIS-PAC5-Niue\operational_v1\tmp/wave_tmp_000.nc';
lat = ncread(nc_tmp_file,'lat');
lon = ncread(nc_tmp_file,'lon');
[lonn,latt] = meshgrid(lon,lat);
xbnd = double([lonn(1:1:end-1,1); lonn(end,1:1:end-1)'; lonn(end:-1:2,end); lonn(1,end:-1:2)']);
ybnd = double([latt(1:1:end-1,1); latt(end,1:1:end-1)'; latt(end:-1:2,end); latt(1,end:-1:2)']);
xbnd = xbnd-360;
figure()
hold on
xlim([min(xbnd)-1 max(xbnd)+1])
ylim([min(ybnd)-1 max(ybnd)+1])
for i = 1:length(xbnd)
    text(xbnd(i), ybnd(i),num2str(i))
end

poly = [xbnd(1),ybnd(1)];
xbnd(1) = nan;
ybnd(1) = nan;
for i =2:length(xbnd);
    did= sqrt((xbnd-poly(end,1)).^2+(ybnd-poly(end,2)).^2);
    [~,posn] = min(did);
    poly = [poly;[xbnd(posn),ybnd(posn)]];
    xbnd(posn) = nan;
    ybnd(posn) = nan;
    
end
    
xbnd = poly(:,1);
ybnd = poly(:,2);


[xbnd, ybnd] = poly2ccw(xbnd, ybnd)

% 


%% load the ADCIRC mesh
[fem,bnd]=read_adcirc_mesh('C:\Users\moritzw\OneDrive - SPC\Documents\Projects\ECIKS\Niue\Forecast\CIS-PAC5-Niue\operational_v1\extras\niue\swan\fort.14');
x1=fem.x(bnd.elev_nodes);
y1=fem.y(bnd.elev_nodes);
[x1, y1] = poly2ccw(x1, y1);%% order the open boundary in counterclockwise direction
 

% figure
% plot(xbnd,ybnd,'.-')
% hold on;
% plot3(x1,y1,repmat(10,numel(x1)),'*')
% hold on;
% plot(x1(1),y1(1),'o')
% plot(x1(posn),y1(posn),'o')
% plot(x1(2),y1(2),'*')

%% look for the closer point to the origin [0,0]
% [~,m] = min(x1.^2+y1.^2);
% did= sqrt((x1-x1(m)).^2+(y1-y1(m)).^2);
% did(m)=nan;
% [~,posn] = min(did);
[~,posn] = min(x1+y1);

% posn=13%% Still no idea why 25 is the first point of the boundary


% 
% Dist=[];
% for i=1:length(x1);
% %     [d1km d2km]=lldistkm([y1(i),x1(i)],[0,0]);
%     d2km=y1(i).^2+x1(i).^2;
%     Dist=[Dist;d2km];
% end
% 
% [kx,posn]=min(Dist)

figure;
plot(x1,y1,'.-')
hold on;
% plot(x1(1),y1(1),'o')
plot(x1(posn),y1(posn),'*')
% 
x1=[x1(posn:end-1);x1(1:posn-1);x1(posn)];
y1=[y1(posn:end-1);y1(1:posn-1);y1(posn)];

%% look for the closest open boundary point

distances=[];
for i=1:length(x1)
    xx=x1(1:i);
    yy=y1(1:i);
    len=0;
    for j=1:i-1;
        len=len+sqrt((xx(j)-xx(j+1)).^2+(yy(j)-yy(j+1)).^2);
    end
    distances=[distances;len];
end



did= sqrt((xbnd-x1(1)).^2+(ybnd-y1(1)).^2);
[len,pose] = min(did);

plot(xbnd,ybnd,'o')
plot(xbnd(pose),ybnd(pose),'ok')
% plot(xbnd(pose+1),ybnd(pose+1),'or')
% plot(xbnd(pose+2),ybnd(pose+2),'og')


xbnd=[xbnd(pose:end);xbnd(1:pose-1)];
ybnd=[ybnd(pose:end);ybnd(1:pose-1)];

figure
plot(xbnd,ybnd,'-.')
hold on;
plot(xbnd(1),ybnd(1),'ok')
plot(xbnd(2),ybnd(2),'or')
plot(xbnd(3),ybnd(3),'og')
plot(xbnd(10),ybnd(10),'om')

%% order the open boundary in counterclockwise direction
distancesE=[];
for i=1:length(xbnd)
%     xx=xbnd(1:i);
%     yy=ybnd(1:i);
    %     len=0;
    if i == 1
        distancesE=[distancesE;len];
    else
%         for j=2:i-1;
            len=len+sqrt((xbnd(i-1)-xbnd(i)).^2+(ybnd(i-1)-ybnd(i)).^2);
%         end
        distancesE=[distancesE;len];
    end
end

distancesE = round(distancesE*100)./100 ;
% 
% 
% 
% xbnew=ones(length(xbnd),1);
% ybnew=xbnew;Len=[];
% for i=1:length(xbnd);
%     [kk,posn]=min(sqrt((x1-xbnd(i)).^2+(y1-ybnd(i)).^2));%% search the closest point of the boundary to each of the hindcast points
%     xbnew(i)=x1(posn);ybnew(i)=y1(posn);%% store the new position to exactly match with one point of the boundary
%     %% Determine the distance of each point respect the first point of the open boundary, as the boundary was rearranged in ccw direction the distances are too
%     xx=x1(1:posn);
%     yy=y1(1:posn);
%     len=0;
%     for j=1:posn-1;
%         len=len+sqrt((xx(j)-xx(j+1)).^2+(yy(j)-yy(j+1)).^2);
%     end
%     
%     Len=[Len;len];
% end
% 
% %% ordering the hindcast points regarding the distance from the first point of the open boundary
% [~,is] = sort(Len);
% 
% xaux=xbnew(is);
% yaux=ybnew(is);
% Len=Len(is);


figure;
hold on;
plot(x1,y1,'.-')
plot(x1(1),y1(1),'ok')
plot(xbnd,ybnd,'og')
% plot(x1(10),y1(10),'or')
% plot(xaux,yaux,'*')
for i=1:length(xbnd);
     text(xbnd(i),ybnd(i), [num2str(i),' ,',num2str(round(distancesE(i),1))],'FontName', 'Bookman Old Style','fontsize',8,'fontweight','bold');
end
grid on;
box on;
grid on; box on;
xlabel('lonÂº');
ylabel('latÂº');
daspect([1 1 1])
legend('open boundary','1st point of the open boundary', 'original hindcast locations');
set(gca, 'FontName', 'Bookman Old Style','fontsize',10);
% file=['D:\Projects_SPC\Kiribati\Figures\BoundaryPointsandDistance'];
% print('-dpng','-r300',file)


data=[xbnd,ybnd,distancesE];
save D:\Projects_SPC\CK_GCF\CK_hindcast\meshes\Rarotonga\LEN_POS_boundary_Raro_operational.txt data -ascii


index = [1:length(data)]';

% Create a table with time and sea level data
Data = table(index, xbnd,ybnd, 'VariableNames', {'index', 'x', 'y'});
% Specify the file name
filename = 'D:\Projects_SPC\CK_GCF\CK_hindcast\meshes\Rarotonga\f14_boundary_points.csv';
% Write the table to a CSV file
writetable(Data, filename);