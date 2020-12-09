function disp_atmos_ifg(a,area_pt,ifg_pt,ifg_x,ifg_y,bar_str,clim)

figure
plot(area_pt(:,1),area_pt(:,2),'b-','linewidth',2);
hold on
plot(ifg_pt(:,1),ifg_pt(:,2),'b','linewidth',2);
axis image

if ~isempty(clim)
    scatter(ifg_x(:),ifg_y(:),5*ones(length(ifg_y(:)),1),a(:),'d','filled');
    caxis(clim)
else
    scatter(ifg_x(:),ifg_y(:),5*ones(length(ifg_y(:)),1),a(:),'d','filled');
end

colorbar
Label_format('Y (km)','X (km)');
axis image

% [n0,m0]=size(a);
m0=max(area_pt(:,1));
n0=max(area_pt(:,2));

m=4;n=4;
set(gca,'Xtick',linspace(0,m0,m+1),'Ytick',linspace(0,n0,n+1));
xstr={};ystr={};
dx=(m0/m/1000);
dy=(n0/n/1000);
% xlim([0,m0]);
% ylim([0,n0]);
for i=1:m+1
    tmp=strcat(num2str(0+dx*(i-1),'%.1f'));
    xstr=[xstr,tmp];
end
for i=1:n+1
    tmp=strcat(num2str(0+dy*(i-1),'%.1f'));
    ystr=[ystr,tmp];
end
set(gca,'Xticklabel',xstr,'Yticklabel',ystr,'YTickLabelRotation',90);
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.6]);

theta1=atan((ifg_pt(2,2)-ifg_pt(1,2))/(ifg_pt(2,1)-ifg_pt(1,1))*4/3);
dx1=0.1*cos(theta1);
dy1=0.1*sin(theta1);
theta2=atan((ifg_pt(4,2)-ifg_pt(1,2))/(ifg_pt(4,1)-ifg_pt(1,1))*4/3);
dx2=0.1*cos(theta2);
dy2=0.1*sin(theta2);

if (ifg_pt(2,1)-ifg_pt(1,1))<0
    theta_azi=theta1+pi;
else
    if (ifg_pt(2,2)-ifg_pt(1,2))<0
        theta_azi=theta1+2*pi;
    else
        theta_azi=theta1;
    end
end

if theta_azi<pi/2
    annotation('textarrow',[0.32 0.32+dx1],[0.22 0.22+dy1],'String','Azimuth','fontsize',16);
    annotation('textarrow',[0.25 0.25-dx2],[0.3 0.3-dy2],'String','Range','fontsize',16);
elseif theta_azi>pi/2 && theta_azi<pi
    annotation('textarrow',[0.62 0.62-dx1],[0.32 0.32-dy1],'String','Azimuth','fontsize',16);
    annotation('textarrow',[0.62 0.62-dx2],[0.22 0.22-dy2],'String','Range','fontsize',16);
elseif theta_azi>pi && theta_azi<pi*3/2
    annotation('textarrow',[0.55 0.55-dx1],[0.88 0.88-dy1],'String','Azimuth','fontsize',16);
    annotation('textarrow',[0.68 0.68+dx2],[0.78 0.78+dy2],'String','Range','fontsize',16);
else
    annotation('textarrow',[0.31 0.31+dx1],[0.73 0.73+dy1],'String','Azimuth','fontsize',16);
    annotation('textarrow',[0.32 0.32+dx2],[0.82 0.82+dy2],'String','Range','fontsize',16);
end

c=colorbar('FontSize',18,'Linewidth',1.5);
c.Label.String = bar_str;
c.FontSize=18;
end