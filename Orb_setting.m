function Orb_setting(input_path)
if nargin<1
    clear,clc
    fclose all;
    close all;
    input_path='D:/Study/TUD_Work/data2/N_delay_out/';
end
% inc_min=23*pi/180;
Rm=6372593.6944;
near_ra=561544.5107;
R2sar=6884407.4799;%For TSX
% Rm=6373480.462756;
% near_ra=901073.190;
% R2sar=7072396.782;%For Sentinel
inc_min = acos((R2sar^2 +near_ra.^2 - Rm^2)./ (2 * R2sar * near_ra));

azi_angle=120*pi/180;
ifg_x_sp=25;
ifg_y_sp=25;

xyz_file=[input_path,'N_xyz.nc'];
xt=ncread(xyz_file,'xt');
yt=ncread(xyz_file,'yt');
zt=ncread(xyz_file,'zt');
t=ncread(xyz_file,'time');
nx=length(xt);
ny=length(yt);
nt=length(t);
%model resolution
atmos_x_sp=xt(2)-xt(1);
atmos_y_sp=yt(2)-yt(1);
atmos_x_ed=nx*atmos_x_sp;
atmos_y_ed=ny*atmos_y_sp;

dh=[diff(zt)];
dh=[dh;dh(end)];


nz=length(zt);
z_max=zt(nz);

%Range offset due to the slant view
Rg_off=ifg_x_sp*(ceil(z_max*tan(inc_min)/ (ifg_x_sp)));

%ifg resampling grid based on the azimuth angle
if azi_angle>=0 && azi_angle<=pi/2
    heading=azi_angle;
    if azi_angle==pi/2
        dx=0;
    else
        dx=atmos_x_ed/(1+tan(heading));
    end
    model_pt=[atmos_x_ed-dx,0;
             atmos_x_ed,atmos_y_ed-dx;
             dx,atmos_y_ed;
             0,dx;
             atmos_x_ed-dx,0];
         
     dx_off=Rg_off*sin(heading);
     dy_off=Rg_off*cos(heading);
     
     ifg_pt=[atmos_x_ed-dx-dx_off,0+dy_off;
             atmos_x_ed-dx_off,atmos_y_ed-dx+dy_off;
             dx,atmos_y_ed;
             0,dx;
             atmos_x_ed-dx-dx_off,0+dy_off];
         
elseif azi_angle >pi/2 && azi_angle <=pi
    heading=azi_angle-pi/2;
    dx=atmos_x_ed/(1+tan(heading));
    model_pt=[atmos_x_ed,atmos_y_ed-dx;
             dx,atmos_y_ed;
             0,dx;
             atmos_x_ed-dx,0;
             atmos_x_ed,atmos_y_ed-dx];
         
     dx_off=Rg_off*cos(heading);
     dy_off=Rg_off*sin(heading);
     ifg_pt=[atmos_x_ed-dx_off,atmos_y_ed-dx-dy_off;
             dx-dx_off,atmos_y_ed-dy_off;
             0,dx;
             atmos_x_ed-dx,0;
             atmos_x_ed-dx_off,atmos_y_ed-dx-dy_off];
elseif azi_angle>pi && azi_angle <=3*pi/2
    heading=azi_angle-pi;
    if azi_angle==3*pi/2
        dx=0;
    else
        dx=atmos_x_ed/(1+tan(heading));
    end
    model_pt=[dx,atmos_y_ed;
             0,dx;
             atmos_x_ed-dx,0;
             atmos_x_ed,atmos_y_ed-dx;
             dx,atmos_y_ed];
         
     dx_off=Rg_off*sin(heading);
     dy_off=Rg_off*cos(heading);
     ifg_pt=[dx+dx_off,atmos_y_ed-dy_off;
             0+dx_off,dx-dy_off;
             atmos_x_ed-dx,0;
             atmos_x_ed,atmos_y_ed-dx;
             dx+dx_off,atmos_y_ed-dy_off];
else
    heading=azi_angle-3*pi/2;
    dx=atmos_x_ed/(1+tan(heading));
    model_pt=[0,dx;
             atmos_x_ed-dx,0;
             atmos_x_ed,atmos_y_ed-dx;
             dx,atmos_y_ed;
             0,dx];
         
     dx_off=Rg_off*cos(heading);
     dy_off=Rg_off*sin(heading);
     ifg_pt=[0+dx_off,dx+dy_off;
             atmos_x_ed-dx+dx_off,0+dy_off;
             atmos_x_ed,atmos_y_ed-dx;
             dx,atmos_y_ed;
             0+dx_off,dx+dy_off];
end

model_range_dis=distance_2D(model_pt(1,:),model_pt(4,:));
model_azi_dis=distance_2D(model_pt(1,:),model_pt(2,:));
model_rg=floor(model_range_dis/ifg_y_sp);
model_azi=floor(model_azi_dis/ifg_x_sp);

ifg_range_dis=distance_2D(ifg_pt(1,:),ifg_pt(4,:));
ifg_azi_dis=distance_2D(ifg_pt(1,:),ifg_pt(2,:));
ifg_rg=floor(ifg_range_dis/ifg_y_sp);
ifg_azi=floor(ifg_azi_dis/ifg_x_sp);
         
model_pt_img=[0,0;
    model_azi*ifg_x_sp,0;
    model_azi*ifg_x_sp,model_rg*ifg_y_sp;...
    0,model_rg*ifg_y_sp];
geo_par=trans_6par(model_pt_img,model_pt(1:4,:));

ifg_pt_img=[0,0;
    ifg_azi*ifg_x_sp,0;
    ifg_azi*ifg_x_sp,ifg_rg*ifg_y_sp;...
    0,ifg_rg*ifg_y_sp];
ifg_geo_par=trans_6par(ifg_pt_img,ifg_pt(1:4,:));

far_ra=near_ra+ifg_x_sp*ifg_rg;
inc_max = acos((R2sar^2 +far_ra.^2 - Rm^2)./ (2 * R2sar * far_ra));

% d_inc=10*pi/180;
% inc_max=inc_min+d_inc;

inc_max0=atan((model_range_dis-ifg_y_sp*2)/z_max);
if inc_max>inc_max0
    inc_max=inc_max0;
end
% 
inc=linspace(inc_min,inc_max,ifg_rg);

orb_slc.model_rg=model_rg;
orb_slc.model_azi=model_azi;
orb_slc.model_pt=model_pt;

orb_slc.ifg_y_sp=ifg_y_sp;
orb_slc.ifg_x_sp=ifg_x_sp;
orb_slc.ifg_rg=ifg_rg;
orb_slc.ifg_azi=ifg_azi;
orb_slc.ifg_pt=ifg_pt;
orb_slc.geo_par=geo_par;
orb_slc.ifg_geo_par=ifg_geo_par;
orb_slc.lam=0.056;

orb_slc.z_max=z_max;
orb_slc.xt=xt;
orb_slc.yt=yt;
orb_slc.zt=zt;
orb_slc.dh=dh;
orb_slc.t=t;
orb_slc.inc=inc;
% 
figure
plot(model_pt(:,1),model_pt(:,2),'b','linewidth',2);
for i=1:size(model_pt,1)-1
    text(model_pt(i,1),model_pt(i,2),num2str(i),'color','b');
end
axis image
hold on
plot(ifg_pt(:,1),ifg_pt(:,2),'r','linewidth',1);
for i=1:size(ifg_pt,1)-1
    text(ifg_pt(i,1),ifg_pt(i,2),num2str(i),'color','r');
end
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
    annotation('textarrow',[0.65 0.65-dx1],[0.32 0.32-dy1],'String','Azimuth','fontsize',16);
    annotation('textarrow',[0.62 0.62-dx2],[0.22 0.22-dy2],'String','Range','fontsize',16);
elseif theta_azi>pi && theta_azi<pi*3/2
    annotation('textarrow',[0.55 0.55-dx1],[0.88 0.88-dy1],'String','Azimuth','fontsize',16);
    annotation('textarrow',[0.68 0.68+dx2],[0.78 0.78+dy2],'String','Range','fontsize',16);
else
    annotation('textarrow',[0.35 0.35+dx1],[0.73 0.73+dy1],'String','Azimuth','fontsize',16);
    annotation('textarrow',[0.35 0.35+dx2],[0.85 0.85+dy2],'String','Range','fontsize',16);
end

saveas(gcf,strcat(input_path,'/SAR_geo.tif'));
close(gcf);

save([input_path,'SAR_geo.mat'],'orb_slc');
end