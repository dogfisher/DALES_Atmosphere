function process_ncfiles(input_path)
%merge all filed*.nc files and divide the file into several sub files using
%epoch.
%Write by Fengming Hu
%01.12.2019 in TJU
%input_path is the location of the dataset. The output of DALES model
%contains profiles.nc (pressure value of each thickness) and field*.nc
%files (coordinate x, y, z,(xt,yt,zt) acquired time (t), liquid water potential
%temperature(thl), specific humidity (qt), liquid water specific
%humidity(ql)).
if nargin<1
    clear,clc
    fclose all;
    close all;
    input_path='D:/Study/TUD_Work/data2/';
end
%After run the script, all dataset are saved in output_path
output_path=[input_path,'N_delay_out'];
if exist(output_path,'dir')
    delete([output_path,'/*']);
else
    mkdir(output_path);
end

output_path=[output_path,'/'];

profile_file0=dir([input_path,'profiles*.nc']);
profile_file=[input_path,profile_file0.name];
% ncdisp(profile_file);
para_str1={'presh','time','zt'};
Ps0=(ncread(profile_file,para_str1{1}));
t_epoch=(ncread(profile_file,para_str1{2}));

ncfile_info=ncinfo(profile_file,para_str1{2});
[~,t_unit,~]=ncfile_info.Attributes.Value;
%Check the unit of the time, convert the value to sec.
if ~strcmp(t_unit,'s')
    t_epoch=round(t_epoch*3600*24);
end
zt0=double(ncread(profile_file,para_str1{3}));

%Get all file name of field*.nc simulated using different CPU 
atmos_file=dir([input_path,'field*.nc']);    
ncore=length(atmos_file);
para_str2={'xt','yt','time','zt','qt','ql','thl'};
xt0=[];yt0=[];%sampling along x and y direction
x_n=0;y_n=0;
core_x=0;core_y=0;%number of blocks along x and y direction
for i=1:ncore
    file_path=[input_path,atmos_file(i).name];
    data=double(ncread(file_path,para_str2{1}));
    xt0=unique([xt0;data]);
    if length(xt0)>x_n
        x_n=length(xt0);
        core_x=core_x+1;
    end
    data=double(ncread(file_path,para_str2{2}));
    yt0=unique([yt0;data]);
    if length(yt0)>y_n
        y_n=length(yt0);
        core_y=core_y+1;
    end
end
ii=0;
%merge all coordinate x and y in different blocks
for i=1:core_y:ncore
    ii=ii+1;
    file_path=[input_path,atmos_file(i).name];
    data=(ncread(file_path,para_str2{1}));
    n=length(data);
    xt(ii*n-n+1:ii*n)=data;
end
ii=0;
for i=1:core_y
    ii=ii+1;
    file_path=[input_path,atmos_file(i).name];
    data=(ncread(file_path,para_str2{2}));
    n=length(data);
    yt(ii*n-n+1:ii*n)=data;
end
t=(ncread(file_path,para_str2{3}));

ncfile_info=ncinfo(file_path,para_str2{3});
[~,t_unit,~]=ncfile_info.Attributes.Value;
if ~strcmp(t_unit,'s')
    t=round(t*3600*24);
end
%find the common time samples between field*.nc and profiles.nc
[~,t_idx]=intersect(t_epoch,t);

zt=double(ncread(file_path,para_str2{4}));
[~,z_idx]=intersect(zt0,zt);
Ps=Ps0(z_idx,t_idx);

nx=length(xt);
ny=length(yt);
nz=length(zt);
nt=length(t);

%parameters
dx=nx/core_x;
dy=ny/core_y;

%parameters for atmosphere
Rv = 461.495;
Rd = 287.0586;
eps=Rd/Rv;
k1 = 77.6e-2;
k2 = 70.4e-2;
k3 = 3.739e+3;

Lv=2.25e6;
Cd=1004;
p0=1e5;

%generate the output nc files.
for j=1:nt
    if ~exist([output_path,'N_epoch',num2str(j,'%.3d'),'.nc'],'file')
        nccreate([output_path,'N_epoch',num2str(j,'%.3d'),'.nc'],'N_delay',...
            'Datatype','single','Dimensions',{'dim1',nx,'dim2',ny,'dim3',nz});
    end
end

for i=1:ncore
    fprintf(['output refractivity block',num2str(i),' \n']);
    file_path=[input_path,atmos_file(i).name];

    i_x=ceil(i/core_y);
    i_y=mod(i-1,core_y);

    for j=1:nt
        qt=(ncread(file_path,para_str2{5},[1,1,1,j],[inf,inf,inf,1]));% specific humidity
        ql=(ncread(file_path,para_str2{6},[1,1,1,j],[inf,inf,inf,1]));% liquid water specific humidity
        Tl=(ncread(file_path,para_str2{7},[1,1,1,j],[inf,inf,inf,1])); %liquid water potential temperature
        wv=1./(1-(qt-ql))-1;% mixing ratio
%         wl=1./(1-ql)-1; % liquid water mixing ratio

        C_possion= 0.2854*(1-0.24*wv); % possion constant
        N_delay=zeros(size(Tl));
        for k=1:nz
            e=wv(:,:,k)./(wv(:,:,k)+eps)*Ps(k,j); % partial water vapor pressure
            T0=Tl(:,:,k)./((p0./Ps(k,j)).^C_possion(:,:,k))+Lv/Cd*ql(:,:,k); %surface temperature
            N_delay(:,:,k)=(k1*Ps(k,j)./T0+(k2-eps*k1)*e./T0+k3*e./T0.^2); 
        end
        %save the refractivity
        ncwrite([output_path,'N_epoch',num2str(j,'%.3d'),'.nc'],'N_delay',N_delay,[dx*i_x-dx+1,dy*i_y+1,1]);
    end
end

if ~exist([output_path,'N_xyz.nc'],'file')
    nccreate([output_path,'N_xyz.nc'],'xt','Dimensions',{'dim1',nx});
    nccreate([output_path,'N_xyz.nc'],'yt','Dimensions',{'dim2',ny});
    nccreate([output_path,'N_xyz.nc'],'zt','Dimensions',{'dim3',nz});
    nccreate([output_path,'N_xyz.nc'],'time','Dimensions',{'dim4',nt});
end
ncwrite([output_path,'N_xyz.nc'],'xt',xt0);
ncwrite([output_path,'N_xyz.nc'],'yt',yt0);
ncwrite([output_path,'N_xyz.nc'],'zt',zt);
ncwrite([output_path,'N_xyz.nc'],'time',t);
fprintf(['Merge files completed! :) \n']);
end