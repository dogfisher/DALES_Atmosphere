function atmos_interf_ts(input_path,orb_slc)
%Generate atmospheric interferograms using slant delay
%Written by Fengming in TJU
if nargin<1
    clear,clc
    fclose all;
    close all;
    input_path='D:/Study/TUD_Work/data2/N_delay_out/Slant_delay_out/';
end
nx=length(orb_slc.xt);
ny=length(orb_slc.yt);
nt=length(orb_slc.t);
%=========
atmos_x_sp=orb_slc.xt(2)-orb_slc.xt(1);
atmos_y_sp=orb_slc.yt(2)-orb_slc.yt(1);
x_max=max(orb_slc.xt);
y_max=max(orb_slc.yt);
x_min=min(orb_slc.xt);
y_min=min(orb_slc.yt);
area_pt=[[x_min x_min x_max x_max x_min]',[y_min y_max y_max y_min y_min]'];
[ifg_grid_img_x,ifg_grid_img_y]=meshgrid(orb_slc.ifg_x_sp*[0.5:orb_slc.ifg_azi-0.5]...
    ,orb_slc.ifg_y_sp*[0.5:orb_slc.ifg_rg-0.5]);
slc_grid=tran_forward([ifg_grid_img_x(:),ifg_grid_img_y(:)],orb_slc.ifg_geo_par);
slc_grid_x=reshape(slc_grid(:,1),[orb_slc.ifg_rg,orb_slc.ifg_azi]);
slc_grid_y=reshape(slc_grid(:,2),[orb_slc.ifg_rg,orb_slc.ifg_azi]);

%%================
%% generate total delay time series, output files are *.fig *.tif and a video
% for dataset 1, 
% total_clim=[2260 2340];
% for dataset 2
%total_clim=[1220 1243];%for azi 30
total_clim=[1225 1250];%for azi 120
%total_clim=[1225 1250];%for azi 210
% total_clim=[1225 1251];%for azi 300

% % for dataset 3
% % total_clim=[1110 1155];
% total_clim=[];
for i=1:nt
    D_delay1=(ncread([input_path,'D_delay.nc'],'D_delay',[1,1,i],[inf,inf,1]));
%     imagesc(D_delay1);
    disp_atmos_ifg(D_delay1*1000,area_pt,orb_slc.ifg_pt,slc_grid_x,slc_grid_y,'Total Delay (mm)',total_clim);
    title(['Epoch',num2str(i),'  ',num2str(orb_slc.t(i)),'s']);
    saveas(gcf,strcat(input_path,'Total_delay_epoch',num2str(i),'.tif'));
    saveas(gcf,strcat(input_path,'Total_delay_epoch',num2str(i),'.fig'));
    close(gcf);
end
%p
myVideo=VideoWriter([input_path,'Total_delay.avi']);
myVideo.FrameRate = 1.5;
myVideo.Quality=100;
open(myVideo);
for i=1:nt
    a=imread(strcat(input_path,'Total_delay_epoch',num2str(i),'.tif'));
    writeVideo(myVideo,a);
end
close(myVideo);

%%
%%generate atmospheric interferograms
for type=0:1
    if type==0
        output_path=[input_path,'TS_phase_single'];
        ph_clim=[-pi pi];
        video_name=[output_path,'/ph_single.avi'];
    else
        output_path=[input_path,'TS_phase_daisy_chain'];
        ph_clim=[-0.5 0.5];
        video_name=[output_path,'/ph_daisy.avi'];
    end
    if exist(output_path,'dir')
        delete([output_path,'/*']);
    else
        mkdir(output_path);
    end
    output_path=[output_path,'/'];

    nccreate([output_path,'D_delay_ph.nc'],'D_delay' ...
        ,'Datatype','single','Dimensions',{'dim1',orb_slc.ifg_rg,'dim2',orb_slc.ifg_azi,'dim3',nt});
    
    for i=1:nt-1
        if type==0
            D_delay0=(ncread([input_path,'D_delay.nc'],'D_delay',[1,1,1],[inf,inf,1]));
        else
            D_delay0=(ncread([input_path,'D_delay.nc'],'D_delay',[1,1,i],[inf,inf,1]));
        end
        D_delay1=(ncread([input_path,'D_delay.nc'],'D_delay',[1,1,i+1],[inf,inf,1]));
        
        d_ph=wrap((D_delay1-D_delay0)*4*pi/orb_slc.lam);
        disp_atmos_ifg(d_ph,area_pt,orb_slc.ifg_pt,...
            slc_grid_x,slc_grid_y,'Phase (rad)',ph_clim);
        
        title(['Epoch',num2str(i),'  ',num2str(orb_slc.t(i)),'s']);
        saveas(gcf,strcat(output_path,'/Phase_epoch',num2str(i),'.tif'));
        saveas(gcf,strcat(output_path,'/Phase_epoch',num2str(i),'.fig'));
        close(gcf);
        ncwrite([output_path,'/D_delay_ph.nc'],'D_delay',d_ph,[1,1,i]);
    end
    
    myVideo=VideoWriter(video_name);
    myVideo.FrameRate = 1.5;
    myVideo.Quality=100;
    open(myVideo);
    for i=1:nt-1
        a=imread(strcat(output_path,'Phase_epoch',num2str(i),'.tif'));
        writeVideo(myVideo,a);
    end
    close(myVideo);
end

end