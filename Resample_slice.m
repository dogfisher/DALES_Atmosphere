function Resample_slice(input_path,orb_slc)
if nargin<1
    clear,clc
    fclose all;
    close all;
    input_path='D:/Study/TUD_Work/data2/N_delay_out/';
end
output_path=[input_path,'Slant_delay_out'];
if exist(output_path,'dir')
    delete([output_path,'/*']);
else
    mkdir(output_path);
end

output_path=[output_path,'/'];

[model_grid_img_x,model_grid_img_y]=...
    meshgrid(orb_slc.ifg_x_sp*[0.5:orb_slc.model_azi-0.5],...
    orb_slc.ifg_y_sp*[0.5:orb_slc.model_rg-0.5]);
slc_grid=tran_forward([model_grid_img_x(:),model_grid_img_y(:)],orb_slc.geo_par);
%ifg grid
slc_grid_x=reshape(slc_grid(:,1),[orb_slc.model_rg,orb_slc.model_azi]);
slc_grid_y=reshape(slc_grid(:,2),[orb_slc.model_rg,orb_slc.model_azi]);

%Orignal model grid
[atmos_grid_x,atmos_grid_y]=meshgrid(orb_slc.xt,orb_slc.yt);

atmos_file=dir([input_path,'N_epoch*.nc']);    
nt=length(atmos_file);
nz=length(orb_slc.zt);
for i=1:nt
    fprintf(['Resample epoch ',num2str(i),'\n']);
    N_model=ncread([input_path,atmos_file(i).name],'N_delay');
    N_slc=zeros(orb_slc.model_rg,orb_slc.model_azi,nz);%Range, Azimuth and height
    for j=1:nz
        Vq = interp2(atmos_grid_x,atmos_grid_y,N_model(:,:,j),slc_grid_x,slc_grid_y);
        [tmp_i,tmp_j]=find(isnan(Vq));
        if ~isempty(tmp_i)
            n_rsp=length(tmp_i);
            slc_rsp_x=zeros(n_rsp,1);
            slc_rsp_y=zeros(n_rsp,1);
            for k=1:n_rsp
                ii=tmp_i(k);
                jj=tmp_j(k);
                slc_rsp_x(k)=slc_grid_x(ii,jj);
                slc_rsp_y(k)=slc_grid_y(ii,jj);
            end
            Vq1=interp2(atmos_grid_x,atmos_grid_y,N_model(:,:,j),...
                slc_rsp_x,slc_rsp_y,'makima');
            
            for k=1:n_rsp
                ii=tmp_i(k);
                jj=tmp_j(k);
                Vq(ii,jj)=Vq1(k);
            end
        end
        N_slc(:,:,j)=Vq;
%         imagesc(N_slc(:,:,j));
%         axis image
%         colorbar
    end
    if ~exist([output_path,'D_delay_rsp',num2str(i,'%.3d'),'.nc'],'file')
        nccreate([output_path,'D_delay_rsp',num2str(i,'%.3d'),'.nc'],'N_delay' ...
            ,'Datatype','single','Dimensions',{'dim1',orb_slc.model_rg,'dim2',orb_slc.model_azi,'dim3',nz});
    end
    ncwrite([output_path,'D_delay_rsp',num2str(i,'%.3d'),'.nc'],'N_delay',N_slc);
end
fprintf(['Resample completed! :) \n']);
end