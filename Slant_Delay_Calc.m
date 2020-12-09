function Slant_Delay_Calc(input_path,orb_slc,is_parallel)
%Calculate the slant delay for all epoch
%Wrtten by Fengming in TJU
if exist([input_path,'D_delay.nc'],'file')
    delete([input_path,'D_delay.nc']);
end

n_patch=orb_slc.ifg_azi;
if is_parallel
    ncores=feature('numCores');
    d_patch=round(n_patch/ncores);
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        parpool('local',ncores);
    end
end

nt=length(orb_slc.t);
nz=length(orb_slc.zt);
nccreate([input_path,'D_delay.nc'],'D_delay' ...
        ,'Datatype','single','Dimensions',{'dim1',orb_slc.ifg_rg,...
        'dim2',orb_slc.ifg_azi,'dim3',nt});
    
[xx,zz]=meshgrid(orb_slc.ifg_x_sp*[0.5:1:orb_slc.model_rg-0.5],orb_slc.zt);
rg_offset=distance_2D(orb_slc.ifg_pt(1,:),orb_slc.model_pt(1,:))/orb_slc.ifg_x_sp;

drg=orb_slc.z_max*tan(orb_slc.inc);

rg_list=zeros(orb_slc.ifg_rg,nz);
for i=1:orb_slc.ifg_rg
    rg_be=rg_offset*orb_slc.ifg_x_sp+(i-0.5)*orb_slc.ifg_x_sp;
    rg_ed=rg_be-drg(i);
    ifg_rg_list=linspace(rg_be,rg_ed,nz);
    rg_list(i,:)=ifg_rg_list;
end

for i=1:nt
    fprintf(['Slant delay calculation for epoch ',num2str(i),'\n']);
    N_slc=ncread([input_path,'D_delay_rsp',num2str(i,'%.3d'),'.nc'],'N_delay');
    D_slant_delay=zeros(orb_slc.ifg_rg,orb_slc.ifg_azi);
    if is_parallel
        spmd(ncores)
            for k=1:ncores
                ibe=(k-1)*d_patch+1;
                if k==ncores
                    ied=n_patch;
                else
                    ied=k*d_patch;
                end             
                if labindex==k
                    i_idx=ibe:ied;
                end
            end
            D_slant_delay0=slice2delay(xx,zz,N_slc(:,i_idx,:),rg_list,orb_slc.zt,orb_slc.dh,orb_slc.inc);  
        end
        for k=1:ncores
            ibe=(k-1)*d_patch+1;
            if k==ncores
                ied=n_patch;
            else
                ied=k*d_patch;
            end
            D1=D_slant_delay0{k};
            D_slant_delay(:,ibe:ied)=D1(:,1:ied-ibe+1);
        end
    else
        D_slant_delay=slice2delay(xx,zz,N_slc,rg_list,orb_slc.zt,orb_slc.dh,orb_slc.inc);
    end
%     imagesc(D_slant_delay);
%     colorbar
%     axis image
    ncwrite([input_path,'D_delay.nc'],'D_delay',D_slant_delay,[1,1,i]);
end
fprintf(['Slant delay calculation completed! :) \n']);
end