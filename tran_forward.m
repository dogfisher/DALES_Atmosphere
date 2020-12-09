function pt1=tran_forward(pt,par)
np=size(pt,1);
B=zeros(np*2,6);
tmp=[ones(np,1) pt(:,1:2)];
B(1:2:end,1:3)=tmp;
B(2:2:end,4:6)=tmp;
y=B*par;
pt1=[y(1:2:end) y(2:2:end)];
end