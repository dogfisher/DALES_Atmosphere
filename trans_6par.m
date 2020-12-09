function par=trans_6par(pt0,pt1)
np=size(pt0,1);
B=zeros(np*2,6);
l=reshape(pt1',[],1);
tmp=[ones(np,1) pt0];
B(1:2:end,1:3)=tmp;
B(2:2:end,4:6)=tmp;
par=(B'*B)\(B'*l);
end