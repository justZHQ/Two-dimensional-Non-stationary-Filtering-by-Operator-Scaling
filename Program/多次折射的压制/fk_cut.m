function spec_cut=fk_cut(spec,f,kx,leftpoint_low,rightpoint_low,leftpoint_heigh,rightpoint_heigh)
[~,nkx]=size(spec);
vector_low=zeros(1,nkx);
vector_heigh=zeros(1,nkx);
df=f(2)-f(1);
for i=1:nkx/2 
    vector_low(i)=leftpoint_low(2)/leftpoint_low(1).*kx(i);
    vector_low(i+nkx/2)=rightpoint_low(2)/rightpoint_low(1).*kx(nkx/2+i);
    vector_heigh(i)=leftpoint_heigh(2)/leftpoint_heigh(1).*kx(i);
    vector_heigh(i+nkx/2)=rightpoint_heigh(2)/rightpoint_heigh(1).*kx(nkx/2+i);
end
spec_cut=spec;
for i=1:nkx
    f_edge_low=floor(vector_low(i)./df);
    f_edge_heigh=floor(vector_heigh(i)./df);
    if f_edge_low==0
        f_edge_low=1;
    end
    if f_edge_heigh==0
        f_edge_heigh=1;
    end
    spec_cut(f_edge_low:f_edge_heigh,i)=0;
end
    
end