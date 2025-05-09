function [h]=Matrix_h(H_reference,length_h_x,v_reference,SCALE,v)

u=round(length_h_x/2);
a=round(SCALE*v_reference/v);
hx=zeros(1,u);
hx(1)=round(a/2);

for i=1:u-1
    hx(i+1)=i*a+hx(1);
end
for i=1:u
    if hx(i)==0
            hx(i)=1;
    end
end
% disp(hx);
h=H_reference(:,hx)*a/SCALE;
h=[fliplr(h),h];   
    
    
    
    
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%