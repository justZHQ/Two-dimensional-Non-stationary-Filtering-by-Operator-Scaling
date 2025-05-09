function [G_reference]=Matrix_G_reference(dt,dx,length_h_t,length_h_x,fn_filter_H,v_reference,SCALE)
filter_dt=dt;
filter_dx=dx/SCALE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=round(length_h_t/2);
U=round(length_h_x/2);
n=-N:1:N;
u=1:1:U*SCALE;
t_filter=filter_dt*n;
x_filter=(2*u-1)/2*filter_dx; 
G_reference=zeros(length(n),length(u));

for i=1:length(n)
    for j=1:length(u)           
        if i~=1                                             
            G_reference(i,j)= sin(2*pi*x_filter(j)*fn_filter_H/v_reference)*sin(2*pi*t_filter(i)*fn_filter_H)/(pi^2*t_filter(i)*x_filter(j)+0.000000000001);     
        else
            G_reference(i,j)=2*fn_filter_H*sin(2*pi*x_filter(j)*fn_filter_H/v_reference)/(pi*x_filter(j));
        end                                  
    end
end
G_reference= G_reference*filter_dt*dx;