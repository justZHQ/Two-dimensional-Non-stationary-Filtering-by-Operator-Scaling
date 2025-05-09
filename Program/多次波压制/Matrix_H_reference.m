function [H_reference]=Matrix_H_reference(dt,dx,length_h_t,length_h_x,fn_filter_L,fn_filter_H,v_reference,SCALE)
filter_dt=dt;
filter_dx=dx/SCALE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=round(length_h_t/2);
U=round(length_h_x/2);
n=-N:1:N;
u=1:1:U*SCALE;
t_filter=filter_dt*n;
x_filter=(2*u-1)/2*filter_dx; 
% H_reference=zeros(length(n),length(u));
H_reference_fh=zeros(length(n),length(u));
H_reference_fl=zeros(length(n),length(u));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for p=1:length(n)
    for q=1:length(u)
        H_reference_fh(p,q)=(1/(pi*x_filter(q)))*((1/(2*pi*(x_filter(q)/(v_reference+0.00000000000001)+t_filter(p))+0.00000000000001)*(1-cos(2*pi*(x_filter(q)/(v_reference+0.00000000000001)+t_filter(p))*fn_filter_H)))+...
            (1/(2*pi*(x_filter(q)/(v_reference+0.00000000000001)-t_filter(p))+0.00000000000001)*(1-cos(2*pi*(x_filter(q)/(v_reference+0.00000000000001)-t_filter(p))*fn_filter_H))));
        H_reference_fl(p,q)=(1/(pi*x_filter(q)))*((1/(2*pi*(x_filter(q)/(v_reference+0.00000000000001)+t_filter(p))+0.00000000000001)*(1-cos(2*pi*(x_filter(q)/(v_reference+0.00000000000001)+t_filter(p))*fn_filter_L)))+...
            (1/(2*pi*(x_filter(q)/(v_reference+0.00000000000001)-t_filter(p))+0.00000000000001)*(1-cos(2*pi*(x_filter(q)/(v_reference+0.00000000000001)-t_filter(p))*fn_filter_L))));
    end
end
H_reference=H_reference_fh-H_reference_fl;
H_reference= H_reference*filter_dt*dx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

