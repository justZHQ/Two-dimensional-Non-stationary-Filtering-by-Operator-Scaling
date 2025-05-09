function [SEISMIC]=Nonstationary_function(seismic,dt,dx,filter_dt,length_h_t,length_h_x,fn_filter_H,fn_filter_L,v_filter_H,v_filter_L)
T=0:dt:(length(seismic(:,1))-1)*dt;
X=0:dx:(length(seismic(1,:))-1)*dx;

if filter_dt==0
    max_fn=max(max(fn_filter_H));
    filter_dt=1/(2*max_fn);
else 
    filter_dt=dt;
    max_fn=1/(2*filter_dt);
end


theory_dx=zeros(size(seismic));
for i=1:length(seismic(:,1))
    for j=1:length(seismic(1,:))
        theory_dx(i,j)=v_filter_L(i,j)/fn_filter_H(i,j);
    end
end
min_dx=min(min(theory_dx));
% max_dx=max(max(theory_dx));
% [min_dx_i,min_dx_j]=find(theory_dx==min_dx,1);
filter_dx=min_dx/2;
% v_filter=v_filter_L(min_dx_i,min_dx_j);
v_filter= floor(max_fn*min_dx);

chazhi_x=X(1):filter_dx:X(end);
chazhi_t=T(1):filter_dt:T(end);
chazhi_seismic=interp2(X,T,seismic,chazhi_x,chazhi_t','spline');
chazhi_v_filter_L=interp2(X,T,v_filter_L,chazhi_x,chazhi_t','nearest');
chazhi_v_filter_H=interp2(X,T,v_filter_H,chazhi_x,chazhi_t','nearest');
chazhi_fn_filter_H=interp2(X,T,fn_filter_H,chazhi_x,chazhi_t','nearest');
chazhi_fn_filter_L=interp2(X,T,fn_filter_L,chazhi_x,chazhi_t','nearest');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=max_fn;
G=1000*max_fn/v_filter;
N=round(length_h_t/2);
U=round(length_h_x/2);
reference_dt=filter_dt/K;
reference_dx=filter_dx/G;
reference_v=v_filter;
reference_fn=max_fn;
half_length_reference_n=N*K;
half_length_reference_u=U*G;
reference_n=0:1:half_length_reference_n-1;
reference_u=1:1:half_length_reference_u;
reference_t_filter=reference_dt*reference_n;
reference_x_filter=(2*reference_u-1)/2*reference_dx;

reference_h=zeros(length(reference_n),length(reference_u));
for i=1:length(reference_n)
    for j=1:length(reference_u)
        reference_h(i,j)=(1/(pi*reference_x_filter(j)))*((1/(2*pi*(reference_x_filter(j)/reference_v+reference_t_filter(i))+0.00000000000001)*(1-cos(2*pi*(reference_x_filter(j)/reference_v+reference_t_filter(i))*reference_fn)))+...
                                                         (1/(2*pi*(reference_x_filter(j)/reference_v-reference_t_filter(i))+0.00000000000001)*(1-cos(2*pi*(reference_x_filter(j)/reference_v-reference_t_filter(i))*reference_fn))));
    end
end
reference_h= reference_h*filter_dt*filter_dx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=zeros(2*N+1,2*U);
chazhi_SEISMIC=zeros(length(chazhi_seismic(:,1))+length(H(:,1))-1,length(chazhi_seismic(1,:))+length(H(1,:))-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 时间——空间域滤波矩阵形式
for i=1:length(chazhi_seismic(:,1))
    for j=1:length(chazhi_seismic(1,:))
%         a_FH_VL=(chazhi_fn_filter_H(i,j)/chazhi_v_filter_L(i,j))/(reference_fn/reference_v);
        a_FH_VL=chazhi_fn_filter_H(i,j)*reference_v*G/(chazhi_v_filter_L(i,j)*reference_fn);
        b_FH_VL=chazhi_fn_filter_H(i,j)*K/reference_fn;
        h_t=zeros(1,N+1);
        h_t(1)=1;
        for q=2:N+1
            h_t(q)=round((q-1)*b_FH_VL);
            if h_t(q)==0
                h_t(q)=1;
            end
        end
        h_x=zeros(1,U);
        h_x(1)=round(a_FH_VL/2);
        for p=2:U
            h_x(p)=round(p*a_FH_VL)-h_x(1);
            if h_x(p)==0
                h_x(p)=1;
            end
        end
        h_FH_VL=a_FH_VL/G*b_FH_VL/K*reference_h(h_t,h_x);
        h=h_FH_VL;
    
        if chazhi_v_filter_H(i,j)~=0
%             a_FH_VH=(chazhi_fn_filter_H(i,j)/chazhi_v_filter_H(i,j))/(reference_fn/reference_v);
            a_FH_VH=chazhi_fn_filter_H(i,j)*reference_v*G/(chazhi_v_filter_H(i,j)*reference_fn);
            b_FH_VH=chazhi_fn_filter_H(i,j)*K/reference_fn;
            h_t=zeros(1,N+1);
            h_t(1)=1;
            for q=2:N+1
                if h_t(q)==0
                    h_t(q)=1;
                end
                h_t(q)=round((q-1)*b_FH_VH);
            end
            h_x=zeros(1,U);
            h_x(1)=round(a_FH_VH/2);
            if h_x(1)==0
                h_x(1)=1;
            end
            for p=2:U
                h_x(p)=round(p*a_FH_VH)-h_x(1);
                if h_x(p)==0
                    h_x(p)=1;
                end
            end    
%             h_t
%             h_x
        h_FH_VH=a_FH_VH/G*b_FH_VH/K*reference_h(h_t,h_x);
        h=h_FH_VL-h_FH_VH;         
        end

        
        if chazhi_fn_filter_L(i,j)~=0
%             a_FL_VL=(chazhi_fn_filter_L(i,j)/chazhi_v_filter_L(i,j))/(reference_fn/reference_v);
            a_FL_VL=chazhi_fn_filter_L(i,j)*reference_v*G/(chazhi_v_filter_L(i,j)*reference_fn);
            b_FL_VL=chazhi_fn_filter_L(i,j)*K/reference_fn;
            h_t=zeros(1,N+1);
            h_t(1)=1;
            for q=2:N+1
                if h_t(q)==0
                    h_t(q)=1;
                end
                h_t(q)=round((q-1)*b_FL_VL);
            end
            h_x=zeros(1,U);
            h_x(1)=round(a_FL_VL/2);
            for p=2:U
                h_x(p)=round(p*a_FL_VL)-h_x(1);
                if h_x(p)==0
                    h_x(p)=1;
                end
            end              
        h_FL_VL=a_FL_VL/G*b_FL_VL/K*reference_h(h_t,h_x);
        h=h_FH_VL-h_FL_VL;             
        end
        
        if chazhi_fn_filter_L(i,j)~=0&&chazhi_v_filter_H(i,j)~=0
%             a_FL_VH=(chazhi_fn_filter_L(i,j)/chazhi_v_filter_H(i,j))/(reference_fn/reference_v);
            a_FL_VH=chazhi_fn_filter_L(i,j)*reference_v*G/(chazhi_v_filter_H(i,j)*reference_fn);
            b_FL_VH=chazhi_fn_filter_L(i,j)*K/reference_fn;
            h_t=zeros(1,N+1);
            h_t(1)=1;
            for q=2:N+1
                if h_t(q)==0
                    h_t(q)=1;
                end
                h_t(q)=round((q-1)*b_FL_VH);
            end
            h_x=zeros(1,U);
            h_x(1)=round(a_FL_VH/2);
            for p=2:U
                h_x(p)=round(p*a_FL_VH)-h_x(1);
                if h_x(p)==0
                    h_x(p)=1;
                end
            end              
        h_FL_VH=a_FL_VH/G*b_FL_VH/K*reference_h(h_t,h_x);        
        h=h_FH_VL-h_FH_VH-h_FL_VL+h_FL_VH;  
        end
        
        h=[fliplr(h),h];
        h=[flipud(h);h(2:end,:)];
        chazhi_SEISMIC(i:i+length(H(:,1))-1,j:j+length(H(1,:))-1)=chazhi_SEISMIC(i:i+length(H(:,1))-1,j:j+length(H(1,:))-1)+h*chazhi_seismic(i,j);
        clear h;
    end
end
chazhi_SEISMIC_mid=chazhi_SEISMIC((length(H(:,1))-1)/2+1:length(chazhi_SEISMIC(:,1))-(length(H(:,1))-1)/2,(length(H(1,:))-1)/2+1:length(chazhi_SEISMIC(1,:))-(length(H(1,:))-1)/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear reference_h;
SEISMIC=interp2(chazhi_x,chazhi_t,chazhi_SEISMIC_mid,X(1:end-1),T(1:end-1)','spline');
x_last=interp1(chazhi_t,chazhi_SEISMIC_mid(:,end),T(1:end-1)','spline');
t_last=interp1(chazhi_x,chazhi_SEISMIC_mid(end,:),X(1:end-1));
SEISMIC=[SEISMIC,x_last;[t_last,chazhi_SEISMIC_mid(end)]];