function [SEISMIC,Time_difference]=fan_function_calculation(seismic,H_reference,length_h_t,length_h_x,v_reference,SCALE,v)
tic  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=2*round(length_h_t/2)+1;
u=2*round(length_h_x/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 时间——空间域滤波矩阵形式
SEISMIC=zeros(length(seismic(:,1))+n-1,length(seismic(1,:))+u-1);
for i=1:length(seismic(:,1))
    for j=1:length(seismic(1,:))
        [h]=Matrix_h(H_reference,length_h_x,v_reference,SCALE,v(i,j));
        SEISMIC(i:i+n-1,j:j+u-1)=SEISMIC(i:i+n-1,j:j+u-1)+h*seismic(i,j);
%         [h]=Matrix_h(H_reference,length_h_x,v_reference,SCALE,v(i,j));
%         [g]=Matrix_g(G_reference,length_h_x,v_reference,SCALE,v(i,j));
%         SEISMIC(i:i+n-1,j:j+u-1)=SEISMIC(i:i+n-1,j:j+u-1)+(g-h)*seismic(i,j);
        clear h;
    end
end
SEISMIC=SEISMIC((n-1)/2+1:length(SEISMIC(:,1))-(n-1)/2,(u)/2:length(SEISMIC(1,:))-(u)/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time_difference=toc;