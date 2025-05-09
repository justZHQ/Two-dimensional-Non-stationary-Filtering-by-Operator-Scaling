clear ;
close all;
clc;

dt=0.001;
dx=30;
h=0:dx:2010;
fs=1/dt;
f0=50;
tmax=1.1;
tau =[0.2,0.42,0.555];
v=[2800,2850,2900];
k=linspace(-1,1,length(h));
amp=[1,1,1];
snr=1.5;
L=4;
[signal,h,t] = hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L);

d_085=zeros(size(signal));
for i=1:length(h)
  tau_085=0.87;
  v_085=3200;
  amp_085=k(i);
  [d_mid,h,t] = hyperbolic_events(dt,f0,tmax,h,tau_085,v_085,amp_085,snr,L);
  d_085(:,i)=d_mid(:,i);
end
signal=signal+d_085;

%% 动校正
max_stretch=1000;
tnmo=[0.2,0.42,0.555,tau_085];
vnmo=[2800,2830,2900,3200];
[signal_nmo,M,ti,vi] = nmo(signal,dt,h,tnmo,vnmo,max_stretch);
signal_nmo(1:260,30:end)=0;

%% 添加干扰
tau_interfere=[0.35,0.5,0.8];
v_interfere=[0.00000003,0.00000004,0.00000002];
amp_interfere=[1,1,1];snr=1000000;
[interfere,h,t,Apparent_velocity] = parabola_events_ApparentVelocityMatrix(dt,f0,tmax,h,tau_interfere,v_interfere,amp_interfere,snr,L);

%%
X=h;T=t;
data=signal_nmo+interfere;
scale=1;
figure;   
wigb(signal,scale,X,T); 
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
txt1 = 'p1';txt2 = 'p2';txt3 = 'p3';txt4 = 'p4';
text(-190,0.12,txt1,'FontSize',16,'color','r');
text(-190,0.32,txt2,'FontSize',16,'color','r');
text(-190,0.52,txt3,'FontSize',16,'color','r');
text(-190,0.88,txt4,'FontSize',16,'color','r');
figure;   
wigb(signal_nmo,scale,X,T); 
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
txt1 = 'p1';txt2 = 'p2';txt3 = 'p3';txt4 = 'p4';
text(-190,0.12,txt1,'FontSize',16,'color','r');
text(-190,0.32,txt2,'FontSize',16,'color','r');
text(-190,0.52,txt3,'FontSize',16,'color','r');
text(-190,0.88,txt4,'FontSize',16,'color','r');

figure;   
wigb(data,scale,X,T); 
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
txt1 = 'm1';txt2 = 'm2';txt3 = 'm3';
text(-210,0.31,txt1,'FontSize',16,'color','r');
text(-210,0.48,txt2,'FontSize',16,'color','r');
text(-210,0.72,txt3,'FontSize',16,'color','r');


length_h_t=100;
length_h_x=60;
v_reference=10000;
SCALE=1000;
fn_filter_L=0;
fn_filter_H=150;
Seismic=data;

width=30;
v=ones(size(Seismic))*1000000;
for i=1:length(tnmo)
v(tnmo(i)/dt-width:tnmo(i)/dt+width,:)=35000;
end

figure;
imagesc(X,T,v);
colormap(gray);
ch=colorbar('YTickLabel',{(0:10000:40000),'inf'});
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
set(get(ch,'title'),'position',[-10 250],'string','Apparent velocity[m/s]');
caxis([0,100000]);


[H_reference]=Matrix_H_reference(dt,dx,length_h_t,length_h_x,fn_filter_L,fn_filter_H,v_reference,SCALE);
[interfere_filter]=fan_function_calculation(Seismic,H_reference,length_h_t,length_h_x,v_reference,SCALE,v);

signal_filter=Seismic-interfere_filter;
figure;   
wigb(interfere_filter,scale,X,T); 
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);

figure;   
wigb(signal_filter,scale,X,T); 
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);

%% 平稳滤波
v=ones(size(Seismic))*35000;

[SEISMIC_stationary]=fan_function_calculation(Seismic,H_reference,length_h_t,length_h_x,v_reference,SCALE,v);
figure;   
wigb(SEISMIC_stationary,scale,X,T); 
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
figure;   
wigb(Seismic-SEISMIC_stationary,scale,X,T); 
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);