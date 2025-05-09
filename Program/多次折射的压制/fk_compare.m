clear ;
close all;
clc;
%% 
dt=0.001;nt=1000;
dx=10;nx=200;
T=0:dt:(nt-1)*dt;
X=-nx/2*dx:dx:(nx/2-1)*dx;

L=30;%子波半长
fmax=145;%最大显示频率/Hz
fa=0;%相位
f0=40;
%% Primary wave
amat=zeros(nt,nx);
tnot_primary=[0.07,0.35,0.5,0.65,0.85];
xnot_primary=X(nx/2);
v_primary=5000;
amp_primary=[1,1,1,1,1];
flag_primary=[1,1,1,1,1];
for i=1:length(tnot_primary)
    [amat]=event_hyp(amat,T,X,tnot_primary(i),xnot_primary,v_primary,amp_primary(i),flag_primary(i));
end
%% Interference wave
x_end_right=600;
x_end_left=-600;
x_point_right=[170,260,340,450];
t_start_right=Initial_point(T,X,x_end_right,x_point_right);
x_point_left=-x_point_right;
t_start_left=Initial_point(T,X,x_end_left,x_point_left);
Slope_parameter=0.2;
for i=1:length(x_point_right)
    amat=event_dip(amat,T,X,[t_start_right(i),t_start_right(i)+Slope_parameter],[x_point_right(i),x_point_right(i)+1000]);
    amat=event_dip(amat,T,X,[t_start_left(i),t_start_left(i)+Slope_parameter],[x_point_left(i),x_point_left(i)-1000]);
end
%% Convolution with Ricker wavelet
[t_wavelet,wavelet,f_wavelet,amplitude_spectrum_wavelet]=Ricker_my(dt,L,f0,fa,fmax);
waveletname='Ricker';
seismic=zeros(size(amat));
for i=1:nx
    seismic_original=conv(amat(:,i),wavelet);
    seismic(:,i)=seismic_original(L+1:end-L);
end
%% f-k spectrum
[spec_seismic,f,kx]=fktran(seismic,T,X);
leftpoint_low=[kx(1),f(220)];
rightpoint_low=[kx(end),f(220)];
leftpoint_heigh=[kx(1),f(330)];
rightpoint_heigh=[kx(end),f(330)];
spec_signal=fk_cut(spec_seismic,f,kx,leftpoint_low,rightpoint_low,leftpoint_heigh,rightpoint_heigh);
[seis,t,x]=ifktran(spec_signal,f,kx);
signal=seis(1:nt,1:nx);
Interference=seismic-signal;
[spec_Interference,~,~]=fktran(Interference,T,X);
mag_spec_seismic=abs(spec_seismic);
mag_spec_Interference=abs(spec_Interference);
mag_spec_signal=abs(spec_signal);
%%
figure;
imagesc(kx,f(1:200),mag_spec_seismic(1:200,:));
set(gca,'YDir','normal')
% title('Synthetic seismic gather');
xlabel('Wavenumber(rad/m)','FontName','Times New Roman');
ylabel('Frequency(HZ)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
colormap(jet);
colorbar;
caxis([0,8]);

figure;
imagesc(kx,f(1:200),mag_spec_Interference(1:200,:));
set(gca,'YDir','normal')
% title('Effective wave');
xlabel('Wavenumber(rad/m)','FontName','Times New Roman');
ylabel('Frequency(HZ)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
colormap(jet);
colorbar;
caxis([0,8]);

figure;
imagesc(kx,f(1:200),mag_spec_signal(1:200,:));
set(gca,'YDir','normal')
% title('Iterference wave');
xlabel('Wavenumber(rad/m)','FontName','Times New Roman');
ylabel('Frequency(HZ)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
colormap(jet);
colorbar;
caxis([0,8]);

reflection_cut=zeros(nt,nx/2);
seismic_cut=zeros(nt,nx/2);
signal_cut=zeros(nt,nx/2);
Interference_cut=zeros(nt,nx/2);
X_cut=zeros(1,nx/2);

for i=1:nx/2
    reflection_cut(:,i)=amat(:,2*i-1);
    seismic_cut(:,i)=seismic(:,2*i-1);
    signal_cut(:,i)=signal(:,2*i-1);
    Interference_cut(:,i)=Interference(:,2*i-1);
    X_cut(i)=X(2*i-1);
end


figure;
wigb(seismic_cut,1,X_cut,T);
% title('Synthetic seismic gather');
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
figure;
wigb(signal_cut,1,X_cut,T);
% title('Effective wave');
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
figure;
wigb(Interference_cut,1,X_cut,T);
% title('Iterference wave');
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);