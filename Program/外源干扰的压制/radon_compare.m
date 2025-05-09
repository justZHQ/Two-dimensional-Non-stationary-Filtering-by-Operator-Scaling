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
tnot_primary=[0.2,0.3,0.5,0.6,0.8];
xnot_primary=X(nx/2);
v_primary=5000;
amp_primary=[1,1,1,1,1];
flag_primary=[1,1,1,1,1];
aper=[inf,inf,inf,inf,500];
m=[X(end)-X(1),X(end)-X(1),0,0,0];
w=[0,0,0,0,10];
kw=[0,0,1500,-500,0];
for i=1:length(tnot_primary)
    [amat]=event_loosehyp(amat,T,X,tnot_primary(i),xnot_primary,v_primary,amp_primary(i),flag_primary(i),aper(i),m(i),w(i),kw(i));
end
%% Interference wave
tnot_interference=[0.15,0.4,0.2,0.75];
xnot_interference=[X(20),X(95),X(150),X(160)];
v_interference=[2300,3000,2200,2000];
amp_interference=[1,1,1,1];
flag_interference=[0,0,0,0];
for i=1:length(tnot_interference)
    [amat]=event_hyp(amat,T,X,tnot_interference(i),xnot_interference(i),v_interference(i),amp_interference(i),flag_interference(i));
end
%% Convolution with Ricker wavelet
[t_wavelet,wavelet,f_wavelet,amplitude_spectrum_wavelet]=Ricker_my(dt,L,f0,fa,fmax);
waveletname='Ricker';
Seismic=zeros(size(amat));
for i=1:nx
    seismic_original=conv(amat(:,i),wavelet);
    Seismic(:,i)=seismic_original(L+1:end-L);
end
%% Radon transformation  qmin=min(radon_p);
radon_dt=dt;
radon_h=X;
radon_flow=2;
radon_fhigh=150;
radon_N=2;
dp=0.02;
radon_p=-80*dp:dp:dp*80;
radon_mu=0.1;
d=Seismic;
[R_seismic]=inverse_radon_freq(d,radon_dt,radon_h,radon_p,radon_N,radon_flow,radon_fhigh,radon_mu,'ls');
qmax=max(radon_p);
qmin=min(radon_p);
nq=length(radon_p);
% q_cutstart=[radon_p(1),radon_p(2);radon_p(1),radon_p(2);radon_p(1),radon_p(2)];
% q_cutend=[radon_p(end),radon_p(end-1);radon_p(end),radon_p(end-1);radon_p(end),radon_p(end-1)];
q_cutstart=[0.75,0;-0.14,0.3;-0.14,1];
q_cutend=[radon_p(end),0;0.18,0.4;0.5,1];
% q_cutstart=[radon_p(2),0;radon_p(2),1];
% q_cutend=[radon_p(end),0;radon_p(end),1];
[signal,m,tau,q]=pradon_demultiple_my(Seismic,dt,radon_h,qmin,qmax,nq,radon_flow,radon_fhigh,radon_mu,q_cutstart,q_cutend);
Interference=Seismic-signal;
[R_interference]=inverse_radon_freq(Interference,radon_dt,radon_h,radon_p,radon_N,radon_flow,radon_fhigh,radon_mu,'ls');
R_signal=R_seismic-R_interference;
%% Sorting gather
reflection_cut=zeros(nt,nx/2);
seismic_cut=zeros(nt,nx/2);
signal_cut=zeros(nt,nx/2);
Interference_cut=zeros(nt,nx/2);
X_cut=zeros(1,nx/2);
for i=1:nx/2
    seismic_cut(:,i)=Seismic(:,2*i-1);
    signal_cut(:,i)=signal(:,2*i-1);
    Interference_cut(:,i)=Interference(:,2*i-1);
    X_cut(i)=X(2*i-1);
end
%%
figure;
wigb(seismic_cut,1,X_cut,T);
title('Synthetic seismic gather');
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
figure;
wigb(signal_cut,1,X_cut,T);
title('Effective wave');
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
figure;
wigb(Interference_cut,1,X_cut,T);
title('Iterference wave');
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);

% figure;
% imagesc(signal_cut);
% colorbar;
% figure;
% imagesc(Interference_cut);
% colorbar;

figure;
imagesc(radon_p,T,R_seismic);
title('Synthetic seismic gather');
colormap(seismic(3));
xlabel('q (Residual Moveout) [s]','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
colorbar;
caxis([-0.1,0.2])

figure;
imagesc(radon_p,T,R_interference);
title('Iterference wave');
colormap(seismic(3));
xlabel('q (Residual Moveout) [s]','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
colorbar;
caxis([-0.1,0.2])

figure;
imagesc(radon_p,T,R_signal);
title('Effective wave');
colormap(seismic(3));
xlabel('q (Residual Moveout) [s]','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
colorbar;
caxis([-0.1,0.2])






