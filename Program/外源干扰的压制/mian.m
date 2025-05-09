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
v_interference=[2300,2500,2200,2000];
amp_interference=[1,1,1,1];
flag_interference=[2,2,2,2];
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
%% Two-dimensional nonstationary filtering (Pie slice filter)
filter_dt=0;
length_h_t=120;
length_h_x=80;
signal=Seismic;
for ii=1:length(tnot_interference)
    Apparent_velocity=zeros(nt,nx);
    [Apparent_velocity]=ApparentVelocityMatrix(T,X,tnot_interference(ii),xnot_interference(ii),v_interference(ii),floor(L),Apparent_velocity);
    fn_filter_H=zeros(size(Seismic));
    fn_filter_L=zeros(size(Seismic));
    v_filter_H=ones(size(Seismic))*10000;
    v_filter_L=ones(size(Seismic))*8000;
    MAX_Apparent_velocity=3500;
    for i=1:length(Seismic(:,1))
        for j=1:length(Seismic(1,:))
            if Apparent_velocity(i,j)==0
                fn_filter_H(i,j)=150;
                fn_filter_L(i,j)=130;
            else
                fn_filter_H(i,j)=150;
                fn_filter_L(i,j)=0;
                if Apparent_velocity(i,j)<=MAX_Apparent_velocity
                    v_filter_H(i,j)=Apparent_velocity(i,j)*1.5;
                    v_filter_L(i,j)=Apparent_velocity(i,j)*0.6;
                elseif Apparent_velocity(i,j)>MAX_Apparent_velocity
                    v_filter_H(i,j)=MAX_Apparent_velocity*0;
                    v_filter_L(i,j)=MAX_Apparent_velocity*0.5;
                end
            end
        end
    end
    [noise]=Nonstationary_function(signal,dt,dx,filter_dt,length_h_t,length_h_x,fn_filter_H,fn_filter_L,v_filter_H,v_filter_L);
    signal=signal-noise;
end
Interference=Seismic-signal;
%% Apparent_velocity of interference wave
Apparent_velocity=zeros(nt,nx);
for ii=1:length(tnot_interference)
    [Apparent_velocity]=ApparentVelocityMatrix(T,X,tnot_interference(ii),xnot_interference(ii),v_interference(ii),floor(L),Apparent_velocity);
end
%% f-k spectrum
[spec_seismic,f,kx]=fktran(Seismic,T,X);
mag_spec_seismic=abs(spec_seismic);
[spec_Interference,~,~]=fktran(Interference,T,X);
mag_spec_Interference=abs(spec_Interference);
[spec_signal,~,~]=fktran(signal,T,X);
mag_spec_signal=abs(spec_signal);
%% Sorting gather
reflection_cut=zeros(nt,nx/2);
seismic_cut=zeros(nt,nx/2);
signal_cut=zeros(nt,nx/2);
Interference_cut=zeros(nt,nx/2);
X_cut=zeros(1,nx/2);

for i=1:nx/2
    reflection_cut(:,i)=amat(:,2*i-1);
    seismic_cut(:,i)=Seismic(:,2*i-1);
    signal_cut(:,i)=signal(:,2*i-1);
    Interference_cut(:,i)=Interference(:,2*i-1);
    X_cut(i)=X(2*i-1);
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
[R_seismic]=inverse_radon_freq(Seismic,radon_dt,radon_h,radon_p,radon_N,radon_flow,radon_fhigh,radon_mu,'ls');
[R_signal]=inverse_radon_freq(signal,radon_dt,radon_h,radon_p,radon_N,radon_flow,radon_fhigh,radon_mu,'ls');
[R_interference]=inverse_radon_freq(Interference,radon_dt,radon_h,radon_p,radon_N,radon_flow,radon_fhigh,radon_mu,'ls');
%%
figure;
imagesc(kx,f(1:200),mag_spec_seismic(1:200,:));
set(gca,'YDir','normal')
title('Synthetic seismic gather');
xlabel('Wavenumber(rad/m)','FontName','Times New Roman');
ylabel('Frequency(HZ)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
colormap(jet);
colorbar;
caxis([0,5])

figure;
imagesc(kx,f(1:200),mag_spec_Interference(1:200,:));
set(gca,'YDir','normal')
title('Effective wave');
xlabel('Wavenumber(rad/m)','FontName','Times New Roman');
ylabel('Frequency(HZ)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
colormap(jet);
colorbar;
caxis([0,5])

figure;
imagesc(kx,f(1:200),mag_spec_signal(1:200,:));
set(gca,'YDir','normal')
title('Iterference wave');
xlabel('Wavenumber(rad/m)','FontName','Times New Roman');
ylabel('Frequency(HZ)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
colormap(jet);
colorbar;
caxis([0,5])

figure;
plot(t_wavelet,wavelet,'b','linewidth',0.7);
title('Ricker wavelet');
xlabel('\itt/s');
set(gca,'FontName','Times New Roman','FontSize',18,'linewidth',1.5);
set(gca,'TickLength',[0 0.0001]);

figure;
imagesc(X,T,Apparent_velocity);
colormap(gray);
ch=colorbar('YTickLabel',{(0:2000:8000),'>10000'});
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
set(get(ch,'title'),'position',[-10 245],'string','Apparent velocity[m/s]');
caxis([0,10000])

figure;
wigb(reflection_cut,1,X_cut,T);
title('reflection coefficient');
xlabel('Offset(m)','FontName','Times New Roman');
ylabel('Time(ms)','FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',18);
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
