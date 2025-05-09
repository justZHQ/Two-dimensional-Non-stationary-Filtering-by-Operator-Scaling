function [t,ricker,f,amplitude_spectrum]=Ricker_my(dt,L,rick_mainf,fa,fmax)
% dt;%ʱ�����/s
% L;%�Ӳ����ȣ�%wavelength=2*L+1
% rick_mainf;%�׿��Ӳ�����Ƶ/HZ
% fa;%��λ
% fmax;%�����ʾƵ��/Hz
% t;%���ʱ������/s
% ricker;%���ʱ�����׿��Ӳ�
% f;%���Ƶ������
% amplitude_spectrum;%����׿��Ӳ��������

t=-L*dt:dt:L*dt;
s1=(1-2*(pi*rick_mainf*t).^2).*exp(-(pi*rick_mainf*t).^2);
s2=hilbert(s1);
ricker=real(s2)*cos(fa)+imag(s2)*sin(fa);
[f,amplitude_spectrum]=Amplitude_spectrum_my(dt,ricker,fmax);