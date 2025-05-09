function [f,amplitude_spectrum]=Amplitude_spectrum_my(dt,signal,fmax)
% dt;%ʱ�����/s
% signal;%����ʱ���ź�
% fmax;%�����ʾƵ��/Hz
% f;%���Ƶ������/Hz
% amplitude_spectrum;%�����ź�signal�������
F=fft(signal,2^ceil(log2(1000)));
N=length(F);
amplitude_spectrum=abs(F);%abs(Cn)=abs(X(k))=sqrt(realX(k)^2+imagX(k)^2)
f=(0:N-1)*(1/(N*dt));%Ƶ�ʲ���������ڻ���Ƶ��f0=1/T;T=dt*N=N/fs
fmax_number=ceil(fmax/(1/(N*dt)));
f=f(1:fmax_number);
amplitude_spectrum=amplitude_spectrum(1:fmax_number)*2/N;%Ƶ������ʵ���֮��Ĺ�ϵ��������ף�����abs(Cn)=abs(Xn)=1/2*An,������ɢʱ�����и���Ҷ�任��X(n)��X(w)ʱ��Ϊ�˷�����X(w)ǰ����1/N;���Ϊ�˻ָ���Ƶ�ʶ�Ӧ���������
                              %Ӧ��mag=abs(y)����2/N; ����ʵ�����У���ɢʱ�����еĸ���Ҷ�任��X(k)����N/2���жԳ���;���ֻҪд��NyquistƵ��֮ǰ��N/2��Ƶ�ʵ����
