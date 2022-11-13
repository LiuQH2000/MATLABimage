%角谱传播的GS算法
%参数初始化
clc;clear;close all;
M=512;N=512;%设置图像的尺寸
lambda=632.8e-6;%设置波长
z=10000;%设置物平面到衍射平面的距离
k=2*pi/lambda;%设置波矢
iterative=100;%设置迭代次数
piesize=8e-3;%像素大小
L=N*piesize;%长宽
%若对图像添加噪声，使用imnoise（I,type,parameters）函数，I为原始图像。
%type为噪声类型（gaussian：高斯噪声，salt & peper：椒盐噪声，poisson 泊松噪声），parameters为噪声系数。
%生成振幅以及读取图像
objectIntensity=im2double(ones(M,N));
objectAmplitude0=abs(sqrt(objectIntensity));%得到幅度值
%给初始强度图像加入椒盐噪声
%objectAmplitude0=imnoise(abs(sqrt(objectIntensity)),'salt & pepper');
%给初始强度图像加入高斯噪声
%objectAmplitude0=imnoise(abs(sqrt(objectIntensity)),'gaussian');
%给初始强度图像加入泊松噪声
objectAmplitude0=imnoise(abs(sqrt(objectIntensity)),'poisson');
figure(1);imshow(objectIntensity,[]);title('原始强度图像');
%读取相位图像
phase=imresize(imread('cameraman.tif'),[M,N]);%
%phase=imnoise(phase,'salt & pepper');%给初始相位图片添加椒盐噪声
%phase=imnoise(phase,'gaussian');%给初始相位图片添加高斯噪声
phase=im2double(phase);%对相位归一化
figure(2);imshow(phase,[M,N]);title('原始相位');
objectAmplitude1=objectAmplitude0.*exp(1i*phase);%得到物平面的复振幅分布
%频域初始化
[x,y,~]=size(objectAmplitude1);
fX=[0:fix(x/2),ceil(x/2)-1:-1:1]./L;
fY=[0:fix(y/2),ceil(y/2)-1:-1:1]./L;
[fx,fy]=meshgrid(fX,fY);
%定义角谱传播函数
H=exp(1i*k*z.*sqrt(1-(lambda*lambda).*(fx.^2+fy.^2)));%得到角谱传递函数
HB=1./H;
%
objectAmplitudeJP=abs(ifft2(fft2(objectAmplitude1).*H));%%得到经过角谱传播之后的像面上复振幅值 已知量
figure(3);imshow(objectAmplitudeJP,[]);title('衍射面图像');
phase=rand(M,N);
objectAmplitude2=objectAmplitude0.*exp(1i.*phase);
%开始迭代
tic
for i=1:iterative
    OA=objectAmplitude2;
    %角谱正衍射
    objectAmplitudeJP2=ifft2(fft2(OA).*H);%得到经过角谱传播之后的像面上复振幅的分布
    phase1=angle(objectAmplitudeJP2);
    objectAmplitudeJP2=objectAmplitudeJP.*exp(1i.*phase1);%得到替换后的像面上的复振幅分布
    %角谱逆衍射
    OA1=ifft2(fft2(objectAmplitudeJP2).*HB);%得到角谱逆衍射的复振幅分布
    phase2=angle(OA1);
    OA1=objectAmplitude0.*exp(1i.*phase2);
    objectAmplitude3=OA1;%迭代完成得到的相位图像
    objectAmplitude2=objectAmplitude3;
    error(i)=RMSE(angle(objectAmplitude3),angle(objectAmplitude1));
end
toc
  huifu_phase=angle(OA1);
  intensity=(abs(objectAmplitude2)).^2;
  figure(4);imshow(huifu_phase,[]);title('迭代完成后图像');
  error_RMSE_I=RMSE(angle(objectAmplitude3),angle(objectAmplitude1));
% RMSE_objectAmplitude=RMSE(angle(objectAmplitude3),angle(objectAmplitude1));
  figure(5);
  x=1:iterative;  
  plot(x,error,'-r','LineWidth', 2);
  legend('角铺迭代');
  xlabel('迭代次数');
  ylabel('RMSE');