clear;close all;clc;
%% pre-stage
dim    = 3;
N      = 2^6;
N_half = N/2;
L      = 2*pi;
Re     = 40; 
nu     = 1/Re;
%% physical space
[X,Y,Z] = meshgrid( (0:N-1)*(L/N), (0:N-1)*(L/N), (0:N-1)*(L/N) );
P       =   zeros(N,N,N);
%F       = force_func(X,Y,Z,epa);
%% spectral space
kx               = [0:N_half-1, 0, -N_half+1:-1]';
K                = { zeros(N,N,N),zeros(N,N,N),zeros(N,N,N) };
[K{1},K{2},K{3}] = meshgrid(kx,kx,kx);
K22              = K{1}.^2 + K{2}.^2 + K{3}.^2;
K22_1            = 1./K22;
K22_1(K22==0)    = 0;
K_over_K22       = { zeros(N,N,N),zeros(N,N,N),zeros(N,N,N) };
for ind =1:dim
    temp              = K{ind}./K22;
    temp(isnan(temp)) = 0;
    K_over_K22(ind)   = mat2cell(temp,N);
end
%% dealiasing matrix
alpha = 36;beta = 36;
dealias = @(x,y,z) exp(-alpha*(x).^beta).*exp(-alpha*(y).^beta).*exp(-alpha*(z).^beta);
dealias = dealias(K{1}/N_half,K{2}/N_half,K{3}/N_half);
dealias(N_half+1,:,:) = 0;
dealias(:,N_half+1,:) = 0;
dealias(:,:,N_half+1) = 0;
clear('kx','alpha','beta');
%% initial conditions
% %---------------------  1. Beltrami flows  ------------------------------------------
% A=1;B=2;C=3;k=2;
% t=0;
% u=@(x,y,z) exp(-nu*k^2*t)*( A*sin(k*z)+C*cos(k*y) );
% v=@(x,y,z) exp(-nu*k^2*t)*( B*sin(k*x)+A*cos(k*z) );
% w=@(x,y,z) exp(-nu*k^2*t)*( C*sin(k*y)+B*cos(k*x) );
% U    = { zeros(N,N,N),zeros(N,N,N),zeros(N,N,N) };
% U{1} = u(X,Y,Z);
% U{2} = v(X,Y,Z);
% U{3} = w(X,Y,Z);
% clear('a','d','t','u','v','w');
% U0_hat = {fftn(U{1}),fftn(U{2}),fftn(U{3})};
% for ind = 1:dim
%     U0_hat{ind} = U0_hat{ind}.*dealias;
%     U0_hat{ind}(1,1,1) = 0;
% end
%---------------------  2. random  ------------------------------------------
U0_hat = { zeros(N,N,N),zeros(N,N,N),zeros(N,N,N) }; 
c1=rand()*2*pi;
c2=rand()*2*pi;
c3=rand()*2*pi;

kf      = 3;
Ek_hat  = @(k1,k2,k3) (9/11/kf)*( (sqrt(k1.^2+k2.^2+k3.^2)<=kf).*(k1.^2+k2.^2+k3.^2)./(kf^2)+ (sqrt(k1.^2+k2.^2+k3.^2)>kf).*((sqrt(k1.^2+k2.^2+k3.^2)/kf+eps).^(-5/3)) );
Ek_hat  = Ek_hat(K{1}(:,:,1:N_half),K{2}(:,:,1:N_half),K{3}(:,:,1:N_half));
U0_hat{1}(:,:,1:N_half) = sqrt(Ek_hat.*K22_1(:,:,1:N_half)/(2*pi)).*   K{3}(:,:,1:N_half)*cos(c1)*exp(1i*c2);
U0_hat{2}(:,:,1:N_half) = sqrt(Ek_hat.*K22_1(:,:,1:N_half)/(2*pi)).*   K{3}(:,:,1:N_half)*sin(c1)*exp(1i*c3);
U0_hat{3}(:,:,1:N_half) = sqrt(Ek_hat.*K22_1(:,:,1:N_half)/(2*pi)).*( -K{1}(:,:,1:N_half)*cos(c1)*exp(1i*c2)-K{2}(:,:,1:N_half)*sin(c1)*exp(1i*c3)  );
for ind=1:dim
    U0_hat{ind} = half2full_func(U0_hat{ind}(:,:,1:N_half),N,N_half);
    U0_hat{ind} = U0_hat{ind}.*dealias;
end
clear('c1','c2','c3','kf','Ek_hat');
%% main dns process
paras = { L,dim,N,Re,K,K22,K_over_K22, dealias};
t0 = 0;
MAX_time = 200;
dt     = 0.01;
epa    = 1;  % equal-strength forces
sample_interval = 0.1;
% [U1,Ek1,U_hat1] = dns_3d_f_1_func(paras,t0,MAX_time,dt,epa,U0_hat);
% [EE1,Ek1,U_hat1] = dns_3d_f_proj_func(paras,t0,MAX_time,dt,epa,U0_hat,QQ,sample_interval);
% [UVW1,Ek1,U_hat1] = dns_3d_f_2_func(paras,t0,MAX_time,dt,epa,U0_hat,sample_interval);
[U_data1,Ek1,U_hat1] = dns_3d_f_full_func(paras,t0,MAX_time,dt,epa,U0_hat,sample_interval);

% %% test with benchmark 
% A=1;B=2;C=3;k=2;
% t=MAX_time;
% u=@(x,y,z) exp(-nu*k^2*t)*( A*sin(k*z)+C*cos(k*y) );
% v=@(x,y,z) exp(-nu*k^2*t)*( B*sin(k*x)+A*cos(k*z) );
% w=@(x,y,z) exp(-nu*k^2*t)*( C*sin(k*y)+B*cos(k*x) );
% U0={zeros(N,N,N),zeros(N,N,N),zeros(N,N,N)};
% U0{1} = u(X,Y,Z);
% U0{2} = v(X,Y,Z);
% U0{3} = w(X,Y,Z);
% norm(reshape(U{1},N^3,1)-reshape(U0{1},N^3,1))+norm(reshape(U{2},N^3,1)-reshape(U0{2},N^3,1))+norm(reshape(U{3},N^3,1)-reshape(U0{3},N^3,1))
%% Saving AVI/MP4-- stream tubes
t_start0 = 0;
fielName = 'U_streamtubes_test';
figure(1);
temp_size = get(0,'screensize');
set(gcf,'outerposition',[temp_size(1) temp_size(1) temp_size(4) temp_size(4)],...
        'color',[1 1 1],...
        'Toolbar','none');
videoMaker = VideoWriter(fielName,'MPEG-4');
open(videoMaker);
formatSpec = '%.1f';
[sx,sy,sz] = meshgrid( [0,(L/N)*(N/6),(L/N)*(N/3),(L/N)*(N/2),(L/N)*(2*N/3),(L/N)*(5*N/6)], ...
                       [0,(L/N)*(N/6),(L/N)*(N/3),(L/N)*(N/2),(L/N)*(2*N/3),(L/N)*(5*N/6)], ...
                       [0,(L/N)*(N/6),(L/N)*(N/3),(L/N)*(N/2),(L/N)*(2*N/3),(L/N)*(5*N/6)] );
LEN = size(U_data{1},4);
for ind=1:LEN
    cla;
    streamtube(X,Y,Z,U_data{1}(:,:,:,ind),U_data{2}(:,:,:,ind),U_data{3}(:,:,:,ind),sx,sy,sz);
    view(3)
    axis tight
    shading interp;
    camlight;
    lighting gouraud
    colorbar;
    set(colorbar,'TickLength',[0]);
    set(gca,'LineWidth',2,...
            'FontName','Times New Roman',...
            'FontSize',22,...
            'TickLength',[0,0]);
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    title(['$t=',num2str(t_start0+sample_interval*(ind-1), formatSpec),'$'],'Interpreter','latex');
    axis([0 2*pi 0 2*pi 0 2*pi]);
    frame = getframe(gcf);
	writeVideo(videoMaker,frame);
end
close(videoMaker);
clear('videoMaker','t_start0','Level','fielName','temp_size','max_W','min_W','interval','h','frame');
close(1);
%% Saving AVI/MP4-- slices of |U|
LEN = size(U_data1{1},4);
UVW = zeros(N,N,N,LEN);
for ind=1:LEN
    UVW(:,:,:,ind)=sqrt( U_data1{1}(:,:,:,ind).^2+U_data1{2}(:,:,:,ind).^2+U_data1{3}(:,:,:,ind).^2 );
end

t_start0 = 0;
fielName = 'U_test';
figure(4);
temp_size = get(0,'screensize');
set(gcf,'outerposition',[temp_size(1) temp_size(1) temp_size(4) temp_size(4)],...
        'color',[1 1 1],...
        'Toolbar','none');
videoMaker = VideoWriter(fielName,'MPEG-4');
open(videoMaker);

xslice = 0;   
yslice = 0;
zslice = (L/N)*(N-1);
max_c = 0.1*round(max(UVW,[],'all')*10);
min_c = 0.1*round(min(UVW,[],'all')*10);
if max_c==min_c
    max_c = max(UVW,[],'all');
    min_c = min(UVW,[],'all');
end
formatSpec = '%.1f';

for ind=1:LEN
    cla;
    slice(X,Y,Z,UVW(:,:,:,ind),xslice,yslice,zslice,'cubic');
    grid off;
    view(3)
    axis tight
    shading interp;
    % camlight;
    lighting gouraud
    colorbar;
    set(colorbar,'TickLength',[0]);
    set(gca,'LineWidth',2,...
            'FontName','Times New Roman',...
            'FontSize',22,...
            'TickLength',[0,0]);
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    title(['$|U|,~t=',num2str(t_start0+sample_interval*(ind-1), formatSpec),'$'],'Interpreter','latex');
    caxis([min_c,max_c]);
    frame = getframe(gcf);
	writeVideo(videoMaker,frame);
end
close(videoMaker);
clear('videoMaker','t_start0','fielName','temp_size','max_W','min_W','interval','h','frame','xslice','yslice','zslice');
close(4);