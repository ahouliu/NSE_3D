%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Fourier-Galerkin, NS on [0,2\pi]^3, tri-periodic.
%    dealiasing: high (36th) order Fourier smoothing 
%       via <<Computing Nearly Singular Solutions Using Pseudo-Spectral Methods>>
%
%    INPUT:
%      paras = { L,dim,N,Re,K,K22,K_over_K22 ,dealias}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,Ek,U_hat] = dns_3d_f_1_func(paras,c1,c2,d,e,U_hat)
tic;
L = paras{1};
N = paras{3};
nu = 1/paras{4};
dim = paras{2};
t0 = c1;
MAX_time = c2;
dt = d;
Ninteration = floor((MAX_time-t0)/dt); 
if (MAX_time-t0-dt*Ninteration)~=0
    error('NOT MATCH?');
end
%% spectral space
epa = e;
K = paras{5};
K22 = paras{6};
K_over_K22 = paras{7};
dealias = paras{8};
% external forces
% F_hat_proj = { zeros(N,N,N), zeros(N,N,N), zeros(N,N,N) };
F_hat_proj = force_proj_func(epa,L,N,dim,K,K_over_K22);
clear('paras');
%% collecting data
Ek = zeros(Ninteration+1,1);
p_DFT = 1/(N^3);
Ek(1) = 0.5*L^2*sum(sum(sum( abs(U_hat{1}).^2 + abs(U_hat{2}).^2 + abs(U_hat{3}).^2 )))*p_DFT;
%% time evoluting 
RK_stage = 5;
RK_A = [            0;              -567301805773/1357537059087; -2404267990393/2016746695238; -3550918686646/2091501179385; -1275806237668/842570457699  ];
RK_B = [1432997174477/9575080441755; 5161836677717/13612068292357; 1720146321549/2090206949498;  3134564353537/4481467310338;  2277821191437/14882151754819];
RK_c = [            0;              1432997174477/9575080441755;  2526269341429/6820363962896;  2006345519317/3224310063776;  2802321613138/2924317926251 ];
q_star = { zeros(N,N,N),zeros(N,N,N),zeros(N,N,N) };
U      = { zeros(N,N,N),zeros(N,N,N),zeros(N,N,N) };
W      = { zeros(N,N,N),zeros(N,N,N),zeros(N,N,N) }; 
UW     = { zeros(N,N,N),zeros(N,N,N),zeros(N,N,N) }; 

for ind = 1:Ninteration

	% (5,4) R-K
    %q_star = q_star*0;
    for jnd=1:RK_stage
    	Hconv_proj = nonProj_func(U_hat,W,UW,K,K_over_K22,N,dim,dealias);
    	t_temp = t0 + dt*( ind-1 + RK_c(jnd) );
    	q_star{1} = RK_A(jnd)*q_star{1} + ...
    							   dt*( Hconv_proj{1} - nu*K22.*U_hat{1} + F_hat_proj{1}*g_func(t_temp) );
    	q_star{2} = RK_A(jnd)*q_star{2} + ...
    							   dt*( Hconv_proj{2} - nu*K22.*U_hat{2} + F_hat_proj{2}*g_func(t_temp) );
    	q_star{3} = RK_A(jnd)*q_star{3} + ...
    							   dt*( Hconv_proj{3} - nu*K22.*U_hat{3} + F_hat_proj{3}*g_func(t_temp) );
                        
        U_hat{1}   = U_hat{1} + RK_B(jnd)*q_star{1};
        U_hat{2}   = U_hat{2} + RK_B(jnd)*q_star{2};
        U_hat{3}   = U_hat{3} + RK_B(jnd)*q_star{3};
    end

    N2 = sum(sum(sum( abs(U_hat{1}).^2 + abs(U_hat{2}).^2 + abs(U_hat{3}).^2 ))); 
    if isnan(N2)
        error('Explode at t=%g (dt=%g).\n ',t0+dt*ind,dt);
    end
    Ek(1+ind) = 0.5*L^2*N2*p_DFT;

end

for ind=1:dim
    U{ind} = real(ifftn(U_hat{ind}));
end

toc
end