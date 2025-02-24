function [F_hat] = force_proj_func(epa,L,N,dim,K,K_over_K22)
%%  physical space
    [X,Y,Z] = meshgrid( (0:N-1)*(L/N), (0:N-1)*(L/N), (0:N-1)*(L/N) );
%%  spectral space  
    F_hat = {zeros(N,N,N),zeros(N,N,N),zeros(N,N,N)};
    forceAmp = (4*sqrt(3)/9)^3;
    k0 = 4;
    for jnd=1:dim
        F_hat{jnd}                   = fftn( epa*forceAmp*sin(k0*X).*sin(k0*Y).*sin(k0*Z).*(1-cos(k0*X)).*(1-cos(k0*Y)).*(1-cos(k0*Z)) );
        F_hat{jnd}(abs(F_hat{jnd})<1e-9) = 0;
        index                        = find(abs(real(F_hat{jnd}))>0);
        F_hat{jnd}(index)            = 1i*imag(F_hat{jnd}(index));
    end
%%  projection process
    P_temp = K_over_K22{1}.*F_hat{1}+K_over_K22{2}.*F_hat{2}+K_over_K22{3}.*F_hat{3};
    for jnd=1:dim
        F_hat{jnd} = F_hat{jnd} - K{jnd}.*P_temp;
    end
end