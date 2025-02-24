%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
%    U_hat(k1,k2,k3), k3=0,1,...,N/2-1 -> U_hat(k1,k3,k3), k3=-N/2+1,...,-1,0,1,...,N/2-1,N/2

function [U1_hat] = half2full_func(U_hat,N,N_half)

    if size(U_hat,3)~=N_half
        error('not match.');
    end

    U1_hat                           = zeros(N,N,N);
    U1_hat(:,:,1:N_half)             = U_hat(:,:,1:N_half);
    U1_hat(1,1,N_half+2:end)         = conj(       flip(  U_hat(1,1,2:N_half),3)  );
    U1_hat(2:end,1,N_half+2:end)     = conj(  flip(flip(  U_hat(2:end,1,2:N_half),3),1)      );
    U1_hat(1,2:end,N_half+2:end)     = conj(  flip(flip(  U_hat(1,2:end,2:N_half),3),2)      );
    U1_hat(2:end,2:end,N_half+2:end) = conj(  flip(rot90( U_hat(2:end,2:end,2:N_half),2),3)  );

end