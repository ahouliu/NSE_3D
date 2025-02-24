%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  high-order filter rule
%
%
function [UW] = nonProj_func( U,W,UW,K,K_over_K22,N,dim,dealias )
N_half = N/2;
%% velocity->vorticity
for ind=1:dim
    index1 = mod(ind,dim)+1;
    index2 = mod(ind+1,dim)+1;
    W{ind} = 1i*( K{index1}.*U{index2}-K{index2}.*U{index1} );
    W{ind} = W{ind}.*dealias;
end
%% spectral space -> physical space
for ind=1:dim
    U{ind} = real(ifftn(U{ind}));
    W{ind} = real(ifftn(W{ind}));
end
%% physical space -> spectral space
for ind=1:dim
    index1 = mod(ind,dim)+1;
    index2 = mod(ind+1,dim)+1;
    temp = fftn( U{index1}.*W{index2} - U{index2}.*W{index1} );
    temp(1,1,1) = 0;
    temp(N_half+1,:,:) = 0;
    temp(:,N_half+1,:) = 0;
    temp(:,:,N_half+1) = 0;
    UW{ind} = temp;
end
%% projection process
P_temp = K_over_K22{1}.*UW{1}+K_over_K22{2}.*UW{2}+K_over_K22{3}.*UW{3};
for ind=1:dim
    UW{ind} = UW{ind} - K{ind}.*P_temp;
end
end