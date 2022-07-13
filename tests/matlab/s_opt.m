%% S-OPT Example for 1D LSROM ---------------------------------------------
%
% Written by Y.Shin
% Last edited: 2022/05/20
% Email: yeonjong_shin@brown.edu
%
% Load data ---------------------------------------------------------------
load('../s_opt_data/1D_LSROM_basis5.mat'); % The basis matrix phi_r is loaded

% QR mode  ----------------------------------------------------------------
QR_mode = 'off' ;   % if 'on' , Q_r is obtained from QR-factorization of phi_r
                    % if 'off', Q_r = phi_r

% Input -------------------------------------------------------------------
num_col     = 5;  % Choose the number of basis to be used; Ncol <= num_col
num_samples = 10; % Choose the number of samples to generate.

% File Name for the S-OPT index -------------------------------------------
switch QR_mode
    case 'on'
        [Q_r, ~] = qr(phi_r(:,1:num_col),0);
        str_save    = strcat('../s_opt_data/1D_LSROM_basis_QRindex_Col',num2str(num_col)); 
    case 'off'
        Q_r = phi_r(:,1:num_col);
        str_save    = strcat('../s_opt_data/1D_LSROM_basis_index_Col',num2str(num_col)); 
end
[num_row,~] = size(Q_r);

% S-OPT calucation --------------------------------------------------------
Q_sampled = s_opt_generator(Q_r,num_samples,[],str_save);
Q_sampled = sort(Q_sampled);
inv_Q = pinv(Q_r(Q_sampled,:))';
