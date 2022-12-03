function [ index ] = s_opt_generator( Vo, M, index, outfile )
%
% s_opt_generator.m - Generating S-optimal Points indices of a matrix
%
% Syntax     s_opt_generator( A, M, index, outfile )
%
% Input:     A       = Candidate Matrix
%            M       = the number of S-optimal points to be computed,
%            number of rows basically
%            idx     = if there is a previous index set, put it in. if not,
%                      just put [].
%            outfile = put 'filename', then 'filename' will be created with
%            the information of index.
% 
% Output:    index = S-optimal points index set, of size (M x 1)
%
% 
% Written by Yeonjong Shin   
% Last edited on 05/20/2022
% Email : shin.481@osu.edu / yeonjong_shin@brown.edu
%--------------------------------------------------------------------------
[N_Bm, N] = size(Vo);
nVo = Vo.*Vo;

inum = length(index);
if  inum == 0
    [~, i_idx] = max(abs(Vo(:,1)));
    index = [];
    index = [index; i_idx];
    inum = inum + 1;
end

for i = inum+1 : M
    if i < N+1
        V1 = Vo(index, 1:i);
        A0 = V1(1:i-1,1:i-1); atA0 = V1(1:i-1,i)'*A0;
        tt = zeros(N_Bm,i-1); tt1 = zeros(i-1,N_Bm);
        ata = sum(V1(1:i-1,i).^2) ; 
        bbb = (A0'*A0)\[atA0; Vo(:,1:i-1)]'; c = bbb(:,2:end);
        b = 1 + sum(Vo(:,1:i-1).*c',2);
        
        for zz = 1 : i-1
            tt(:,zz) = Vo(:,zz).*Vo(:,i);
        end
        g1 = repmat(atA0,N_Bm,1) + tt; g1 = g1';
        g2 = bbb(:,1);
        oneprc = 1 + sum(Vo(:,1:i-1).*c',2);
        g3 = sum(c'.*g1',2)./oneprc;
        for zz = 1 : i-1
            tt1(zz,:) = c(zz,:).*(Vo(:,i) - g3)';
        end
        
        GG = repmat(g2,1,N_Bm) + tt1; 

        A = ata + Vo(:,i).^2 - sum(g1'.*GG',2);
        A(A<0) = 0;
        nV = sum(nVo(index,1:i),1);
        noM = sum(log(repmat(nV,N_Bm,1) + nVo(:,1:i)),2);        
        
        A = log(abs(A)) + log(b) - noM;
    else
        V1 = Vo(index, :);
        b = (V1'*V1)\Vo';
        nV = sum(nVo(index,:));  
        noM = sum(log(repmat(nV,N_Bm,1) + nVo),2);
        A = log(1+sum(Vo.*b',2)) - noM;
    end
    A(index) = -inf;
    [~, ggg] = max(A);
    index = [index; ggg(1)];

    if i == floor(M/4) || i == floor(M/2) || i == floor(3/4*M)
        fprintf('i = %d, index = %d, N = %d, M = %d\n', i, ggg(1), N, M);
    elseif i == M
        fprintf('Done\n');
    end
    
end

if isempty(outfile) == 0
    sfile = sprintf('%s%s',outfile,'.mat');
    save(sfile,'index');
end
    
end
