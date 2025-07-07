clear all
close all
clc

tic

nit_sim = 5400/0.3;

Lyapbasis = load("U_D_collection.txt","-ascii");
Lyapbasis = Lyapbasis(1:(nit_sim*168506),:);

Rfull = load("R_full.txt","-ascii");
Rfull = Rfull(1:(nit_sim*6),:);

% NonL = load("NonL.txt","-ascii");
% NonL = NonL(1:(nit_sim*168506));

%% Backward evolution to form C_N

N = 84253*2;
m = 6;
stept = 0.3;
% nit_sim = size(NonL,1)/(N); % Based on data.mat
T = nit_sim*stept;
K = T/abs(stept);
k1 = 2100/0.3;
k2 = 4200/0.3;


% Intialisation of C_K
detC = 0;
while detC == 0
    C = triu(rand(m,m));
    for i = 1:m
        C(i,i) = 1;
    end
    detC = det(C);
end
C_n = zeros(m,m,T/abs(stept));
C_n(:,:,end) = C;
for i=1:m
    C_n(:,i,K) = C_n(:,i,K)/norm(C_n(:,i,K)); % normalise C by columns
end

% Refining Rfull
for(i=1:nit_sim)

    Rfull((i-1)*6 + 2,1) = 0;
    Rfull((i-1)*6 + 3,1) = 0;
    Rfull((i-1)*6 + 4,1) = 0;
    Rfull((i-1)*6 + 5,1) = 0;
    Rfull((i-1)*6 + 6,1) = 0;

    Rfull((i-1)*6 + 3,2) = 0;
    Rfull((i-1)*6 + 4,2) = 0;
    Rfull((i-1)*6 + 5,2) = 0;
    Rfull((i-1)*6 + 6,2) = 0;

    Rfull((i-1)*6 + 4,3) = 0;
    Rfull((i-1)*6 + 5,3) = 0;
    Rfull((i-1)*6 + 6,3) = 0;    

    Rfull((i-1)*6 + 5,4) = 0;
    Rfull((i-1)*6 + 6,4) = 0; 

    Rfull((i-1)*6 + 6,5) = 0;
end    

for i=1:nit_sim
    R_full(:,:,i) = Rfull(1+(i-1)*m:m*i,:);
    Lyap_b(:,:,i) = Lyapbasis(1+(i-1)*N:N*i,:);
end


% Forming the matrices of coefficients C
for i=1:(K-k1)
    C_n(:,:,end-i) = R_full(:,:,end-(i-1))\C_n(:,:,end-(i-1));
    C_n(:,:,end-i) = normalize(C_n(:,:,end-i),1,'norm');
end

% Forming the CLVs now
v = zeros(N,m,K);
for t = k1:k2
    for i=1:m
        for j=1:i
            v(:,i,t) = v(:,i,t) + C_n(j,i,t)*Lyap_b(:,j,t);
        end
    end
end

v1 = v(:,1,:); 
CLV1 = reshape(v1,[N nit_sim]);
CLV1x = CLV1(1:2:N,:);
CLV1y = CLV1(2:2:N,:);
writematrix(CLV1x,'CLV1x.txt','Delimiter','tab');
writematrix(CLV1y,'CLV1y.txt','Delimiter','tab');
clear CLV1x CLV1y;

v2 = v(:,2,:); 
CLV2 = reshape(v2,[N nit_sim]);
CLV2x = CLV2(1:2:N,:);
CLV2y = CLV2(2:2:N,:);
writematrix(CLV2x,'CLV2x.txt','Delimiter','tab');
writematrix(CLV2y,'CLV2y.txt','Delimiter','tab');
clear CLV2x CLV2y;

v3 = v(:,3,:); 
CLV3 = reshape(v3,[N nit_sim]);
CLV3x = CLV3(1:2:N,:);
CLV3y = CLV3(2:2:N,:);
writematrix(CLV3x,'CLV3x.txt','Delimiter','tab');
writematrix(CLV3y,'CLV3y.txt','Delimiter','tab');
clear CLV3x CLV3y;

% v4 = v(:,4,:); 
% CLV4 = reshape(v4,[N nit_sim]);
% CLV4x = CLV4(1:2:N,:);
% CLV4y = CLV4(2:2:N,:);
% writematrix(CLV4x,'CLV4x.txt','Delimiter','tab');
% writematrix(CLV4y,'CLV4y.txt','Delimiter','tab');
% clear CLV4x CLV4y;
% 
% v5 = v(:,5,:); 
% CLV5 = reshape(v5,[N nit_sim]);
% CLV5x = CLV5(1:2:N,:);
% CLV5y = CLV5(2:2:N,:);
% writematrix(CLV5x,'CLV5x.txt','Delimiter','tab');
% writematrix(CLV5y,'CLV5y.txt','Delimiter','tab');
% clear CLV5x CLV5y;

v6 = v(:,6,:); 
CLV6 = reshape(v6,[N nit_sim]);
CLV6x = CLV6(1:2:N,:);
CLV6y = CLV6(2:2:N,:);
writematrix(CLV6x,'CLV6x.txt','Delimiter','tab');
writematrix(CLV6y,'CLV6y.txt','Delimiter','tab');
clear CLV6x CLV6y;

toc
