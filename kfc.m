function [A, B]=kfc(y, u, t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman Filter Based Online System Identification and Denoising          %
%                                                                         %
% Modify Q (line:42) and R (line:47) for your model to get better results %
%                                                                         %
% Author: Reza Samarghandi                                                %
%                                                                         %
% Email: Rezasamargandi@yahoo.com                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



n=size(y,1);
m=size(u,1);


if t==0
    I0=eye(n,n+m);
    It0=I0.';
    I0v = It0(:);
    theta=I0v; %state variable vector 
    gam=1000; %large value
    Pk=gam.*eye(n*(n+m)); %covariance matrix
end




X=[y; u];

Xh=X.';


thetak1=theta;



H= kron(eye(n),Xh);

Q=H.'*H;


Pk1=Pk+Q;

R=cov((y-H*thetak1).');

S=H*Pk1*H.'+R; %noise covariance matrix

Kk=Pk1*H.'*S^-1; %Kalman gain


deltheta=Kk*(y-H*thetak1); %modification vector for state variables 





theta=thetak1+deltheta;

Pk=(eye(n*(n+m))-Kk*H)*Pk1;


Ap=reshape(theta,[n+m n]).'; %[A B]


A=Ap(:,1:n); %A Matrix of State Space 
B=Ap(:,n+1:end); %B Matrix of State Space
