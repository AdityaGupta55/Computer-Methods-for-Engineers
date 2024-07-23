% A-->B-->C
% B-->D


k1 =2; %rate constant of A-->B reaction                     %initial conditions
k2=0.5;%rate constant of B-->C reaction
k3=0.3;%rate constant of B-->D reaction
h=0.1;
t_final=10;%final time of the set we are taking to get readings
t_initial=0;%initial time

%%EXPLICITING1 PART1
N_explicit= round((t_final-t_initial)/h);%number of point taken in grid here is N_explicit
A_explicit =zeros(N_explicit+1,1);  %Storing value of Ca_explicit          %store concentrations
B_explicit =zeros(N_explicit+1,1);   %Storing value of Cb_explicit
C_explicit =zeros(N_explicit+1,1);%Storing value of Cc_explicit
D_explicit =zeros(N_explicit+1,1);%Storing value of Cd_explicit
t =zeros(N_explicit+1,1);%Matrix storing time of reaction

%%intial conditions 
A_explicit(1)=1;
B_explicit(1)=0;
C_explicit(1)=0;
D_explicit(1)=0;

for i=1:N_explicit                %explicit euler algorithm
    t(i+1)=t(i)+h;
    A_explicit(i+1)=A_explicit(i) + h*(-k1*A_explicit(i));%applying explicit y(i+1)=y(i)+h*F(xi,yi)
    B_explicit(i+1)=B_explicit(i) + h*(k1*A_explicit(i)-k2*B_explicit(i)-k3*B_explicit(i));
    C_explicit(i+1)=C_explicit(i) + h*(B_explicit(i)*k2);
    D_explicit(i+1)=D_explicit(i) + h*(B_explicit(i)*k3);
end
figure(1)      
plot(t,A_explicit,"-r")%plotting the graphs between time and concentration of A,B,C and D for following reaction
hold on
plot(t,B_explicit,"-b",t,C_explicit,"-g",t,D_explicit,"-m")
title('explicit method')          %graph for explicit (part1)



%%implicit (using jacobian )  PART2
N_implicit= ceil((t_final-t_initial)/h);%N_implicit is no. of grid points in case of implicit
t=zeros(N_implicit+1,1);
Y=zeros(4,N_implicit+1);

%%intial conditions 
Y(1,1)=1;         %Y(1,:) is showing conc of A at diff i
Y(2,1)=0;         %Y(2,:) is showing conc of b at diff i
Y(3,1)=0;         %Y(3,:) is showing conc of c at diff i   Y----> implicit
Y(4,1)=0;         %Y(4,:) is showing conc of d at diff i 
X=zeros(4,1);     %Matrix for storing 
t(1)=0;

for i =2:N_implicit+1                   %implicit euler algorithm
     Y(:,i)=Y(:,i-1);
     t(i)=t(i-1)+h;
    for k=1:10%for implicit functiion y(i+1)=y(i)+h*f(x(i+1),y(i+1))    
        j=give_jacobian(h,k1,k2,k3);%calling the funtion of jacobian
        r=give_residual(Y,i,h,k1,k2,k3);%calling the function of residual
        Y(:,i) = Y(:,i)- inv(j)*(r);%applying newton raphson method for Y
    end 
     
end



figure(2)%plotting impicit concentration vs time graph for a,b,c and d.
plot(t',Y(1,:),"-r")      %graph for implict part2
hold on 
plot(t',Y(2,:),"-b",t',Y(3,:),"-g",t',Y(4,:),"-m")
title('Implicit method')




% %FUNCTIONS
% For your easy reference written function used here 
% It is commented here but written in code below to avoid syntax error

% function R = give_residual(Y,i,h,k1,k2,k3)
% 
%   R=[Y(1,i)-Y(1,i-1)+h*(k1)*Y(1,i);
%      Y(2,i)-Y(2,i-1)-h*(k1*Y(1,i)-k2*Y(2,i)-k3*Y(2,i));  
%      Y(3,i)-Y(3,i-1)-h*(k2)*Y(2,i);
%      Y(4,i)-Y(4,i-1)-h*k3*Y(2,i)];
% 
% end
% 
% function j = give_jacobian(h,k1,k2,k3)
% 
%  j=[1+(k1)*h,0,0,0;
%     -h*k1,1+k2*h+k3*h,0,0;
%     0,-h*k2,1,0;
%     0,-h*k3,0,1];
% end


%%rk4 part 3
N_RK4=ceil(t_final-t_initial)/h;
y=zeros(4,N_RK4+1); 
y(1,1)=1;           %y(1,:) showing conc of a at different i
y(2,1)=0;           %y(2,:) showing conc of b at different i
y(3,1)=0;           %y(3,:) showing conc of c at different i   y---->rk4
y(4,1)=0;           %y(4,:) showing conc of d at different i


for i=1:N_RK4           %RK4 algorithm
    y(1,i+1)=y(1,i)+ h*f_A(i,h,y);
    y(2,i+1)=y(2,i)+ h*f_B(i,h,y);
   
    y(3,i+1)=y(3,i)+ h*f_C(i,h,y);
    y(4,i+1)=y(4,i)+ h*f_D(i,h,y);
end


figure(3)%Rk_4 figure 
plot(t',y(1,:),"-r")     %graph for rk4
hold on 
plot(t',y(2,:),"-b",t',y(3,:),"-g",t',y(4,:),"-m")
title('RK4')


%function of A 
function k_A = f_A(i,h,y)
    v1_A= g_A(y(1,i));%for Rk_4 y(i+1)=y(i)+h/6(k1+2*k2+2*k3+k4)
    v2_A= g_A(y(1,i)+h/2*v1_A);%k1=f(x(i),y(i)),k2=f(x(i)+h/2,y(i)+h*k1/2)
    v3_A= g_A(y(1,i)+h/2*v2_A);%k3=f(x(i)+h/2,y(i)+h*k2/2)
    v4_A= g_A(y(1,i)+h/2*v3_A);%k4=f(x(i)+h,y(i)+h*k3)
    k_A = (v1_A+2*v2_A+2*v3_A+v4_A)/6;
end

function x_A=  g_A(A)
    k1=2;
    x_A = -k1*(A);
end

% function of b
function k_B = f_B(i,h,y)

    v1_B= g_B(y(1,i),y(2,i));
    v2_B= g_B((y(1,i)),(y(2,i)+h/2*v1_B));
    v3_B= g_B((y(1,i)),(y(2,i)+h/2*v2_B));
    v4_B= g_B((y(1,i)),(y(2,i)+h/2*v3_B));
    k_B = (v1_B+2*v2_B+2*v3_B+v4_B)/6;

end

function x_B=  g_B(A,B)
    k1=2;
    k2=0.5;
    k3=0.3;
    x_B = k1*(A)-k2*(B)-k3*(B);
end

%function of c
    function k_C = f_C(i,h,y)
    v1_C= g_C(y(2,i));
    v2_C= g_C(y(2,i));
    v3_C= g_C(y(2,i));
    v4_C= g_C(y(2,i));
    k_C = (v1_C+2*v2_C+2*v3_C+v4_C)/6;
    end


function x_C=  g_C(A)
    k2=0.5;
    x_C= k2*(A);
end


%function of d
function k_D = f_D(i,h,y)
    v1_D= g_D(y(2,i));
    v2_D= g_D(y(2,i));
    v3_D= g_D(y(2,i));
    v4_D= g_D(y(2,i));
    k_D = (v1_D+2*v2_D+2*v3_D+v4_D)/6;
end

function x_D=  g_D(A)
    k3=0.3;
    x_D= k3*(A);
end


%functions for implicit part2

function R = give_residual(Y,i,h,k1,k2,k3)%function for making matrix of the residuals
%Residual is difference between response data and fix to response data at
%each predictor value
  R=[Y(1,i)-Y(1,i-1)+h*(k1)*Y(1,i);
     Y(2,i)-Y(2,i-1)-h*(k1*Y(1,i)-k2*Y(2,i)-k3*Y(2,i));
     Y(3,i)-Y(3,i-1)-h*(k2)*Y(2,i);
     Y(4,i)-Y(4,i-1)-h*k3*Y(2,i)];

end

function j = give_jacobian(h,k1,k2,k3)%function form for writing jacobian in the form of 4by4 matrix

 j=[1+(k1)*h,0,0,0;
    -h*k1,1+k2*h+k3*h,0,0;
    0,-h*k2,1,0;
    0,-h*k3,0,1];
end