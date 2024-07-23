function omega = flow_around_cylinder_unsteady
Re=60;
%%%%% define the grid %%%%%
n=101; m=202; % number of grid points
N=n-1; M=m-2; % number of grid intervals: 2 ghost points, theta=-h,2*pi
h=2*pi/M; % grid spacing based on theta variable
xi=(0:N)*h; theta=(-1:M)*h; % xi and theta variables on the grid
%%%%% time-stepping parameters %%%%%
t_start=0; t_end=0.5; %vortex street starts at around t=1000
tspan=[t_start t_end]; 
%%%%% Initialize vorticity field %%%%%
omega=zeros(n,m);
diagonals = [4*ones(n*m,1), -ones(n*m,4)];
A=spdiags(diagonals,[0 -1 1 -n n], n*m, n*m); 
%set up all the boundary conditions in the matrix
I=speye(n*m);
index=1:n:1+(m-1)*n; %rows for xi=0
A(index,:)=I(index,:); 
index=n:n:m*n; %rows for xi=xi_max
A(index,:)=I(index,:); 
index1=1:n; %indices for theta = -h
index2=1+(m-2)*n:(m-1)*n; %indices for theta = 2*pi-h
A(index1,:)=I(index1,:)-I(index2,:); %periodic: j=1 same as j=m-2
index1=1+(m-1)*n:m*n;%indices for theta=2*pi 
index2=1+n:2*n; %indices for theta = 0
A(index1,:)=I(index1,:)-I(index2,:); %periodic j=m same as j=2
%%%%% Find the LU decomposition %%%%%
[L,U]=lu(A); clear A;
%%%%% Compute any time-independent constants %%%%%
fac1=h^2*repmat(exp(2*xi'),1,m); %psi equation
fac2=(1/(4*h^2))*repmat(exp(-2*xi'),1,m); %omega equation (first term)
fac3= (2/(h^2*Re))*repmat(exp(-2*xi'),1,m); %omega equation (second term)   
%%%%% advance solution using ode23 %%%%%
options=odeset('RelTol', 1.e-03);
omega=omega(2:n-1,2:m-1); % strip boundary values for ode23
omega=omega(:); % make a column vector
[t,omega]=ode23...
      (@(t,omega)omega_eq(omega,n,m,xi,theta,h,fac1,fac2,fac3,L,U),...
                                             tspan, omega, options);%%%%% expand omega to include boundaries %%%%%
temp=zeros(n,m);
temp(2:n-1,2:m-1)=reshape(omega(end,:),n-2,m-2);
omega=temp; clear temp;
%%%%% compute stream function (needed for omega boundary values) %%%%%
b=zeros(n,m); 
b(2:n-1,2:m-1)=fac1(2:n-1,2:m-1).*omega(2:n-1,2:m-1);
b(n,2:m-1)=exp(xi(n))*sin(theta(2:m-1)); %free stream
b=reshape(b,n*m,1);
psi=reshape(U\(L\b),n,m);
%%%%% set omega boundary conditions %%%%%
omega(1,:)=(1/(2*h^2))*(psi(3,:)-8*psi(2,:));
omega(n,:)=0;
omega(:,1)=omega(:,m-1);
omega(:,m)=omega(:,2); 
close; 
%%%%% plot scalar vorticity field %%%%%
plot_Re60(omega,t_end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_omega_dt=omega_eq(omega,n,m,xi,theta,h,fac1,fac2,fac3,L,U)
%%%%% expand omega to include boundary points %%%%%
temp=zeros(n,m);
index1=2:n-1; index2=2:m-1;
temp(index1,index2)=reshape(omega,n-2,m-2);
omega=temp; clear temp;
%%%%% compute stream function %%%%%
b=zeros(n,m); %rhs of matrix equation for psi
b(index1,index2)=fac1(index1,index2).*omega(index1,index2);  
%inhomogeneous boundary conditions 
b(n,index2)=exp(xi(n))*sin(theta(index2)); %free stream
b=reshape(b,n*m,1);
psi=reshape(U\(L\b),n,m);
%%%%% compute derivatives of omega %%%%%
% boundary conditions
omega(1,:)=(1/(2*h^2))*(psi(3,:)-8*psi(2,:));
    omega(n,:)=0;
    omega(:,1)=omega(:,m-1);
    omega(:,m)=omega(:,2);
    %derivative
    d_omega_dt=...
       fac2(index1,index2).*((psi(index1+1,index2)-psi(index1-1,index2))...
         .*(omega(index1,index2+1)-omega(index1,index2-1))...
          -(psi(index1,index2+1)-psi(index1,index2-1))...
         .*(omega(index1+1,index2)-omega(index1-1,index2)))...
    +fac3(index1,index2).*(omega(index1+1,index2)+omega(index1-1,index2)...
          +omega(index1,index2+1)+omega(index1,index2-1)...
          -4*omega(index1,index2));
     d_omega_dt=d_omega_dt(:);
end
