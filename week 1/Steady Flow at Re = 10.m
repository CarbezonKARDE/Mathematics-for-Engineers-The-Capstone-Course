function [psi, omega] = flow_around_cylinder_steady
Re=10; 
%%%%% define the grid %%%%%
n=101; m=101; % number of grid points
N=n-1; M=m-1; % number of grid intervals
h=pi/M; % grid spacing based on theta variable
xi=(0:N)*h; theta=(0:M)*h; % xi and theta variables on the grid
%%%%% Initialize the flow fields %%%%%
psi=zeros(n,m);
omega=zeros(n,m);
psi(n,:)=exp(xi(n))*sin(theta); %free stream bc
%%%%% Set relax params, tol, extra variables %%%%%
r_psi=1.8; r_omega=0.9; %relaxation parameters
delta=1.e-08; % error tolerance
error=2*delta; % initialize error variable
%%%%% Add any additional variable definitions here %%%%%
f=zeros(n,m);
fac=h^2*exp(2*xi);
%%%%% Main SOR Loop %%%%%
while (error > delta)
    psi_old = psi; omega_old = omega;
    for i=2:n-1
        for j=2:m-1
            psi(i,j)=(1-r_psi)*psi(i,j)+(r_psi/4)*...
            (psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1)+fac(i)*omega(i,j));
        end
    end
    error_psi=max(abs(psi(:)-psi_old(:)));
    omega(1,:)=(1/(2*h^2))*(psi(3,:)-8*psi(2,:)); %Second order bc
    for i=2:n-1
        for j=2:m-1
                f(i,j)=...
                (psi(i+1,j)-psi(i-1,j))*(omega(i,j+1)-omega(i,j-1))...
               -(psi(i,j+1)-psi(i,j-1))*(omega(i+1,j)-omega(i-1,j));
            omega(i,j)=(1-r_omega)*omega(i,j)+(r_omega/4)*...
             (omega(i+1,j)+omega(i-1,j)+omega(i,j+1)+omega(i,j-1)...
             +(Re/8)*f(i,j));
        end
    end
    error_omega=max(abs(omega(:)-omega_old(:)));
    error=max(error_psi, error_omega);
end
plot_Re10(psi);
