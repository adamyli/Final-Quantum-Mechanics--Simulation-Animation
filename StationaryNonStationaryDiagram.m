%Stationary state plot%
L=10;
n=1;
x = linspace(0,L,1000); %1000 = number of points in between
m=1
Psi_n = sqrt(2/L)* sin((pi/L)*n*x)
E_n =((n*pi)^2 / (2*m*L));
Psi_sq = abs(Psi_n).^2

%figure
%for t=0:0.05:1000
    %Psi_Re= Psi_n*cos((E_n)*t);
    %Psi_Im = Psi_n*sin(-(E_n)*t);
    %plot(x,Psi_Re ,'b', x,Psi_Im, 'r',x,abs(Psi_Re+1i*Psi_Im).^2,'g')
    %axis([0,L, -1,1])
    %pause(0.01)
%end

%Non stationary state plot%
%We have coeffcients C_n to find now
%C_n = Integral of [(Psi_n^*)(Psi_n)] with bounds 0 to L
%Solve this definite integral with trapizoidal rule
dx = x(2)-x(1)
%Create a vector for our dx to be multipled into
once = dx*ones([1,1000])
once(1) = 0.5
once(1000) = 0.5

%Taking a non-stationary function and then normalizing it to find the
%normalization constant A
Psi_non = cos(x)
A = 1/sqrt((Psi_non.*once)*Psi_non'); %Here it is multipled by its complex con
Psi_non = A.*Psi_non;


n_matrix = 1:2;
Psi_basis = sqrt(2/L)* sin((pi/L)*n_matrix'*x);

C_n = (Psi_basis)*(Psi_non.*once)';

E_nmatrix =(((n_matrix).^2*(pi)^2) /(2*m*L));

figure
for t=0:0.05:1000
    E_Re = cos((E_nmatrix)*t);
    E_Im = sin(-(E_nmatrix)*t);
    PsiNon_Re = (C_n'.*E_Re)*Psi_basis;
    PsiNon_Im = (C_n'.*E_Im)*Psi_basis;
    plot(x,PsiNon_Re,'b',x,PsiNon_Im,'r',x,abs(PsiNon_Re+1i*PsiNon_Im).^2, 'g')
    axis([0,L, -1,1])
    pause(0.1)
end



