function qtoneatonefl_qubit(tf, ntraj)
% tf = 70/10;
% ntraj = 100;

% % Define Pauli matrices
%rng('shuffle')
sX = [0 1; 1 0];
sY = [0 -1i; 1i 0];
sZ = [1 0; 0 -1];
sX_1 = sX;
sZ_1 = sZ;
A_1 = sZ_1;
A = A_1;

Hp = sZ_1;
plus = [1/sqrt(2); 1/sqrt(2)];
plusy = [1/sqrt(2); 1i/sqrt(2)];

psi0 = sparse(plus);
rho0 = psi0*psi0';
hbar = 1;



B = 0*[0;0;1];
tic
% Quantum trajectory by using waiting time distribution
dt = tf/100;

counter = 0;

Tmatrix = cell(1,ntraj);
Mxmatrix = cell(1,ntraj);
Mymatrix = cell(1,ntraj);
Mzmatrix = cell(1,ntraj);
mpmatrix = cell(1,ntraj);
Bbmatrix = cell(1,ntraj);
psimatrix = cell(1,ntraj);
fidelitymatrix = cell(1,ntraj);
% b = 5.*[0;0;1];
dlm = dlmread('DW2000_parameters.txt');
%dlm = dlmread('DW1_parameters.txt');
slist = dlm(:,1).';
A_s = dlm(:,2).';
B_s = dlm(:,3).';
A_sp1 = @(s)interp1(slist,A_s,s);
B_sp1 = @(s)interp1(slist,B_s,s);

for n = 1:ntraj
    %b = -5.*[0;0;1];
    if mod(n,2) == 1
        b = 2.5*[0;0;1];
    else
        b = -2.5*[0;0;1];
    end
    psi0 = sparse(plus);
    psi = psi0;
    %gamma = 16;
    gamma = 2;
    T = 0;
    M0 = [1;0;0];
    M = M0.';    

    Bb = B + b;
    Q = -log(rand(1))/(gamma/2);
    tstep = [0 Q];
    %tstep = 0:dt:Q;
    if Q > tf
       tstep = [0 tf];
       %tstep = 0:dt:tf;
    end

    for it = 1:100  %it_max needs to be large for large gamma
       %[T1,M1] = ode45(@(t,M)precession(B+b,M),tstep,M0);
       options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    
       [T1,psi1] = ode23tb(@(t,psi)dpsi(b,psi,t,A_sp1,B_sp1,tf),tstep,psi0, options);
       norm(psi1(end,:))
       Mx1 = psi1(:,1).*conj(psi1(:,2)) + conj(psi1(:,1)).*psi1(:,2);
       My1 = 1i*psi1(:,1).*conj(psi1(:,2)) - 1i*conj(psi1(:,1)).*psi1(:,2);
       Mz1 = psi1(:,1).*conj(psi1(:,1)) - psi1(:,2).*conj(psi1(:,2));
       
       
       M1 = [Mx1 My1 Mz1];

       M0 = M1(end,:).';
       psi0 = psi1(end,:).';
       
       norm2 = M0(1)^2 + M0(2)^2 + M0(3)^2;
       for k = 2:length(T1)
           Bb = [Bb B + b];
       end
       b = -b;
       T = [T;T1(2:end)];
       M = [M;M1(2:end,:)];
       psi = [psi psi1(2:end,:).'];
       
       
       if T1(end) == tf
           break
       end
       
       Q = T1(end)-log(rand(1))/(gamma/2);
       tstep = [T1(end) Q];
       if Q > tf
           tstep = [T1(end) tf];
           %tstep = T1(end):dt:tf;
       end
    end

    norm3 = M(:,1).^2 + M(:,2).^2 + M(:,3).^2;
    norm4 = sqrt(norm3);
    mp = M(:,1)+1i*M(:,2);

    T(end);
    Tmatrix{n} = T;
    Mxmatrix{n} = M(:,1);
    Mymatrix{n} = M(:,2);
    Mzmatrix{n} = M(:,3);
    mpmatrix{n} = mp;
    Bbmatrix{n} = Bb(3,:);
    psimatrix{n} = psi;
    
    
    fidelitylist = zeros(1, length(T));
    for jj = 1:length(T) 
        Hs = -((A_sp1(T(jj)./tf))).*(2*pi/2).*(sX_1) + ((B_sp1(T(jj)./tf))).*(2*pi/2).*sZ_1;
        [V,D] = eig(Hs);
        [D,I] = sort(diag(D));
        D = diag(D);
        V = V(:, I);
        
        v0 = V(:,1);
        v1 = V(:,2);
        e0 = D(1,1);
        e0list(1, jj) = e0;
        e1 = D(2,2);
        e1list(1, jj) = e1;    
        v0 = sparse(V(:,1));
        
        psi(:,jj);
        fidelity = abs(v0'* psi(:,jj))^2;
        fidelitylist(1,jj) = fidelity;    
    end
    fidelitymatrix{n} = fidelitylist;
   
end

avglength = 0;
for n = 1:ntraj
    avglength = avglength + length(Tmatrix{n});
end
avglength = floor(avglength/ntraj)

% interpolation
Ti = linspace(0, tf, avglength);    
for n = 1:ntraj
    Mxmatrix{n} = interp1(Tmatrix{n}, Mxmatrix{n}, Ti, 'linear', 'extrap');
    Mymatrix{n} = interp1(Tmatrix{n}, Mymatrix{n}, Ti, 'linear', 'extrap');
    Mzmatrix{n} = interp1(Tmatrix{n}, Mzmatrix{n}, Ti, 'linear', 'extrap');
    mpmatrix{n} = interp1(Tmatrix{n}, mpmatrix{n}, Ti, 'linear', 'extrap');
    Bbmatrix{n} = interp1(Tmatrix{n}, Bbmatrix{n}, Ti, 'linear', 'extrap');
    fidelitymatrix{n} = interp1(Tmatrix{n}, fidelitymatrix{n}, Ti, 'linear', 'extrap');
end

Mxmatrix_mean = 0;
Mymatrix_mean = 0;
Mzmatrix_mean = 0;
mpmatrix_mean = 0;
Bbmatrix_mean = 0;
fidelitymatrix_mean = 0;
for n = 1:ntraj
    Mxmatrix_mean = Mxmatrix_mean + Mxmatrix{n};
    Mymatrix_mean = Mymatrix_mean + Mymatrix{n};
    Mzmatrix_mean = Mzmatrix_mean + Mzmatrix{n};
    mpmatrix_mean = mpmatrix_mean + mpmatrix{n};
    Bbmatrix_mean = Bbmatrix_mean + Bbmatrix{n};
    fidelitymatrix_mean =  fidelitymatrix_mean + fidelitymatrix{n};
end
Mxmatrix_mean = Mxmatrix_mean./ntraj;
Mymatrix_mean = Mymatrix_mean./ntraj;
Mzmatrix_mean = Mzmatrix_mean./ntraj;
mpmatrix_mean = mpmatrix_mean./ntraj;
Bbmatrix_mean = Bbmatrix_mean./ntraj;
fidelitymatrix_mean =  fidelitymatrix_mean./ntraj;

env = abs(hilbert(Mxmatrix_mean));            % Calculate Envelope
pf = polyfit(Ti, log(env), 1); 

toc

figure
%plot(Tmatrix{1},M(:,i),'-')
plot(Ti,Mxmatrix_mean,'-','LineWidth',2)
hold on
plot(Ti,Mymatrix_mean,'-','LineWidth',2)
hold on
plot(Ti,Mzmatrix_mean,'-','LineWidth',2)
hold on
% plot(Ti,env,'-')
% text(0.5, 0.75, sprintf('Decay Rate = %.3f', pf(1)))
legend('Mx','My','Mz','Location','best')
figure
plot(Ti,Bbmatrix_mean,'-')
axis tight
figure
plot(Ti,real(mpmatrix_mean),'-')
hold on
plot(Ti,sign(real(mpmatrix_mean)).*(real(mpmatrix_mean).^2),'-')
%plot(Ti,real(mpmatrix_mean).*abs(real(mpmatrix_mean)),'-')
axis tight


figure
%plot(Tmatrix{1},M(:,i),'-')
plot(Ti,fidelitymatrix_mean,'-','LineWidth',2)

axis tight


header1 = 'Ti';
header2 = 'Mx_mean';
header3 = 'My_mean';
header4 = 'Mz_mean';
fid=fopen('blochb03_testing.txt','w');
fprintf(fid, [ header1 ' ' header2 ' ' header3 ' ' header4 '\n']);
fprintf(fid, '%f %f %f %f\n', [Ti.' Mxmatrix_mean.' Mymatrix_mean.' Mzmatrix_mean.'].');
fclose(fid);



fid=fopen('fidb03.txt','w');
%fprintf(fid, [ header1 ' ' header2 '\n']);
fprintf(fid, '%f %f\n', [Ti.' Mxmatrix_mean.'].');
fclose(fid);


header1 = 'Ti';
header2 = 'df';
fid=fopen('dfb03.txt','w');
fprintf(fid, [ header1 ' ' header2 '\n']);
fprintf(fid, '%f %f\n', [Ti.' log(abs(mpmatrix_mean)).'].');
fclose(fid);


header1 = 'Ti';
header2 = 'gsp';
fid=fopen('gsp.txt','w');
fprintf(fid, [ header1 ' ' header2 '\n']);
fprintf(fid, '%f %f\n', [Ti.' fidelitymatrix_mean.'].');
fclose(fid);



function Mdot = precession(B,M)
   Mdot = zeros(3,1);
   %Mdot = cross(B,M)
   Mdot(1) = B(2)*M(3) - B(3)*M(2);
   Mdot(2) = B(3)*M(1) - B(1)*M(3);
   Mdot(3) = B(1)*M(2) - B(2)*M(1);

   
function psidot = dpsi(b,psi,t,A_sp1,B_sp1,tf)
sX = [0 1; 1 0];
sY = [0 -1i; 1i 0];
sZ = [1 0; 0 -1];
sX_1 = sX;
sZ_1 = sZ;
hbar = 1;
%Hs = -((A_sp1(t./tf))).*(2*pi/2).*(sX_1) + ((B_sp1(t./tf))).*(2*pi/2).*sZ_1;
Hs = -((A_sp1(t./tf))).*(2*pi/2).*(sX_1) + ((B_sp1(t./tf))).*(2*pi/2).*sZ_1;
H_eff = Hs + sZ.*(b(3)); %b(3) for stochastic b-field in z-axis
psidot = -1i*hbar*H_eff*psi;
  
  
