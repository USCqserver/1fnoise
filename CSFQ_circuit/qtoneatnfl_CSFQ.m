function qtoneatnfl_CSFQ(tf, ntraj, nd, dec, bmean, bvariance)
% tf = 6;
% ntraj = 1000;
% nd = 25;
% dec = 2;
% bmean = 0.5*9.2e7; %need to divide by two to get to the qubit version
% bmean = 0.5*9.2e7./1e8;
% bvariance = (0.2*bmean)^2;

constants.H_PLANCK = 6.62607004e-34;
constants.PHI_0 = 2.067833831e-15;
constants.E_CHARGE = 1.6021766208e-19;
constants.BOLTZ = 20836612000.0;

hbar = 1;

%external parameters
i_c = 190;
c_shunt = 45;
alpha = 0.46;
d = 0.1;
amp = 0.326;

%E_J, E_C and r from external parameters
E_J = get_E_J(i_c, constants);
E_C = get_E_C(c_shunt, constants);
r = get_r(E_J, E_C);

%temperature
T = 20;
beta = convert_T_2_beta(T, constants);

%Temperature and coupling parameters
wc = 8*pi;
gsq = 1e-4;
gsq2pi = 1e-4*2*pi;
%betainv = 2.6;
betainv = 1/beta;

%default qmax=60 => n, n2, d1, d2
qmax = 60;
nmax = 2*qmax + 1;
n = diag(-qmax:qmax);
n2 = n^2;
d1 = diag(ones(1, nmax-1), -1);
d2 = diag(ones(1, nmax - 2), -2);

%phi_x, z_grid => phi_z
phi_x = @(s)(0.96 * s+ 1.04)*pi;
z_grid = linspace(-1.27,-1.20,100);
z = -1.21;
phi_z = @(s, z)(amp * s + z) * pi;
%neval and nevaltruc
neval = nmax;
nevaltruc = 5;

gaplist = [];
e1list = [];
e2list = [];
e3list = [];
e4list = [];
e5list = [];

% B = 0*[0;0;1];
B = 0;
tic
% Quantum trajectory by using waiting time distribution
counter = 0;

Tmatrix = cell(1,ntraj);
% Mxmatrix = cell(1,ntraj);
% Mymatrix = cell(1,ntraj);
% Mzmatrix = cell(1,ntraj);
% mpmatrix = cell(1,ntraj);
Bbmatrix = cell(1,ntraj);
psimatrix = cell(1,ntraj);
fidelitymatrix = cell(1,ntraj);
firstexcitedpopmatrix = cell(1,ntraj);
secondexcitedpopmatrix = cell(1,ntraj);
thirdexcitedpopmatrix = cell(1,ntraj);
fourthexcitedpopmatrix = cell(1,ntraj);

t1listmatrix = cell(1,ntraj);
t2listmatrix = cell(1,ntraj);
% b = 5.*[0;0;1];

for n = 1:ntraj
%     if mod(n,2) == 1
%         b = 2.5*[0;0;1];
%     else
%         b = -2.5*[0;0;1];
%     end
    b = bmean + sqrt(bvariance)*randn(1, nd*dec);
    x = 2*(rand(1, nd*dec, 1)<=.54)-1   %0.08 equilibrium
    %x = 2*(rand(1, nd*dec, 1)<=.5)-1   %0.0 equilibrium
    %x = ones(1, nd*dec)
    b = b.*x;
    ra = rand([ 1 2 ]); 
    
    Hs = sparse(get_H(phi_x(0), phi_z(0, z), r, n2, alpha, d, d1, d2, E_J));


    [V,D] = eig(full(Hs));
    if ~issorted(diag(D))
        [V,D] = eig(full(Hs));
        [D,I] = sort(diag(D));
        D = diag(D);
        V = V(:, I);
        fprintf('sorted');
    end

    e = zeros(1,nevaltruc);
    v = sparse(V(:,1:nevaltruc));
    for i = 1:nevaltruc
        e(i) = sparse(D(i,i));
    end
    Z = 0;
    for i = 1:nevaltruc
        Z = Z + exp(-beta*e(i));
    end


    v0 = v(:,1);
    %Pure state initialization
    psi0 = sparse(v0);

    psi = psi0;
    %range of frequency
    K = linspace(log10(0.01), log10(1), nd*dec);
%     gammalist = exp(K);
    gammalist = 10.^(K);
    length(gammalist)
    gamma = sum(gammalist);    
    
    
    Bb = B + sum(b);
    T = 0;
    Q = -log( ra(1) ) / ( gamma/2 )；
    tstep = [0 Q];
    if Q > tf
       tstep = [0 tf];
    end

    for it = 1:1000000  %it_max needs to be large for large gamma
       options = odeset('RelTol',1e-6,'AbsTol',1e-9);
       [T1,psi1] = ode23tb(@(t,psi)dpsi(sum(b),psi,t,tf,z),tstep,psi0, options);
       psi0 = psi1(end,:).'
       
       %norm2 = M0(1)^2 + M0(2)^2 + M0(3)^2;
       for k = 2:length(T1)
           Bb = [Bb B + sum(b)];
       end
       condition = zeros( 1, length( gammalist ) );
       cumsumlist = cumsum( gammalist );
       for m = 1 : length( gammalist )
           condition(m) = (ra(2) < cumsumlist(m) / gamma);
       end
       k = find( condition, 1, 'first' );
       %flip 
       b(k) = -b(k);
       %b = -b;
       ra = rand([1 2]); % Draw new random numbers   
       T = [T;T1(2:end)];
       %M = [M;M1(2:end,:)];
       psi = [psi psi1(2:end,:).'];
       if T1(end) == tf
           break
       end
       
       Q = T1(end)-log(ra(1))/(gamma/2);
       tstep = [T1(end) Q];
       if Q > tf
           tstep = [T1(end) tf];
           %tstep = T1(end):dt:tf;
       end
    end

    T(end);
    Tmatrix{n} = T;
    Bbmatrix{n} = Bb(1,:);
    psimatrix{n} = psi;

    fidelitylist = zeros(1, length(T));
    firstexcitedpoplist = zeros(1, length(T));
    secondexcitedpoplist = zeros(1, length(T));
    thirdexcitedpoplist = zeros(1, length(T));
    fourthexcitedpoplist = zeros(1, length(T));
    t1list = zeros(1, length(T));
    %t2list = zeros(1, length(T));
    ok = 1
    for jj = 1:length(T) 
        Hs = sparse(get_H(phi_x(T(jj)./tf), phi_z(T(jj)./tf, z), r, n2, alpha, d, d1, d2, E_J));
        [V,D] = eig(full(Hs));
        [D,I] = sort(diag(D));
        D = diag(D);
        V = V(:, I);
        e0 = D(1,1);
        e0list(1, jj) = e0;
        e1 = D(2,2);
        e1list(1, jj) = e1;    
        v0 = sparse(V(:,1));
        v1 = sparse(V(:,2));
        v2 = sparse(V(:,3));
        v3 = sparse(V(:,4));
        v4 = sparse(V(:,5));     
        v5 = sparse(V(:,6));
        fidelity = abs(v0'* psi(:,jj))^2;
        firstexcitedpop = abs(v1'* psi(:,jj))^2;
        secondexcitedpop = abs(v2'* psi(:,jj))^2;
        thirdexcitedpop = abs(v3'* psi(:,jj))^2;
        fourthexcitedpop = abs(v4'* psi(:,jj))^2;
        
        %v0'*v0 + v0'*v1 + v0'*v2 + v0'*v3 + v0'*v4
        %fidelity+firstexcitedpop+secondexcitedpop+thirdexcitedpop+fourthexcitedpop

        t1com = abs(v1'* psi(:,jj))^2;
        %t2com = abs(plus'* psi(:,jj))^2;
        fidelitylist(1,jj) = fidelity; 
        firstexcitedpoplist(1,jj)= firstexcitedpop;
        secondexcitedpoplist(1,jj)= secondexcitedpop;
        thirdexcitedpoplist(1,jj)= thirdexcitedpop;
        fourthexcitedpoplist(1,jj)= fourthexcitedpop;
        t1list(1,jj) = t1com;
        %t2list(1,jj) = t2com;
    end
    fidelitymatrix{n} = fidelitylist;
    firstexcitedpopmatrix{n} = firstexcitedpoplist;
    secondexcitedpopmatrix{n} = secondexcitedpoplist;
    thirdexcitedpopmatrix{n} = thirdexcitedpoplist;
    fourthexcitedpopmatrix{n} = fourthexcitedpoplist;
    t1listmatrix{n} = t1list;
    %t2listmatrix{n} = t2list;
    
end
Tmatrix{1};

avglength = 0;
for n = 1:ntraj
    avglength = avglength + length(Tmatrix{n});
end
avglength = floor(avglength/ntraj);

% interpolation
Ti = linspace(0, tf, avglength);    
for n = 1:ntraj
    Bbmatrix{n} = interp1(Tmatrix{n}, Bbmatrix{n}, Ti, 'linear', 'extrap');
    fidelitymatrix{n} = interp1(Tmatrix{n}, fidelitymatrix{n}, Ti, 'linear', 'extrap');
    firstexcitedpopmatrix{n} = interp1(Tmatrix{n}, firstexcitedpopmatrix{n}, Ti, 'linear', 'extrap');
    secondexcitedpopmatrix{n} = interp1(Tmatrix{n}, secondexcitedpopmatrix{n}, Ti, 'linear', 'extrap');
    thirdexcitedpopmatrix{n} = interp1(Tmatrix{n}, thirdexcitedpopmatrix{n}, Ti, 'linear', 'extrap');
    fourthexcitedpopmatrix{n} = interp1(Tmatrix{n}, fourthexcitedpopmatrix{n}, Ti, 'linear', 'extrap');
    t1listmatrix{n} = interp1(Tmatrix{n}, t1listmatrix{n}, Ti, 'linear', 'extrap');
    %t2listmatrix{n} = interp1(Tmatrix{n}, t2listmatrix{n}, Ti, 'linear', 'extrap');
end


Bbmatrix_mean = 0;
fidelitymatrix_mean = 0;
firstexcitedpopmatrix_mean = 0;
secondexcitedpopmatrix_mean = 0;
thirdexcitedpopmatrix_mean = 0;
fourthexcitedpopmatrix_mean = 0;
t1listmatrix_mean = 0;
%t2listmatrix_mean = 0;
for n = 1:ntraj
    Bbmatrix_mean = Bbmatrix_mean + Bbmatrix{n};
    fidelitymatrix_mean =  fidelitymatrix_mean + fidelitymatrix{n};
    firstexcitedpopmatrix_mean = firstexcitedpopmatrix_mean + firstexcitedpopmatrix{n};
    secondexcitedpopmatrix_mean = secondexcitedpopmatrix_mean + secondexcitedpopmatrix{n};
    thirdexcitedpopmatrix_mean = thirdexcitedpopmatrix_mean + thirdexcitedpopmatrix{n};
    fourthexcitedpopmatrix_mean = fourthexcitedpopmatrix_mean + fourthexcitedpopmatrix{n};
    t1listmatrix_mean = t1listmatrix_mean + t1listmatrix{n};
    %t2listmatrix_mean = t2listmatrix_mean + t2listmatrix{n};
end

Bbmatrix_mean = Bbmatrix_mean./ntraj;
fidelitymatrix_mean =  fidelitymatrix_mean./ntraj;
firstexcitedpopmatrix_mean = firstexcitedpopmatrix_mean./ntraj;
secondexcitedpopmatrix_mean = secondexcitedpopmatrix_mean./ntraj;
thirdexcitedpopmatrix_mean = thirdexcitedpopmatrix_mean./ntraj;
fourthexcitedpopmatrix_mean = fourthexcitedpopmatrix_mean./ntraj;
t1listmatrix_mean = t1listmatrix_mean./ntraj;
%t2listmatrix_mean = t2listmatrix_mean./ntraj;

% env = abs(hilbert(Mxmatrix_mean));            % Calculate Envelope
% pf = polyfit(Ti, log(env), 1); 

toc

figure(1)
plot(Ti, fidelitymatrix_mean,'-b','LineWidth',2);
hold on
plot(Ti, firstexcitedpopmatrix_mean,'-','color', [0 0.5 0], 'LineWidth',2);
hold on
plot(Ti, secondexcitedpopmatrix_mean,'-','color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
hold on
plot(Ti, thirdexcitedpopmatrix_mean,'-','color',[0.9290, 0.6940, 0.1250], 'LineWidth',2);
hold on
plot(Ti, fourthexcitedpopmatrix_mean,'-', 'color',[0.50, 0.1850, 0.5580], 'LineWidth',2);
hold on 
plot(Ti, fidelitymatrix_mean+firstexcitedpopmatrix_mean+secondexcitedpopmatrix_mean+thirdexcitedpopmatrix_mean+fourthexcitedpopmatrix_mean,'-k', 'LineWidth',2);
xlabel('$t$ (ns)','Interpreter','latex','FontSize',25)
ylabel('Eigenstate population','FontSize',25)
ylim([0 1])
set(gca,'FontSize',20)
legend('0thpop','1stpop','2ndpop','3rdpop','4thpop', 'location', 'west')
title(['$T =$ ' num2str(20), 'mK,    $\eta = 1e^{-4}$', ',   ntraj= ' num2str(ntraj), ',   $\phi_z(0)/\pi =$ ' num2str(z)],'Interpreter','latex','FontSize', 17)
print -dpdf scurvezgridqtfidelity_z-121_1f


figure(8)
h2 =  plot(Ti*1e9,fidelitymatrix_mean,'-','LineWidth',2);
set(h2(1), 'color','blue');
ax = ancestor(h2, 'axes');
xrule = ax.XAxis;
xrule.FontSize = 18;
yrule = ax.YAxis;
yrule.FontSize = 18;
xlabel('time (ns)','FontSize',23)
ylabel('Instantaneous GS Population','FontSize',23)
annotation('textbox',[.75 .6 .3 .3],'String','$T\rightarrow \infty$','Interpreter', 'latex','FitBoxToText','on','FontSize',25, 'EdgeColor', 'none');
title('Instantaneous Ground State Population vs time', 'Interpreter', 'latex','FontSize',15)
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
ax.TickLength = [0.02 0.035]
ax.YAxis.MinorTickValues = -0.1:0.1:1;
print -dpdf gspplot






figure(9)
h3 =  plot(Ti*1e9,t1listmatrix_mean,'-','LineWidth',2);
set(h3(1), 'color','[0 0.5 0]');
ax = ancestor(h3, 'axes');
xrule = ax.XAxis;
xrule.FontSize = 18;
yrule = ax.YAxis;
yrule.FontSize = 18;
xlabel('time (ns)','FontSize',23)
ylabel('Off-diagonal elements', 'FontSize',23)
title({'(Sum of) Off-diagonal elements (abs value) of';'density matrix in computational basis'},'Interpreter', 'latex','FontSize',15)
annotation('textbox',[.75 .6 .3 .3],'String','$T\rightarrow \infty$','Interpreter', 'latex','FitBoxToText','on','FontSize',25, 'EdgeColor', 'none');
set(gca,'XMinorTick','on');
set(gca,'YMinorTick','on');
ax.TickLength = [0.02 0.035]
ax.YAxis.MinorTickValues = -0.1:0.1:1;
print -dpdf t1plot


header1 = 'Ti';
header2 = 'gsp';
fid=fopen('gsp.txt','w');
fprintf(fid, [ header1 ' ' header2 '\n']);
fprintf(fid, '%f %f\n', [Ti.' fidelitymatrix_mean.'].');
fclose(fid);



header1 = 'Ti';
header2 = 't1mean';
fid=fopen('t1mean.txt','w');
fprintf(fid, [ header1 ' ' header2 '\n']);
fprintf(fid, '%f %f\n', [Ti.' t1listmatrix_mean.'].');
fclose(fid);

% 
% exit(0)


function psidot = dpsi(b,psi,t,tf,z)
sX = [0 1; 1 0];
sZ = [1 0; 0 -1];


constants.H_PLANCK = 6.62607004e-34;
constants.PHI_0 = 2.067833831e-15;
constants.E_CHARGE = 1.6021766208e-19;
constants.BOLTZ = 20836612000.0;

hbar = 1;
%hbar = 2*pi;

%external parameters
i_c = 190;
c_shunt = 45;
alpha = 0.46;
d = 0.1;
amp = 0.326;

%E_J, E_C and r from external parameters
E_J = get_E_J(i_c, constants);
E_C = get_E_C(c_shunt, constants);
r = get_r(E_J, E_C);

%temperature
T = 20;
beta = convert_T_2_beta(T, constants);

wc = 8*pi;
gsq = 1e-4;
gsq2pi = 1e-4*2*pi;
%betainv = 2.6;
betainv = 1/beta;


%default qmax=60 => n, n2, d1, d2
qmax = 60;
nmax = 2*qmax + 1;
n = diag(-qmax:qmax);
n2 = n^2;
d1 = diag(ones(1, nmax-1), -1);
d2 = diag(ones(1, nmax - 2), -2);



phi_x = @(s)(0.96 * s+ 1.04)*pi;
%z_grid = z_grid.*pi
%z_grid = linspace(-1.25,-1.18,100);
phi_z = @(s, z)(amp * s + z) * pi;
%neval and nevaltruc

Ip = get_Ip(d, phi_x(t./tf), phi_z(t./tf, z), alpha, E_J, d2, constants);
size(Ip);
Hs = sparse(get_H(phi_x(t./tf), phi_z(t./tf, z), r, n2, alpha, d, d1, d2, E_J));
size(Hs);
H_eff = Hs + Ip.*(b/2); %b(3) for stochastic b-field in z-axis

psidot = -1i*hbar*H_eff*psi;
  
  
function inverse = convert_T_2_beta(T, constants)
inverse = 1e12 / (2 * pi * T * constants.BOLTZ);

function E_J =  get_E_J(i_c, constants)
E_J = constants.PHI_0 * (i_c * 1e-9) / constants.H_PLANCK / 1e9;

function E_C = get_E_C(c_shunt, constants)
E_C = (2 * pi) * constants.E_CHARGE.^2 / (c_shunt * 1e-15) / constants.H_PLANCK / 1e9;

function r = get_r(E_J, E_C)
r = sqrt(E_C / (2 * E_J));


function phi = get_phi(nmax, qmax)
phi_hat = zeros(nmax, nmax); 
for q = -qmax : qmax 
    for q_p = -qmax : q - 1
        temp = (qmax * sin(2*pi*(qmax+1)*(q-q_p)/(2*qmax+1)) ...
            - (qmax+1)*sin(2*pi*qmax*(q-q_p)/(2*qmax+1)))...
            / (1i*((2*qmax+1)^2)*(1-cos(2*pi*(q-q_p)/(2*qmax+1))));
        phi_hat(q, q_p) = temp;
    end
end
phi_hat = phi_hat - phi_hat';
phi = phi_hat*2*pi;

function Ip = get_Ip(d, phi_x, phi_z, alpha, E_J, d2, constants)
phi_0 = atan(d*tan(phi_x/2));
Ip_hat = 2*alpha*E_J*(cos(phi_x/2)*sqrt(1+ d^2*tan(phi_x/2)^2)...
    * -1/2/1i...
    * (expm(1i * (phi_z - phi_0)) * d2 - expm(-1i * (phi_z - phi_0)) * d2'));
Ip = Ip_hat*(constants.H_PLANCK)/(constants.PHI_0)*1e18;


function H = get_H(phi_x, phi_z, r, n2, alpha, d, d1, d2, E_J)
phi_0 = atan(d.*tan(phi_x/2));

H = (r^2).*n2 - 2*0.5*(d1 + d1') + 2*alpha*cos(phi_x/2) .*...
    sqrt(1 + (d.^2).*((tan(phi_x/2))^2)) .*0.5*(expm(1i*(phi_z-phi_0))*d2 + expm(-1i*(phi_z-phi_0))*d2');

H = E_J.* H;

function [e10, abssin2phi_10] = get_params(phi_x, phi_z, r, n2, alpha, d2)
sin2phi = 1/(2i)*(d2 - d2')
Hs = get_H(phi_x, phi_z, r, n2, alpha);
[V,D] = eig(full(H));
if ~issorted(diag(D))
    [V,D] = eig(full(Hs));
    [D,I] = sort(diag(D));
    D = diag(D);
    V = V(:, I);
    fprintf('sorted');
end
[D,I] = sort(diag(D));
v = sparse(neval,nevaltruc);
e = sparse(1,nevaltruc);

for ii = 1:nevaltruc
    v(:,ii) = sparse(V(:,ii));
    e(ii) = sparse(D(ii,ii));
end
Hsd = v'*Hs*v;
v0 = v(:,1);
neval = 1;
rho = [1 1;1 1];
rhom = reshape(rho,[neval,neval]);
%rhomcb = v'*rhom*v;
rhomcb = sparse(v'*rhom*v);
e10 = e(1) - e(0);
sin2phi_10 = v(:,1)'*sin2phi*v(:,1);
abssin2phi_10 = abs(sin2phi_10);


