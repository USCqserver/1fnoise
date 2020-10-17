function qtoneatnfl_TRANSMON_new_slowonly_dd(E_J, E_C, tf, ntraj, nd, dec, bmean, bvariance)
% E_J, E_C and r from external parameters
%%E_J = get_E_J(i_c, constants);
% E_J = 13.4;
%%E_C = get_E_C(c_shunt, constants);
% E_C = 0.287;
% tf = 2845*5;
% ntraj = 300;
% nd = 10;
% dec = 2;
% bmean = 9.2e7./6.5e10;
% bvariance = 0;

sX = [0 1; 1 0];
sY = [0 -1i; 1i 0];
sZ = [1 0; 0 -1];

constants.H_PLANCK = 6.62607004e-34;
constants.PHI_0 = 2.067833831e-15;
constants.E_CHARGE = 1.6021766208e-19;
constants.BOLTZ = 20836612000.0;

hbar = 1;
i_c = 190;
c_shunt = 45;
alpha = 0.46;
d = 0.03;
amp = 0.326;

%E_J, E_C and r from external parameters
%E_J = get_E_J(i_c, constants);
%E_J = 13.4;
%E_C = get_E_C(c_shunt, constants);
%E_C = 0.287;
r = get_r(E_J, E_C);

%temperature
T = 20;
beta = convert_T_2_beta(T, constants);

%beta = (1/2.6)
%Temperature and coupling parameters
wc = 8*pi;
gsq2pi = 0.08*2*pi;
%betainv = 2.6;
betainv = 1/beta;


qmax = 10;
nmax = 2*qmax + 1;
n = diag(-qmax:qmax);
n2 = n^2;
d1 = diag(ones(1, nmax-1), -1);
d2 = diag(ones(1, nmax - 2), -2);


phi_x = @(s_tmp)(0.5)*pi;
phi_x = @(s_tmp)(0.04)*pi;
neval = nmax;
nevaltruc = 4;

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
Bbmatrix = cell(1,ntraj);
psimatrix = cell(1,ntraj);
fidelitymatrix = cell(1,ntraj);
firstexcitedpopmatrix = cell(1,ntraj);
secondexcitedpopmatrix = cell(1,ntraj);
thirdexcitedpopmatrix = cell(1,ntraj);
fourthexcitedpopmatrix = cell(1,ntraj);
pluspopmatrix = cell(1,ntraj);

t1listmatrix = cell(1,ntraj);
t2listmatrix = cell(1,ntraj);
% b = 5.*[0;0;1];




numsteps = 1000;
dt_qt = tf/numsteps;
tstep_qt = 0:dt_qt:tf;


for n = 1:ntraj
    b = bmean + sqrt(bvariance)*randn(1, nd*dec);
    x = 2*(rand(1, nd*dec, 1)<=.54)-1;   %0.08 equilibrium
    b = b.*x;
    ra = rand([ 1 2 ]); 
    Ip = get_Ip(d, phi_x(0), E_J, d2, constants);
    size(Ip);
    Ip = Ip*1e9;
    Hs = sparse(get_H(phi_x(0), n2, d, d1, E_C, E_J));

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

    plus = (v(:,1)+v(:,2))/sqrt(2);
    psi0 = plus;
    psi = psi0;

    
    K = linspace(log10(0.01), log10(0.1), nd*dec);   

    gammalist = 10.^(K);
    length(gammalist)
    gamma = sum(gammalist);    
        
    Bb = B + sum(b);
    T = 0;


    Q = -log( ra(1) ) / ( gamma/2 );
    tstep = [0 Q];
    %tstep = 0:dt:Q;
    if Q > tf
       tstep = [0 tf];
       %tstep = 0:dt:tf;
    end
    
    fidelitylist = zeros(1, length(tstep_qt));
    firstexcitedpoplist = zeros(1, length(tstep_qt));
    secondexcitedpoplist = zeros(1, length(tstep_qt));
    thirdexcitedpoplist = zeros(1, length(tstep_qt));
    fourthexcitedpoplist = zeros(1, length(tstep_qt));
    pluspoplist = zeros(1, length(tstep_qt)); 
    t1list = zeros(1, length(tstep_qt));

    Hs = sparse(get_H(phi_x(0), n2, d, d1, E_C, E_J));
    [V,D] = eig(full(Hs));
    [D,I] = sort(diag(D));
    D = diag(D);
    V = V(:, I);

    e0 = D(1,1);
    e1 = D(2,2); 
    v0 = sparse(V(:,1));
    v1 = sparse(V(:,2));
    v2 = sparse(V(:,3));
    v3 = sparse(V(:,4));
    v4 = sparse(V(:,5));
    v5 = sparse(V(:,6));
    
    wd = e1 - e0;
    num_op_trunc = diag([0:nevaltruc-1]);
    

    for index = 1:numel(tstep_qt)
        plus = (v0+v1)/sqrt(2);
        fidelity = abs(v0'* psi)^2;
        firstexcitedpop = abs(v1'* psi)^2;
        secondexcitedpop = abs(v2'* psi)^2;
        thirdexcitedpop = abs(v3'* psi)^2;
        fourthexcitedpop = abs(v4'* psi)^2;
        pluspop = abs(plus'* psi)^2;
        

        t1com = abs(v1'* psi)^2;
        fidelitylist(1,index) = fidelity; 
        firstexcitedpoplist(1,index)= firstexcitedpop;
        secondexcitedpoplist(1,index)= secondexcitedpop;
        thirdexcitedpoplist(1,index)= thirdexcitedpop;
        fourthexcitedpoplist(1,index)= fourthexcitedpop;
        pluspoplist(1,index) = pluspop;
        t1list(1,index) = t1com;        
        

        for ii = 1:nevaltruc
            v(:,ii) = sparse(V(:,ii));
            e(ii) = sparse(D(ii,ii));
        end
        
        Hsd = v'*Hs*v;
        psicb = v'*psi;   
        Ipd = v'*Ip*v;
        
        %
        AA = v0*v0' - v1*v1';
        AA = sparse(eye(nmax));
        AA(1,1) = 1;
        AA(2,2) = -1;
        %
        
        
        Hsd = Hsd - e0.*eye(nevaltruc);
        H_eff = Hsd + Ipd.*(sum(b)/2);
        
        
        U_rotated = expm(1i*wd.*num_op_trunc.*tstep_qt(index));
        H_eff2 = Hsd - wd.*num_op_trunc + U_rotated*Ipd*U_rotated'.*(sum(b)/2);
        U_eff = expm(-1i*dt_qt*H_eff2/hbar);
        psicb = U_eff*psicb;
        psi = v*psicb;
        
        if tstep_qt(index) > Q
           index;
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
           Q = tstep_qt(index) -log(ra(1))/(gamma/2);
           if Q > tf
               Q = tf
               %tstep = T1(end):dt:tf;
           end
        end
    end

    fidelitymatrix{n} = fidelitylist;
    firstexcitedpopmatrix{n} = firstexcitedpoplist;
    secondexcitedpopmatrix{n} = secondexcitedpoplist;
    thirdexcitedpopmatrix{n} = thirdexcitedpoplist;
    fourthexcitedpopmatrix{n} = fourthexcitedpoplist;
    
    pluspopmatrix{n} = pluspoplist;
    t1listmatrix{n} = t1list;
    
end


Bbmatrix_mean = 0;
fidelitymatrix_mean = 0;
firstexcitedpopmatrix_mean = 0;
secondexcitedpopmatrix_mean = 0;
thirdexcitedpopmatrix_mean = 0;
fourthexcitedpopmatrix_mean = 0;
pluspopmatrix_mean = 0;
t1listmatrix_mean = 0;
for n = 1:ntraj
    Bbmatrix_mean = Bbmatrix_mean + Bbmatrix{n};
    fidelitymatrix_mean =  fidelitymatrix_mean + fidelitymatrix{n};
    firstexcitedpopmatrix_mean = firstexcitedpopmatrix_mean + firstexcitedpopmatrix{n};
    secondexcitedpopmatrix_mean = secondexcitedpopmatrix_mean + secondexcitedpopmatrix{n};
    thirdexcitedpopmatrix_mean = thirdexcitedpopmatrix_mean + thirdexcitedpopmatrix{n};
    fourthexcitedpopmatrix_mean = fourthexcitedpopmatrix_mean + fourthexcitedpopmatrix{n};
    pluspopmatrix_mean = pluspopmatrix_mean + pluspopmatrix{n}; 
    t1listmatrix_mean = t1listmatrix_mean + t1listmatrix{n};
end

Bbmatrix_mean = Bbmatrix_mean./ntraj;
fidelitymatrix_mean =  fidelitymatrix_mean./ntraj;
firstexcitedpopmatrix_mean = firstexcitedpopmatrix_mean./ntraj;
secondexcitedpopmatrix_mean = secondexcitedpopmatrix_mean./ntraj;
thirdexcitedpopmatrix_mean = thirdexcitedpopmatrix_mean./ntraj;
fourthexcitedpopmatrix_mean = fourthexcitedpopmatrix_mean./ntraj;

pluspopmatrix_mean = pluspopmatrix_mean./ntraj;
t1listmatrix_mean = t1listmatrix_mean./ntraj;


toc
% exit(0)



% function Mdot = precession(B,M)
%    Mdot = zeros(3,1);
%    %Mdot = cross(B,M)
%    Mdot(1) = B(2)*M(3) - B(3)*M(2);
%    Mdot(2) = B(3)*M(1) - B(1)*M(3);
%    Mdot(3) = B(1)*M(2) - B(2)*M(1);

   
function psidot = dpsi(b,psi,t,tf)
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
E_J = 15.3;
E_J = 13.4;
E_C = get_E_C(c_shunt, constants);
E_C = 0.25;
E_C = 0.287;
r = get_r(E_J, E_C);

%temperature
T = 20;
beta = convert_T_2_beta(T, constants);

wc = 8*pi;
%gsq2pi = (1.2)*1e-4;
%gsq2pi = (1.27)*1e-4*2*pi;
gsq = 1e-4;
gsq2pi = 1e-4*2*pi;
gsq2pi = 0.008*2*pi;
%betainv = 2.6;
betainv = 1/beta;

%default qmax=60 => n, n2, d1, d2
qmax = 10;
nmax = 2*qmax + 1;
n = diag(-qmax:qmax);
n2 = n^2;
d1 = diag(ones(1, nmax-1), -1);
d2 = diag(ones(1, nmax - 2), -2);


phi_x = @(s_tmp)(0.5)*pi;


Ip = get_Ip(d, phi_x(0), E_J, d2, constants);
size(Ip);
Ip =Ip*1e9;
Hs = sparse(get_H(phi_x(0), n2, d, d1, E_C, E_J));
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


function Ip = get_Ip(d, phi_x, E_J, d1, constants)
cosphi = 0.5*(d1+d1');
sinphi = 0.5/1i*(d1-d1');
% Ip_hat = 2*alpha*E_J*(cos(phi_x/2)*sqrt(1+ d^2*tan(phi_x/2)^2)...
%     * -1/2/1i...
%     * (expm(1i * (phi_z - phi_0)) * d2 - expm(-1i * (phi_z - phi_0)) * d2'));
Ip_hat = E_J*(sin(phi_x/2)*cosphi - d*cos(phi_x/2)*sinphi);
Ip = Ip_hat*(constants.H_PLANCK)/(constants.PHI_0)*1e9*pi;


function H = get_H(phi_x, n2, d, d1, E_C,E_J)
phi_0 = atan(d.*tan(phi_x/2));
%cosphi = 0.5*(expm(1i*(-phi_0))*d1 + expm(-1i*(-phi_0))*d1');
cosphi = 0.5*(expm(1i*(-phi_0))*d1 + expm(-1i*(-phi_0))*d1');
%cosphi = 0.5*(expm(1i*(phi_0))*d1 + expm(-1i*(phi_0))*d1');
H = 4*E_C.*n2 - E_J*cos(phi_x/2) .*...
    sqrt(1 + (d.^2).*((tan(phi_x/2))^2)) .*cosphi;

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



  
  





