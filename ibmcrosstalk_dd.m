function ibmcrosstalk_dd(choice, frame, ddnumber, omega_q1, omega_q2, J, omega_zz)
% choice = 2;
% frame = 3;
% ddnumber = 100;
sX = [0 1; 1 0];
sY = [0 -1i; 1i 0];
sZ = [1 0; 0 -1];
unit = speye(2);
natom = 2;

% omega_q1 = 5; %5GHz
% omega_q2 = 4.8; %4.8GHz
% J = 100; %100KHz
% omega_zz = 400; %400KHz


neval = 2^natom;
nevaltruc = neval;

sX_1 = kron(sX, speye(2^(natom-1)));
sX_2 = kron(speye(2^(natom-1)),sX);

sZ_1 = kron(sZ, speye(2^(natom-1)));
sZ_2 = kron(speye(2^(natom-1)),sZ);

sY_1 = kron(sY, speye(2^(natom-1)));
sY_2 = kron(speye(2^(natom-1)),sY);

sII = kron(speye(2^(natom-1)), speye(2^(natom-1)));

A_1 = sX_1;
A_2 = sX_2;

A = sZ_1+sZ_2;

%A= sY_1+sY_2;
sZsZ = kron(sZ,sZ);
plus = [1/sqrt(2); 1/sqrt(2)];

zero0 = [1; 0];
one1 = [0; 1];
%Temperature and coupling parameters
beta = (1/2.6)*1e-9;
%beta = (1/1.56)*1e-9;
beta = (1/1.57)*1e-9;
%wc = 8*pi*1e9;
wc = 1*2*pi*1e12;
%gsq2pi = (1.27)*1e-4*2*pi; %g^2*pi
% gsq2pi = (1)*1e-3*2*pi;
gsq2pi = (0.0025)*1e-3*2*pi;

%gsq2pi = (0)*1e-3*2*pi;
betainv = 2.6*1e9; %1/beta
%betainv = 1.56*1e9;
betainv = 1.57*1e9;
%Pure state initialization
plus = [1/sqrt(2); 1/sqrt(2)];
zero0 = [1; 0];
one1 = [0; 1];

perturbation = 0.000000;

Hs = -1e9.*(omega_q1)./2.*sZ_1 +...
    -1e9.*(omega_q2)./2.*sZ_2 +...
    1e3.*J.*(sZsZ);

Hs = Hs - Hs(1,1)*eye(nevaltruc)

num_op_trunc = diag([0:nevaltruc-1]);
modified_num_op_trunc = diag([0 1 1 2]);
wd = 1e9.*omega_q1;



if frame == 1
    wd = 1e9.*omega_q1 - 2*1e3.*J;
elseif frame == 2
    wd = 1e9.*omega_q1 + 2*1e3.*J;
elseif frame == 3
    wd = 1e9.*omega_q1;
end
    

H_rotated = Hs - (wd).*modified_num_op_trunc;

[V,D] = eig(full(Hs));
if ~issorted(diag(D))
    [V,D] = eig(full(Hs));
    [D,I] = sort(diag(D));
    D = diag(D);
    V = V(:, I);
    fprintf('sorted');
end



v = sparse(neval,nevaltruc);
e = sparse(1,nevaltruc);

for ii = 1:nevaltruc
    v(:,ii) = sparse(V(:,ii));
    e(ii) = sparse(D(ii,ii));
end

if choice == 1
    psi0 = kron(plus, zero0);
elseif choice == 2
    psi0 = kron(plus, one1);
elseif choice == 3
    psi0 = kron(plus, plus);
end
    
rho0 = psi0*psi0';
%rho = rho0(:);
rho = rho0;


tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%timestep between each ode executions
step = 1000;

tf = 500e-7;
dt_me = tf/step;
tstep_me = 0:dt_me:tf;

fidelitylist_me = zeros(1, numel(tstep_me));
firstexcitedpoplist_me = zeros(1, numel(tstep_me));
secondexcitedpoplist_me = zeros(1, numel(tstep_me));
thirdexcitedpoplist_me = zeros(1, numel(tstep_me));


pluspoplist_me = zeros(1, numel(tstep_me));
% fourthexcitedpoplist_me = zeros(1, numel(tstep_me));



ablist = zeros(1, numel(tstep_me));
hbar = 1;
psi = psi0;


A_sp = []
B_sp = []

count = 0;
for index = 1:numel(tstep_me)
    Hsd = H_rotated;

    rhom = reshape(rho,[neval,neval]);
    rhomcb = rhom;
    
    
    if ddnumber == 5
        if mod(index, 200) == 0
            index
            count = count +1;
            rhomcb = sX_2*rhomcb*sX_2';
        end 
    elseif ddnumber == 10
         if mod(index, 100) == 0
            index
            count = count +1;
            rhomcb = sX_2*rhomcb*sX_2';
         end 
    elseif ddnumber == 20
         if mod(index, 50) == 0
            index
            count = count +1;
            rhomcb = sX_2*rhomcb*sX_2';
         end  
    elseif ddnumber == 100
         if mod(index, 10) == 0
            index
            count = count +1;
            rhomcb = sX_2*rhomcb*sX_2';
         end  
    end
        
    U_rotated = expm(1i*(wd).*modified_num_op_trunc.*tstep_me(index));
    fidelity = rhomcb(1,1);
    fidelity = zero0'*ptrace(rhomcb,2,[2 2])*zero0;
    firstexcitedpop = rhomcb(2,2);
    secondexcitedpop = rhomcb(3,3);
    thirdexcitedpop = rhomcb(4,4);
    pluspop = plus'*ptrace(rhomcb,2,[2 2])*plus;

    fidelitylist_me(1, index) = fidelity;   
    firstexcitedpoplist_me(1, index) = firstexcitedpop;
    secondexcitedpoplist_me(1, index) = secondexcitedpop;
    thirdexcitedpoplist_me(1, index) = thirdexcitedpop;
    
    pluspoplist_me(1,index) = pluspop;

    A = U_rotated*A*U_rotated';
    alindblad  = @(t, rho)lindblad(t, rho, Hsd, natom, gsq2pi, beta, betainv, wc, A,v,e,nevaltruc);
    t    = [tstep_me(index), tstep_me(index) + dt_me];
    %options = odeset('RelTol',1e-3,'AbsTol',1e-6);
    options = odeset('RelTol',1e-3,'AbsTol',1e-6);
    %options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    %[~, RHO] = ode45(alindblad, t, rhomcb(:));
    [~, RHO] = ode45(alindblad, t, rhomcb(:),options);
    % Final value of rho is initial value for next step:
    rho = RHO(end, :); 
    %rho = v*reshape(rho,[nevaltruc,nevaltruc])*v';   %in computational basis
    
    rho = rho(:);
end

eptime = toc


function drhodt = lindblad(~, rho, Hsd, natom, gsq2pi, beta, betainv, wc,A,v,e,nevaltruc)
neval = 2^natom;
rhomcb = sparse(reshape(rho, [nevaltruc,nevaltruc]));

drhodt = -1i*(Hsd*rhomcb - rhomcb*Hsd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = bsxfun(@minus, e.', e);

[sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),0.1,'DataScale',1);%AbsTol
length(sortedOutput(sortedOutput>0));
count = 0;
for w = sortedOutput(sortedOutput>0)
    gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));
    if isnan(gamma) || isinf(gamma)
        gamma = gsq2pi*betainv;
    end
    [b, a] = ind2sub(size(X), find(abs(X - w)<= 0.1));% AbsTol
    count = count+length(b);

    Lpcomponents1 = sparse(nevaltruc,nevaltruc); 

    for s = 1:length(b)
      matrixelement1 = v(:,a(s))'*A*v(:,b(s));     %j<->a, i<->b in paper
      Lpcomponent1 = matrixelement1*sparse(a(s),b(s),1,nevaltruc,nevaltruc);

      Lpcomponents1 = Lpcomponents1 + Lpcomponent1;
    end

    Lncomponents1 = Lpcomponents1';
    
    Lp1 = sqrt(gamma)*Lpcomponents1;
    Ln1 = sqrt(gamma*exp(-beta*w))*Lncomponents1;
end

gamma0 = gsq2pi*betainv;
%[b0, a0] = ind2sub(size(X), find(X == 0));
[b0, a0] = ind2sub(size(X), find(abs(X - 0)<=0.1));%AbsTol

L01 = sparse(nevaltruc,nevaltruc);

for s = 1:length(b0)
    matrixelement1 = v(:,a0(s))'*A*v(:,b0(s));
    L0component1 = matrixelement1*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);

    L01 = L01 + L0component1;

end

L01 = sqrt(gamma0)*L01;

%%%%%%%%%Superoperator2%%%%%%%%%    
drhodt = drhodt + (L01*rhomcb*(L01)'-0.5*(L01)'*L01*rhomcb - 0.5*rhomcb*(L01)'*L01);


%%%%%%%%%Superoperator2%%%%%%%%%
% drhodt = L*rhovcb;
%%%%%%%%%Superoperator2%%%%%%%%%
%whos


drhodt = drhodt(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = ptrace(p, traceout, dims)
	% check arguments
	if any(traceout > length(dims)) || any(traceout < 0)
		error('Invalid subsystem in traceout')
	end
	if (length(dims) == 1 && mod(length(p)/dims,1) ~= 0) || length(p) ~= prod(dims)
		error('Size of state p inconsistent with dims');
    end

	% remove singleton dimensions
	traceout = setdiff(traceout,find(dims == 1));
	dims = dims(dims ~= 1);


	% calculate systems, dimensions, etc.
	n = length(dims);
	rdim = dims(end:-1:1);
	keep = 1:n;
	keep(traceout) = [];
	dimtrace = prod(dims(traceout));
	dimkeep = length(p)/dimtrace;


	if any(size(p) == 1)
		% state vector
		if size(p,1) == 1

		end

		perm = n+1-[keep(end:-1:1),traceout];
		x = reshape(permute(reshape(p,rdim),perm),[dimkeep,dimtrace]);
		x = x*x';

    else
		perm = n+1-[keep(end:-1:1),keep(end:-1:1)-n,traceout,traceout-n];
		x = reshape(permute(reshape(p,[rdim,rdim]),perm), [dimkeep,dimkeep,dimtrace^2]);
		x = sum(x(:,:,[1:dimtrace+1:dimtrace^2]),3);

	end



