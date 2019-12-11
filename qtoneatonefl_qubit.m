% (1) First, I think it would be good to reproduce the established fluctuator results using the trajectories approach.  
% What I mean by this is the following.  Take a single qubit in a z-pointing magnetic field.  
% The evolution is just the evolution via the Schrodinger equation.  Now let us include a single fluctuator just as Ka-Wa presented
% , meaning that the evolution is now stochastic 
% because the fluctuator flipping is a Poisson process.  Let us do/show the simulations that reproduce the analytical results.
function qtoneatonefl_qubit
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

% tf = 17/10;
tf = 70/10;

B = 0*[0;0;1];
tic
% Quantum trajectory by using waiting time distribution
dt = tf/100;

% fidelitylistmat = zeros(ntraj, numel(tstep));
% unpsinormlistmat = zeros(ntraj, numel(tstep));
% jlistmat = zeros(ntraj, 1000);

counter = 0;

ntraj = 1;

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
    %M0 = 1/sqrt(2)*[1;0;1];
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
    
    
    %Hs = -1e9.*A_sp1(1).*(2*pi/2).*(sX_1) + 1e9.*B_sp1(1).*(2*pi/2).*Hp;
    %Hs = -A_sp1(1).*(2*pi/2).*(sX_1) + B_sp1(1).*(2*pi/2).*Hp;
    
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
    %mp = mp./norm;
    
    
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
        %Hs = 1e9*(1/2).*sX.*((-1)+T(jj)/tf)+ 1e9*(-1/2).*sZ.*T(jj)/tf;
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
% plot(Ti,env,'-')
% text(0.5, 0.75, sprintf('Decay Rate = %.3f', pf(1)))
%legend('Mx','My','Mz','Location','best')

axis tight



% figure
% plot(10*Ti,-4/1*10*10*log(abs(mpmatrix_mean)),'-')
% hold on
% plot(10*Ti,-4/1*10*10*log(abs(Mxmatrix_mean)),'-')
% ylim([-8 0])



% length(T1)
% length(M1)
% 
% [T2,M2] = ode15s(@(t,M)precession(-B,M),tstep,M0);
% 
% length(T2)
% length(M2)
% 
% 
% T = [T1;T2(2:end)];
% M = [M1;M2(2:end,:)];
% length(T)
% length(M)

%{
for i = 1:ntraj
  figure
  %plot(Tmatrix{1},M(:,i),'-')
  plot(Tmatrix{i},Mxmatrix{i},'-')
  hold on
  plot(Tmatrix{i},Mymatrix{i},'-')
  hold on
  plot(Tmatrix{i},Mzmatrix{i},'-')
  hold on
  legend('Mx','My','Mz','Location','best')
  figure
  plot(Tmatrix{i},Bbmatrix{i},'-')
  axis tight
  figure
  plot(Tmatrix{i},abs(mpmatrix{i}),'-')
  axis tight
end

%}


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
  
  

%{
for n = 1:ntraj
%for n = ntrajmin:ntraj
psi = psi0;
unpsi = psi0;
%dp = zeros(1,(neval^2 - neval + 1));
%dp = zeros(1,8^2 - 10 + 1);   
% generate random number 
%r1 = rand([1 1]);
r = rand([1 2]);
%fidelitylist = [];
fidelitylist_qt = zeros(1, numel(tstep_qt));
v = zeros(neval,nevaltruc);
%v = sparse(neval,nevaltruc);
e = zeros(1,nevaltruc);
%e = sparse(1,nevaltruc);
%w = zeros(1,nevaltruc);
%w = sparse(1,nevaltruc);
%gamma = zeros(1,nevaltruc);
jlist_qt = zeros(1, 1000);
jorder = 1;

%FID 
for index = 1:numel(tstep_qt)
    %clear dp
    %Hs = -1e9.*A_sp1(tstep_qt(index)./tf).*(sX_1) + 1e9.*B_sp1(tstep_qt(index)./tf).*((-1).*((1/4).*sZ_1));
    Hs = -1e9.*(sX_1);
    [V,D] = eig(full(Hs));
    %[V,D] = eigs(Hs,18,'sr',opts);
    %[V,D] = eigs(Hs,18,'sa',opts);
    if ~issorted(diag(D))
        %[V,D] = eigs(Hs,18,'sa');
        [V,D] = eig(full(Hs));
        [D,I] = sort(diag(D));
        D = diag(D);
        V = V(:, I);
        sorted = 1
    end
    for ii = 1:nevaltruc
        v(:,ii) = sparse(V(:,ii));
        e(ii) = sparse(D(ii,ii));
    end
    Hsd = v'*Hs*v;
    psicb = v'*psi;
    unpsicb = v'*unpsi;
    norm(psicb);
    %Fidelity
    fidelity = psicb(1)*psicb(1)';
    %fidelity = (v(:,1)'* psi) * (psi' * v(:,1));
    %fidelitylist = cat(2, fidelitylist, fidelity);
    fidelitylist_qt(1, index) = fidelity; 

    X = bsxfun(@minus, e.', e);   %matrix of \omega_{ba}
    %[sortedOutput,~,~] = unique(reshape(X,[1 nevaltruc^2]));
    [sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),0.1,'DataScale',1); %AbsTol
    length(sortedOutput(sortedOutput>0));
    w_unique = length(sortedOutput);

    dp = zeros(1, w_unique*8);

    %gamma0 = 2.*g.^2.*pi.*beta.^(-1);
    gamma0 = gsq2pi*betainv;
    %[b0, a0] = ind2sub(size(X), find(X == 0));
    [b0, a0] = ind2sub(size(X), find(abs(X - 0)<=0.1));%AbsTol
    L01 = sparse(nevaltruc,nevaltruc);
    for s = 1:length(b0)
        matrixelement1 = v(:,a0(s))'*A_1*v(:,b0(s));      %j<->a, i<->b in paper
        L0component1 = matrixelement1*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
        L01 = L01 + L0component1;
    end
    L01 = sqrt(gamma0)*L01;
    dp(1) = (psicb'*L01')*(L01*psicb)*dt_qt;

    H_eff = Hsd - (1i*hbar/2)*(L01'*L01);
    pdx = 1+1;

    %numbz = length(X( X(:)>0 )) 
    count = 0;
    for w = sortedOutput(sortedOutput>0)
        %pdx = pdx + 1;
        %gamma = (2*pi*g^2*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));  
        gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));   
        if isnan(gamma) || isinf(gamma)
            %gamma = 2.*g.^2.*pi.*beta.^(-1);
            gamma = gsq2pi*betainv;
        end
        %[b, a] = ind2sub(size(X), find(X == w));
        [b, a] = ind2sub(size(X), find(abs(X - w)<= 0.1));% AbsTol
        count = count+length(b);
        Lpcomponents1 = sparse(nevaltruc,nevaltruc); 
        for s = 1:length(b)
            matrixelement1 = v(:,a(s))'*A_1*v(:,b(s));
            Lpcomponent1 = matrixelement1*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
          %L = L + sqrt(natom*gamma)*nfactor*Lcomponent;
            Lpcomponents1 = Lpcomponents1 + Lpcomponent1;
        end
        Lncomponents1 = Lpcomponents1';

        Lp1 = sqrt(gamma)*Lpcomponents1;

        Ln1 = sqrt(gamma*exp(-beta*w))*Lncomponents1;

        dp(pdx) = (psicb'*Lp1')*(Lp1*psicb)*dt_qt;

        dp(pdx+1) = (psicb'*Ln1')*(Ln1*psicb)*dt_qt;

        pdx = pdx + 2;
        H_eff = H_eff - (1i*hbar/2)*(Lp1'*Lp1) - (1i*hbar/2)*(Ln1'*Ln1);
    end
    count;
    pdx;
    length(dp);

    dp;
    %nonzeros(dp)
    dpj = sum(dp);
    %H_eff = sparse(H_eff);
    U_eff = expm(-1i*dt_qt*H_eff/hbar);

    dp0 = 1 - dpj;

    %unpsi_prev = unpsi;
    unpsi_prev = unpsicb;
    norm2_prev = norm(unpsi_prev)^2;
    %unpsi = U_eff*unpsi;     %Evolve until the norm of unpsi becomes r1
    unpsicb = U_eff*unpsicb;
    %%%%%%%

    %unpsi0 = unpsi;
    %[T2, UNPSI] = ode45(@(t,unpsi)dunpsi(t,unpsi,hbar,tau),[tstep_qt(index) tstep_qt(index)+dt_qt],unpsi0);
    %[T2, UNPSI] = ode45(@(t,unpsi)dunpsifast(unpsi,hbar,H_eff),[tstep_qt(index) tstep_qt(index)+dt_qt],unpsi0);
    %unpsi = UNPSI(end,:).';
    %toc
    %%%%%%
    %unpsinorm = sqrt(unpsi'*unpsi)

    %unpsinorm = norm(unpsi);
    unpsinorm = norm(unpsicb);

    unpsinormlist_qt(index) = unpsinorm;

    norm2_unpsi = unpsinorm^2;
    r1 = r(1);
    %if unpsi'*unpsi > r1
    if unpsicb'*unpsicb > r1
        %psi = unpsi/unpsinorm;
        psicb = unpsicb/unpsinorm;
        %Change back to comp.basis
        unpsi = v*unpsicb;
        psi = v*psicb;
        %Change back to comp.basis
    else % Rigorously should have implemented backtrack:
         % collapse has occured:
         % find collapse time to within specified tolerance
         % ------------------------------------------------
         % Rigorously should have implemented backtrack:
         % collapse has occured:
         % find collapse time to within specified tolerance
         % ------------------------------------------------
        t_prev = tstep_qt(index)
        t_final = tstep_qt(index) + dt_qt
        %r1;
        ii = 0;
        while ii < 5
            ii = ii + 1;
            t_guess = t_prev + (log(norm2_prev/r1)/log(norm2_prev/norm2_unpsi)) * (t_final - t_prev)
            %t_guess - t_prev
            unpsi_guess = expm(-1i*(t_guess - t_prev)*H_eff/hbar)*unpsi_prev;
            %norm2_guess = norm(unpsi_prev)^2
            norm2_guess = norm(unpsi_guess)^2;
            if abs(r1 - norm2_guess) < 0.001*r1  %error tolerance
                break
            elseif (norm2_guess < r1)
                t_final = t_guess;
                norm2_unpsi = norm2_guess;
            else
                t_prev = t_guess;
                unpsi_prev = unpsi_guess;
                norm2_prev = norm2_guess;
            end
        end
        %r2 = rand([1 1]);
        r2 = r(2);
        %condition = zeros(1, (nevaltruc^2 - nevaltruc + 1));
        condition = zeros(1, length(dp));
        cumsumlist = cumsum(dp);
        for m = 1:length(dp)
            condition(m) = (r2 < cumsumlist(m)/dpj);
        end
        k = find(condition,1,'first');
        Lk = lindbladsearch(k,natom,v,e,nevaltruc);
        psicb = Lk*psicb/norm(Lk*psicb);
        %jlist = [jlist k];
        jlist_qt(jorder) = k;
        jorder = jorder + 1;
        %%%%%%%%%%%%%%%%%%%%%%%%
%             if isempty(k) == true
%                 %psi = L{nLind}*psi/norm(L{nLind}*psi);
%                 Lend = lindbladsearch(pdx,natom,v,e);
%                 %psi = Lend*psi/norm(Lend*psi);
%                 psicb = Lend*psicb/norm(Lend*psicb);
%                 jlist = [jlist pdx]
%             else
%                 %psi = L{k}*psi/norm(L{k}*psi);
%                 Lk = lindbladsearch(k,natom,v,e);
%                 %psi = Lk*psi/norm(Lk*psi);
%                 psicb = Lk*psicb/norm(Lk*psicb);
%                 jlist = [jlist k]
%             end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %t_guess;

        if ~(t_guess >= tstep_qt(index) && t_guess <= tstep_qt(index) + dt_qt)
            t_guess = tstep_qt(index) + dt_qt
        end
%             if isnan(t_guess) || isinf(t_guess)
%                 t_guess = tstep_qt(index) + dt_qt;
%             end
        if t_guess > tf
            t_guess = tf;
        end
        t_guess;
        %Change back to comp.basis
        %unpsi = v*unpsicb;
        psi = v*psicb;
        %Change back to comp.basis

        Hs = -1e9.*A_sp1(t_guess./tf).*(sX_1) + 1e9.*B_sp1(t_guess./tf).*((-1).*((1/4).*sZ_1));
        [V,D] = eig(full(Hs));
        %[V,D] = eigs(Hs,18,'sa');
        if ~issorted(diag(D))
            %[V,D] = eigs(Hs,18,'sa');
            [V,D] = eig(full(Hs));
            [D,I] = sort(diag(D));
            D = diag(D);
            V = V(:, I);
            ok = 1
        end
        ii = 0;
        for ii = 1:nevaltruc
            v(:,ii) = sparse(V(:,ii));
            e(ii) = sparse(D(ii,ii));
        end

        Hsd = v'*Hs*v;
        psicb = v'*psi;

        X = bsxfun(@minus, e', e);
        %[sortedOutput,~,~] = unique(reshape(X,[1 nevaltruc^2]));
        [sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),0.1,'DataScale',1);%AbsTol
        w_unique = length(sortedOutput);
        dp = zeros(1, w_unique*8);

        %gamma0 = 2.*g.^2.*pi.*beta.^(-1);
        gamma0 = gsq2pi*betainv;
        %[b0, a0] = ind2sub(size(X), find(X == 0));
        [b0, a0] = ind2sub(size(X), find(abs(X - 0)<=0.1));%AbsTol
        L01 = sparse(nevaltruc,nevaltruc);
        for s = 1:length(b0)
            matrixelement1 = v(:,a0(s))'*A_1*v(:,b0(s));
            L0component1 = matrixelement1*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
            L01 = L01 + L0component1;
        end
        L01 = sqrt(gamma0)*L01;
        dp(1) = (psicb'*L01')*(L01*psicb)*dt_qt;

        H_eff = Hsd - (1i*hbar/2)*(L01'*L01);
        pdx = 1+1;
        for w = sortedOutput(sortedOutput>0)
            %pdx = pdx + 1;
            %gamma = (2*pi*g^2*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));   
            gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));    
            if isnan(gamma) || isinf(gamma)
                %gamma = 2.*g.^2.*pi.*beta.^(-1);
                gamma = gsq2pi*betainv;
            end
            %[b, a] = ind2sub(size(X), find(X == w));
            [b, a] = ind2sub(size(X), find(abs(X - w)<= 0.1));% AbsTol
            Lpcomponents1 = sparse(nevaltruc,nevaltruc);  
            for s = 1:length(b)
                matrixelement1 = v(:,a(s))'*A_1*v(:,b(s));

                Lpcomponent1 = matrixelement1*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
              %L = L + sqrt(natom*gamma)*nfactor*Lcomponent;
                Lpcomponents1 = Lpcomponents1 + Lpcomponent1;
            end
            Lncomponents1 = Lpcomponents1';

            Lp1 = sqrt(gamma)*Lpcomponents1;

            Ln1 = sqrt(gamma*exp(-beta*w))*Lncomponents1;

            dp(pdx) = (psicb'*Lp1')*(Lp1*psicb)*dt_qt;

            dp(pdx+1) = (psicb'*Ln1')*(Ln1*psicb)*dt_qt;

            pdx = pdx + 2;
            %dp(pdx) = (psi'*Lp')*(Lp*psi)*dt_qt;
            H_eff = H_eff - (1i*hbar/2)*(Lp1'*Lp1) - (1i*hbar/2)*(Ln1'*Ln1);
        end
        pdx;

        dp;
        dpj = sum(dp);

        U_eff = expm(-1i*((tstep_qt(index) + dt_qt)-t_guess)*H_eff/hbar);
        %unpsi = U_eff*psi;
        unpsicb = U_eff*psicb;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %unpsi0 = psi;
        %[~,UNPSI] = ode45(@(t,unpsi)dunpsifast(unpsi,hbar,H_eff),[t_guess tstep_qt(index)+dt_qt],unpsi0);
        %unpsi = UNPSI(end,:).';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %psi = unpsi/norm(unpsi);
        psicb = unpsicb/norm(unpsicb);


        %Change back to comp.basis
        unpsi = v*unpsicb;
        psi = v*psicb;
        %Change back to comp.basis   

        %r1 = rand([1 1]); % Draw a new random number    
        r = rand([1 2]); % Draw new random numbers   
    end
end
    
%fidelitylistall = cat(1, fidelitylistall, fidelitylist);
fidelitylistmat_qt(n, :) = fidelitylist_qt;
jlistmat_qt(n,:) = jlist_qt;
unpsinormlistmat_qt(n,:) = unpsinormlist_qt;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ntrajmin = n+1;  %If we load this checkpoint, we want to start on the next iteration
%     fprintf('Completed iteration %d \n',n);
%     if mod(n,5)==0
%         %save the state of the random number genertor
%          stream = RandStream.getGlobalStream;
%          savedState = stream.State;
%          tic
%          %save checkpoint
%          fprintf('Saving checkpoint\n');
%          save('checkpoint.mat');     
%          timing = toc;
%          fprintf('Checkpoint save took %f seconds\n',timing);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %counter = counter + 1
end

if ntraj == 1
    fidelitylistavg_qt = fidelitylistmat_qt;
else
    fidelitylistavg_qt = mean(fidelitylistmat_qt);
end   
eptime = toc


graphlist = [1 2];
txtlist = [1 2 3 5];

if ismember(1,graphlist)
    figure(1)
    h = plot(tstep_qt, fidelitylistavg_qt,'-');
    xlabel('$t$','Interpreter','latex')
    xlim([0 tf])
    ylabel('$\rho_{--}(t)$','Interpreter','latex')
    title(['tf: ' num2str(tf) ', Number of trajectories: ' num2str(ntraj)])
    print -dpdf qt8atoms_sparsemkron_fc_check_vectorized_omega_cb_new1
end
if ismember(2,graphlist)
    figure(2)
    h = plot(tstep_qt./tf, fidelitylistavg_qt,'-');
    xlabel('$s$','Interpreter','latex')
    ylabel('$\rho_{--}(s)$','Interpreter','latex')
    title(['tf: ' num2str(tf) ', Number of trajectories: ' num2str(ntraj)])
    print -dpdf qt8atoms_sparsemkron_fc_check_vectorized_omega_cb_new2
end
if ismember(3,graphlist)
    figure(3)
    h = plot(tstep_qt./tf, unpsinormlistmat_qt,'-');
    xlabel('$s$','Interpreter','latex')
    ylabel('unpsinorm','Interpreter','latex')
    title(['tf: ' num2str(tf) ', Number of trajectories: ' num2str(ntraj)])
    print -dpdf qt8atoms_sparsemkron_fc_check_vectorized_omega_cb_new3
end
if ismember(1,txtlist)
    txt1 = sprintf('fidelity_tf%d_ntraj%d.txt',tf,ntraj);
    fid1 = fopen(txt1,'w');
    fprintf(fid1,'%13d %8d\n',[tstep_qt;fidelitylistavg_qt]);
    fclose(fid1);
end    
if ismember(2,txtlist)
    fid2 = fopen('qt8atoms_fc_cb_new_eptime.txt','w');
    fprintf(fid2,'%d\n',eptime);
    fclose(fid2);
end
if ismember(3,txtlist)
    fid3 = fopen('qt8atoms_jlistmat_cb.txt','w');
    % fprintf(fid3,'%d\n',jlistmat);
    % fclose(fid3);
    [~, cols] = size(jlistmat_qt');
    x = repmat('%d\t',1,(cols-1));
    fprintf(fid3,[x,'%d\n'],jlistmat_qt);
    fclose(fid3);
end
if ismember(5,txtlist)
    fid5 = fopen('qt8atoms_unpsinormlistmat_cb.txt','w');
    [~, cols] = size(unpsinormlistmat_qt');
    x = repmat('%d\t',1,(cols-1));
    fprintf(fid5,[x,'%d\n'],unpsinormlistmat_qt);
    fclose(fid5);
end
delete(gcp('nocreate'))


% Subfunction for single-trajectory of deterministic evolution ODE 
% Fast: Constant Hamiltonian
function unpsidot = dunpsifast(unpsi,hbar,H_eff)
% No jump
unpsidot = -1i*hbar*H_eff*unpsi;




function L = lindbladsearch(k,natom,v,e,nevaltruc)
sZ = [1 0; 0 -1];
neval = 2^natom;
z = ceil(k/8);

g = sqrt((1.2/(2*pi))*1e-4);
beta = (1/2.6)*1e-9;
wc = 8*pi*1e9;
gsq2pi = (1.2)*1e-4;
gsq2pi = (1.27)*1e-4*2*pi;
betainv = 2.6*1e9;

sZ_1 = sZ;

A_1 = sZ_1;


X = bsxfun(@minus, e.', e);
%[sortedOutput,~,~] = unique(reshape(X,[1 nevaltruc^2]));
[sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),0.1,'DataScale',1);%AbsTol

% if z == 1
%    column = ceil(length(sortedOutput)/2);   
% elseif logical(mod(z,2))
%    column = length(sortedOutput) - (z-3)/2;
% else
%    column = z/2;
% end
if z == 1
   column = ceil(length(sortedOutput)/2);   
elseif logical(mod(z,2))
   column = ceil(length(sortedOutput)/2) - (z-1)/2;
else
   column = ceil(length(sortedOutput)/2) + z/2;
end



w = sortedOutput(column);
%gamma = (2*pi*g^2*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));
gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));
if isnan(gamma) || isinf(gamma)
    %gamma = 2.*g.^2.*pi.*beta.^(-1);
    gamma = gsq2pi*betainv;
end
%[b, a] = ind2sub(size(X), find(X == w));
[b, a] = ind2sub(size(X), find(abs(X - w)< 0.1));% AbsTol
L = zeros(nevaltruc,nevaltruc); 
if mod(k,1) == 1
else
    A = A_1;
end

for s = 1:length(b)
  matrixelement = v(:,a(s))'*A*v(:,b(s));
  Lcomponent = matrixelement*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
  L = L + Lcomponent;
end
L = sqrt(gamma)*L;


%}
