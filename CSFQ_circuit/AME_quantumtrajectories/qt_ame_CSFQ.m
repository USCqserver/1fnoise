function qt_ame_CSFQ(tf, ntraj)
% tf = 60;
% ntraj = 1000;

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
T = 20
beta = convert_T_2_beta(T, constants)

%beta = (1/2.6)
wc = 8*pi;

gsq = 1e-4;
gsq2pi = 1e-4*2*pi;
%betainv = 2.6;
betainv = 1/beta;

qmax = 60;
nmax = 2*qmax + 1;
n = diag(-qmax:qmax);
n2 = n^2;
d1 = diag(ones(1, nmax-1), -1);
d2 = diag(ones(1, nmax - 2), -2);


phi_x = @(s)(0.96 * s+ 1.04)*pi;
z_grid = linspace(-1.27,-1.20,100);

z = -1.25;
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

Ip = get_Ip(d, phi_x(0), phi_z(0, z), alpha, E_J, d2, constants);
size(Ip);
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
psi0 = v0;
rho = psi0*psi0';


dt_qt = tf/100000;
tstep_qt = 0:dt_qt:tf;
fidelitylistmat_qt = zeros(ntraj, numel(tstep_qt));
firstexcitedpoplistmat_qt = zeros(1, numel(tstep_qt));
secondexcitedpoplistmat_qt = zeros(1, numel(tstep_qt));
thirdexcitedpoplistmat_qt = zeros(1, numel(tstep_qt));
fourthexcitedpoplistmat_qt = zeros(1, numel(tstep_qt));
unpsinormlistmat_qt = zeros(ntraj, numel(tstep_qt));
jlistmat_qt = zeros(ntraj, numel(tstep_qt));


tic
for n = 1:ntraj
    psi = psi0;
    unpsi = psi0;
    ra = rand([1 2]);
    %fidelitylist = [];
    fidelitylist_qt = zeros(1, numel(tstep_qt));
    firstexcitedpoplist_qt = zeros(1, numel(tstep_qt));
    secondexcitedpoplist_qt = zeros(1, numel(tstep_qt));
    thirdexcitedpoplist_qt = zeros(1, numel(tstep_qt));
    fourthexcitedpoplist_qt = zeros(1, numel(tstep_qt));
    v = zeros(neval,nevaltruc);
    e = zeros(1,nevaltruc);
    jlist_qt = zeros(1, numel(tstep_qt));
    jorder = 1;
    unpsinormlist_qt = zeros(1, numel(tstep_qt));     
    
    
    for index = 1:numel(tstep_qt)
        Ip = get_Ip(d, phi_x(tstep_qt(index)./tf), phi_z(tstep_qt(index)./tf, z), alpha, E_J, d2, constants);
        A_1 = Ip;
        Hs = sparse(get_H(phi_x(tstep_qt(index)./tf), phi_z(tstep_qt(index)./tf, z), r, n2, alpha, d, d1, d2, E_J));
        [V,D] = eig(full(Hs));
        if ~issorted(diag(D))
            [V,D] = eig(full(Hs));
            [D,I] = sort(diag(D));
            D = diag(D);
            V = V(:, I);
            fprintf('sorted');
        end

        for ii = 1:nevaltruc
            v(:,ii) = sparse(V(:,ii));
            e(ii) = sparse(D(ii,ii));
        end
        Hsd = v'*Hs*v;
        psicb = v'*psi;
        unpsicb = v'*unpsi;
        norm(psicb);
        
        fidelity = psicb(1)*psicb(1)';
        firstexcitedprob = psicb(2)*psicb(2)';
        secondexcitedprob = psicb(3)*psicb(3)';
        thirdexcitedprob = psicb(4)*psicb(4)';
        fourthexcitedprob = psicb(5)*psicb(5)';
        
        fidelitylist_qt(1, index) = fidelity; 
        firstexcitedpoplist_qt(1, index) = firstexcitedprob;
        secondexcitedpoplist_qt(1, index) = secondexcitedprob;
        thirdexcitedpoplist_qt(1, index) = thirdexcitedprob;
        fourthexcitedpoplist_qt(1, index) = fourthexcitedprob;
        
        X = bsxfun(@minus, e.', e);
        [sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),(0.1)/1e9,'DataScale',1);%AbsTol
        length(sortedOutput(sortedOutput>0));
        w_unique = length(sortedOutput);
        dp = zeros(1, w_unique*1);

        gamma0 = gsq2pi*betainv;
        [b0, a0] = ind2sub(size(X), find(abs(X - 0)<=(0.1)/1e9));%AbsTol
        L01 = sparse(nevaltruc,nevaltruc);

        for s = 1:length(b0)
            matrixelement1 = v(:,a0(s))'*A_1*v(:,b0(s));
            %L0component1 = v'*(matrixelement1*v(:,i0(s))*v(:,j0(s))')*v
            L0component1 = matrixelement1*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
            L01 = L01 + L0component1;
        end
        L01 = sqrt(gamma0)*L01;
        dp(1) = (psicb'*L01')*(L01*psicb)*dt_qt;
        H_eff = Hsd - (1i*hbar/2)*(L01'*L01);
        pdx = 1+1;
        
        count = 0;
        for w = sortedOutput(sortedOutput>0)
            %gamma = (2*pi*g^2*w*exp(-abs(w)/wc))/(1 - exp(-beta*w)); 
            gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));
            if isnan(gamma) || isinf(gamma)
                %gamma = 2.*g.^2.*pi.*beta.^(-1);
                gamma = gsq2pi*betainv;
            end
            [b, a] = ind2sub(size(X), find(abs(X - w)<= (0.1)/1e9));% AbsTol
            count = count+length(b);
            Lpcomponents1 = sparse(nevaltruc,nevaltruc);

            for s = 1:length(b)
              %matrixelement = v(:,i(s))'*A*v(:,j(s));
              matrixelement1 = v(:,a(s))'*A_1*v(:,b(s));     %j<->a, i<->b in paper

              Lpcomponent1 = matrixelement1*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
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

        unpsinorm = norm(unpsicb);
        unpsinormlist_qt(index) = unpsinorm;        
        norm2_unpsi = unpsinorm^2;
        r1 = ra(1);
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
            t_prev = tstep_qt(index);
            t_final = tstep_qt(index) + dt_qt;
            %r1;
            ii = 0;
            while ii < 5
                ii = ii + 1;
                t_guess = t_prev + (log(norm2_prev/r1)/log(norm2_prev/norm2_unpsi)) * (t_final - t_prev);
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
            r2 = ra(2);
            %condition = zeros(1, (nevaltruc^2 - nevaltruc + 1));
            condition = zeros(1, length(dp));
            cumsumlist = cumsum(dp);
            for m = 1:length(dp)
                condition(m) = (r2 < cumsumlist(m)/dpj);
            end
            k = find(condition,1,'first');
            Lk = lindbladsearch(k,A_1,v,e,nevaltruc,beta);
            psicb = Lk*psicb/norm(Lk*psicb);
            %jlist = [jlist k];
            jlist_qt(jorder) = k;
            jorder = jorder + 1;


            if ~(t_guess >= tstep_qt(index) && t_guess <= tstep_qt(index) + dt_qt)
                t_guess = tstep_qt(index) + dt_qt;
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
                       
            Ip = get_Ip(d, phi_x(t_guess./tf), phi_z(t_guess./tf, z), alpha, E_J, d2, constants);
            A_1 = Ip;
            Hs = sparse(get_H(phi_x(t_guess./tf), phi_z(t_guess./tf, z), r, n2, alpha, d, d1, d2, E_J));
            [V,D] = eig(full(Hs));
            if ~issorted(diag(D))
                [V,D] = eig(full(Hs));
                [D,I] = sort(diag(D));
                D = diag(D);
                V = V(:, I);
                fprintf('sorted');
            end
            ii = 0;
            for ii = 1:nevaltruc
                v(:,ii) = sparse(V(:,ii));
                e(ii) = sparse(D(ii,ii));
            end
            Hsd = v'*Hs*v;
            psicb = v'*psi;
           

            X = bsxfun(@minus, e.', e);
            [sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),(0.1)/1e9,'DataScale',1);%AbsTol
            length(sortedOutput(sortedOutput>0));
            w_unique = length(sortedOutput);
            dp = zeros(1, w_unique*1);

            gamma0 = gsq2pi*betainv;
            [b0, a0] = ind2sub(size(X), find(abs(X - 0)<=(0.1)/1e9));%AbsTol
            L01 = sparse(nevaltruc,nevaltruc);

            for s = 1:length(b0)
                matrixelement1 = v(:,a0(s))'*A_1*v(:,b0(s));
                %L0component1 = v'*(matrixelement1*v(:,i0(s))*v(:,j0(s))')*v
                L0component1 = matrixelement1*sparse(a0(s),b0(s),1,nevaltruc,nevaltruc);
                L01 = L01 + L0component1;
            end
            L01 = sqrt(gamma0)*L01;
            dp(1) = (psicb'*L01')*(L01*psicb)*dt_qt;
            H_eff = Hsd - (1i*hbar/2)*(L01'*L01);
            pdx = 1+1;
            for w = sortedOutput(sortedOutput>0)
                %gamma = (2*pi*g^2*w*exp(-abs(w)/wc))/(1 - exp(-beta*w)); 
                gamma = (gsq2pi*w*exp(-abs(w)/wc))/(1 - exp(-beta*w));
                if isnan(gamma) || isinf(gamma)
                    %gamma = 2.*g.^2.*pi.*beta.^(-1);
                    gamma = gsq2pi*betainv;
                end
                [b, a] = ind2sub(size(X), find(abs(X - w)<= (0.1)/1e9));% AbsTol
                count = count+length(b);
                Lpcomponents1 = sparse(nevaltruc,nevaltruc);

                for s = 1:length(b)
                  matrixelement1 = v(:,a(s))'*A_1*v(:,b(s));     %j<->a, i<->b in paper
                  Lpcomponent1 = matrixelement1*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
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
            pdx;
            
            dp;
            dpj = sum(dp);
            U_eff = expm(-1i*((tstep_qt(index) + dt_qt)-t_guess)*H_eff/hbar);
            unpsicb = U_eff*psicb;

            psicb = unpsicb/norm(unpsicb);
            
            
            %Change back to comp.basis
            unpsi = v*unpsicb;
            psi = v*psicb;
            %Change back to comp.basis   
            
            %r1 = rand([1 1]); % Draw a new random number    
            ra = rand([1 2]); % Draw new random numbers   
        end  
    end
    %fidelitylistall = cat(1, fidelitylistall, fidelitylist);
    fidelitylistmat_qt(n, :) = fidelitylist_qt;
    firstexcitedpoplistmat_qt(n, :) = firstexcitedpoplist_qt;
    secondexcitedpoplistmat_qt(n, :) = secondexcitedpoplist_qt;
    thirdexcitedpoplistmat_qt(n, :) = thirdexcitedpoplist_qt;
    fourthexcitedpoplistmat_qt(n, :) = fourthexcitedpoplist_qt;

    jlistmat_qt(n,:) = jlist_qt;
    unpsinormlistmat_qt(n,:) = unpsinormlist_qt;
end    
    
if ntraj == 1
    fidelitylistavg_qt = fidelitylistmat_qt;
    firstexcitedpoplistavg_qt = firstexcitedpoplistmat_qt;
    secondexcitedpoplistavg_qt = secondexcitedpoplistmat_qt;
    thirdexcitedpoplistavg_qt = thirdexcitedpoplistmat_qt;
    fourthexcitedpoplistavg_qt = fourthexcitedpoplistmat_qt;
else
    fidelitylistavg_qt = mean(fidelitylistmat_qt);
    firstexcitedpoplistavg_qt = mean(firstexcitedpoplistmat_qt);
    secondexcitedpoplistavg_qt = mean(secondexcitedpoplistmat_qt);
    thirdexcitedpoplistavg_qt = mean(thirdexcitedpoplistmat_qt);
    fourthexcitedpoplistavg_qt = mean(fourthexcitedpoplistmat_qt);
end   
eptime = toc

figure(1)
plot(tstep_qt, fidelitylistavg_qt,'-b','LineWidth',2);
hold on
plot(tstep_qt, firstexcitedpoplistavg_qt,'-','color', [0 0.5 0], 'LineWidth',2);
hold on
plot(tstep_qt, secondexcitedpoplistavg_qt,'-','color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
hold on
plot(tstep_qt, thirdexcitedpoplistavg_qt,'-','color',[0.9290, 0.6940, 0.1250], 'LineWidth',2);
hold on
plot(tstep_qt, fourthexcitedpoplistavg_qt,'-', 'color',[0.50, 0.1850, 0.5580], 'LineWidth',2);

xlabel('$t$ (ns)','Interpreter','latex','FontSize',25)
ylabel('Eigenstate population','FontSize',25)
ylim([0 1])
set(gca,'FontSize',20)
legend('0thpop','1stpop','2ndpop','3rdpop','4thpop', 'location', 'west')
title(['$T =$ ' num2str(T), 'mK,    $\eta = 1e^{-4}$', ',   ntraj= ' num2str(ntraj), ',   $\phi_z(0)/\pi =$ ' num2str(z)],'Interpreter','latex','FontSize', 17)
print -dpdf scurvezgridqtfidelity_z-121



txt1 = sprintf('scurvezgridqtfidelity_T%d_g2%.6f_z%f_tf%d.txt',T,gsq,z,tf);
fid1 = fopen(txt1,'w');
fprintf(fid1,'%13d %8d %8d %8d %8d %8d\n',[tstep_qt;fidelitylistavg_qt;firstexcitedpoplistavg_qt;secondexcitedpoplistavg_qt;thirdexcitedpoplistavg_qt;fourthexcitedpoplistavg_qt]);
fclose(fid1);

fid2 = fopen('qtscurve_fc_cb_new_eptime.txt','w');
fprintf(fid2,'%d\n',eptime);
fclose(fid2);

fid3 = fopen('qtscurve_jlistmat_cb.txt','w');
% fprintf(fid3,'%d\n',jlistmat);
% fclose(fid3);
[~, cols] = size(jlistmat_qt');
x = repmat('%d\t',1,(cols-1));
fprintf(fid3,[x,'%d\n'],jlistmat_qt);
fclose(fid3);

fid5 = fopen('qtscurve_unpsinormlistmat_cb.txt','w');
[~, cols] = size(unpsinormlistmat_qt');
x = repmat('%d\t',1,(cols-1));
fprintf(fid5,[x,'%d\n'],unpsinormlistmat_qt);
fclose(fid5);

delete(pool)



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

(r^2).*n2;
- (d1 + d1'); 
2*alpha*cos(phi_x/2);
sqrt(1 + (d.^2).*((tan(phi_x/2))^2)).*0.5*(expm(1i*(phi_z-phi_0))*d2 + expm(-1i*(phi_z-phi_0))*d2');

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


function L = lindbladsearch(k,A_1,v,e,nevaltruc,beta)
z = ceil(k/1);

T = 20;
wc = 8*pi;
gsq = 1e-4;
gsq2pi = 1e-4*2*pi;
betainv = 1/beta;
   
X = bsxfun(@minus, e.', e);
%[sortedOutput,~,~] = unique(reshape(X,[1 nevaltruc^2]));
[sortedOutput,~,~] = uniquetol(full(reshape(X,[1 nevaltruc^2])),(0.1)/1e9,'DataScale',1);%AbsTol

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
[b, a] = ind2sub(size(X), find(abs(X - w)< (0.1)/1e9));% AbsTol
L = zeros(nevaltruc,nevaltruc); 
% if mod(k,1) == 1
%     A = A_1;
% end
A = A_1;

for s = 1:length(b)
  matrixelement = v(:,a(s))'*A*v(:,b(s));
  Lcomponent = matrixelement*sparse(a(s),b(s),1,nevaltruc,nevaltruc);
  L = L + Lcomponent;
end
L = sqrt(gamma)*L;








