%kawayip@usc.edu 
function onefbfield
tf = 1e-9;
B = 0;
dt = tf / 10000;
Tfixed = 0: dt: tf;
ntraj = 2000;

Mxmatrixsum = zeros(1, numel(Tfixed));
Mymatrixsum = zeros(1, numel(Tfixed));
Mzmatrixsum = zeros(1, numel(Tfixed));
%Mx + iMy
mpmatrixsum = zeros(1, numel(Tfixed));
Bbmatrixsum = zeros(1, numel(Tfixed));

nd = 1000;
dec = 12;
bmean = 0.5*9.2e7;
bvariance = (0.2*bmean)^2;

parfor n = 1 : ntraj
    b = bmean + sqrt(bvariance)*randn(1, nd*dec);
    x = 2*(rand(1, nd*dec, 1)<=.54)-1;    %0.08 equilibrium
    b = b.*x;
    r = rand([ 1 2 ]);
    
%%%%%%%%%%%%
%1/f noise distribution
%%%%%%%%%%%%
%     Y = linspace(2*pi*10^0,2*pi*10^12,nd*dec);
%     Yf = 1./Y;
%     Yf = Yf./sum(Yf);
%     gammalist = randsample(Y,nd*dec,true,Yf);
%     length(gammalist);
%     gamma = sum(gammalist);


%%%%%%%%%%%%
%sample fixed number of fluctuators per decade 
%%%%%%%%%%%
%     Y1 = linspace(2*pi*10^11,2*pi*10^12,nd);
%     Yf1 = 1./Y1;
%     Yf1 = Yf1./sum(Yf1);
%     gammalist1 = randsample(Y1,nd,true,Yf1);
%     Y2 = linspace(2*pi*10^10,2*pi*10^11,nd);
%     Yf2 = 1./Y2;
%     Yf2 = Yf2./sum(Yf2);
%     gammalist2 = randsample(Y2,nd,true,Yf2);
%     Y3 = linspace(2*pi*10^9,2*pi*10^10,nd);
%     Yf3 = 1./Y3;
%     Yf3 = Yf3./sum(Yf3);
%     gammalist3 = randsample(Y3,nd,true,Yf3);    
%     Y4 = linspace(2*pi*10^8,2*pi*10^9,nd);
%     Yf4 = 1./Y4;
%     Yf4 = Yf4./sum(Yf4);
%     gammalist4 = randsample(Y4,nd,true,Yf4); 
%     Y5 = linspace(2*pi*10^7,2*pi*10^8,nd);
%     Yf5 = 1./Y5;
%     Yf5 = Yf5./sum(Yf5);
%     gammalist5 = randsample(Y5,nd,true,Yf5); 
%     Y6 = linspace(2*pi*10^6,2*pi*10^7,nd);
%     Yf6 = 1./Y6;
%     Yf6 = Yf6./sum(Yf6);
%     gammalist6 = randsample(Y6,nd,true,Yf6); 
%     Y7 = linspace(2*pi*10^5,2*pi*10^6,nd);
%     Yf7 = 1./Y7;
%     Yf7 = Yf7./sum(Yf7);
%     gammalist7 = randsample(Y7,nd,true,Yf7); 
%     Y8 = linspace(2*pi*10^4,2*pi*10^5,nd);
%     Yf8 = 1./Y8;
%     Yf8 = Yf8./sum(Yf8);
%     gammalist8 = randsample(Y8,nd,true,Yf8);
%     Y9 = linspace(2*pi*10^3,2*pi*10^4,nd);
%     Yf9 = 1./Y9;
%     Yf9 = Yf9./sum(Yf9);
%     gammalist9 = randsample(Y9,nd,true,Yf9);    
%     Y10 = linspace(2*pi*10^2,2*pi*10^3,nd);
%     Yf10 = 1./Y10;
%     Yf10 = Yf10./sum(Yf10);
%     gammalist10 = randsample(Y10,nd,true,Yf10);        
%     Y11 = linspace(2*pi*10^1,2*pi*10^2,nd);
%     Yf11 = 1./Y11;
%     Yf11 = Yf11./sum(Yf11);
%     gammalist11 = randsample(Y11,nd,true,Yf11);    
%     Y12 = linspace(2*pi*10^0,2*pi*10^1,nd);
%     Yf12 = 1./Y12;
%     Yf12 = Yf12./sum(Yf12);
%     gammalist12 = randsample(Y12,nd,true,Yf12);        
%     
%     gammalist = [gammalist1 gammalist2(1:end) gammalist3(1:end) gammalist4(1:end) gammalist5(1:end)...
%         gammalist6(1:end) gammalist7(1:end) gammalist8(1:end) gammalist9(1:end) gammalist10(1:end)...
%         gammalist11(1:end) gammalist12(1:end)];
%     length(gammalist);
%     gamma = sum(gammalist);
%     
%%%%%%%%%%%%
%     Y1 = linspace(10^9,10^12,nd);
%     Yf1 = 1./Y1;
%     Yf1 = Yf1./sum(Yf1);
%     gammalist1 = randsample(Y1,3*nd,true,Yf1);
% 
%     Y5 = linspace(10^6,10^9,nd);
%     Yf5 = 1./Y5;
%     Yf5 = Yf5./sum(Yf5);
%     gammalist5 = randsample(Y5,3*nd,true,Yf5); 
%     
%     Y9 = linspace(10^3, 10^6,nd);
%     Yf9 = 1./Y9;
%     Yf9 = Yf9./sum(Yf9);
%     gammalist9 = randsample(Y9,3*nd,true,Yf9);    
%   
%     Y12 = linspace(10^0, 10^3,nd);
%     Yf12 = 1./Y12;
%     Yf12 = Yf12./sum(Yf12);
%     gammalist12 = randsample(Y12,3*nd,true,Yf12);        
%     
%     gammalist = [gammalist1 gammalist5(1:end) gammalist9(1:end) gammalist12(1:end)];
%     length(gammalist);
%     gamma = sum(gammalist);
%uniform fixed number of fluctuators per decade
%%%%%%%%%%%%
%agree above

    %range of frequency
    K = linspace(log(10^(0)), log(10^(12)), nd*dec);
    gammalist = exp(K);
    length(gammalist)
    gamma = sum(gammalist)
    
    Bb = B + sum(b);
    T = 0;
    M0 = [ 1; 0];
    M = M0;
    %waiting time
    Q = -log( r(1) ) / ( gamma/2 );
    tstep = [ 0 Q ];
    if Q > tf
       tstep = [0 tf];
    end
    for it = 1 : 1000000000000000
       [T1, M1] = myfun(B+sum(b), tstep, M0, dt);
       M0 = M1(:, end);
       for k = 2 : length(T1)
           Bb = [Bb B+sum(b)];
       end
       condition = zeros( 1, length( gammalist ) );
       cumsumlist = cumsum( gammalist );
       for m = 1 : length( gammalist )
           condition(m) = (r(2) < cumsumlist(m) / gamma);
       end
       k = find( condition, 1, 'first' );
       %flip 
       b(k) = -b(k);
       
       r = rand([1 2]); % Draw new random numbers   
       T = [T T1(2:end)];
       M = [M M1(:, 2:end)];
       if T1(end) == tf
           break
       end
       Q = T1(end) - log( r(1) ) / ( gamma/2 );
       tstep = [T1(end) Q];
       if Q > tf
           tstep = [T1(end) tf];
       end
    end
    %Mx + iMy
    mp = M(1, :) + 1i * M(2, :);
    Mxmatrixsum = Mxmatrixsum + interp1( T, M(1, :), Tfixed, 'linear', 'extrap' ); 
    Mymatrixsum = Mymatrixsum + interp1( T, M(2, :), Tfixed, 'linear', 'extrap' );
    mpmatrixsum = mpmatrixsum + interp1( T, mp, Tfixed, 'linear', 'extrap' );
    Bbmatrixsum = Bbmatrixsum + interp1( T, Bb(1, : ), Tfixed, 'nearest', 'extrap' );
end

Mxmatrix_mean2 = Mxmatrixsum./ntraj;
Mymatrix_mean2 = Mymatrixsum./ntraj;
Mzmatrix_mean2 = Mzmatrixsum./ntraj;
mpmatrix_mean2 = mpmatrixsum./ntraj;
Bbmatrix_mean2 = Bbmatrixsum./ntraj;


figure(5)
plot( Tfixed, Mxmatrix_mean2, '-')
hold on
plot( Tfixed, Mymatrix_mean2, '-')
hold on
plot( Tfixed, Mzmatrix_mean2, '-')
hold on
legend('Mx','My','Mz','Location','best')
title('Mx My Mz from matrix sum (memory)')
xlabel('$t$','Interpreter','latex')
ylabel('M component')
print -dpdf oneatnflfig5

figure(6)
plot( Tfixed, Bbmatrix_mean2, '-' )
legend('Mean magnetic field','Location','best')
title('Mean magnetic field from matrix sum')
xlabel('$t$','Interpreter','latex')
ylabel('$\left< b \right>$','Interpreter','latex')
axis tight
print -dpdf oneatnflfig6

figure(7)
plot(Tfixed, Mxmatrix_mean2,'-','Color','cyan')
hold on
plot(Tfixed, abs(mpmatrix_mean2),'-','Color',[0 0 0.5])
hold on
plot(Tfixed,sign(real(mpmatrix_mean2)).*(real(mpmatrix_mean2).^2),'--','Color',[0.7 0.7 0.7])
legend('Mxmatrix mean','abs(mpmatrix mean)','testing','Location','best')
xlabel('$t$','Interpreter','latex')
ylabel('$\left<m_+ \right>$','Interpreter','latex')
title('Time dependence of the FID amplitude from matrix sum.')
axis tight
print -dpdf oneatnflfig7


header1 = 'Tfixed';
header2 = 'Mx_mean';
header3 = 'My_mean';
header4 = 'Mz_mean';
fid=fopen('blochb03.txt','w');
fprintf(fid, [ header1 ' ' header2 ' ' header3 ' ' header4 '\n']);
%fprintf(fid, '%f %f %f %f\n', [Ti.' Mxmatrix_mean.' Mymatrix_mean.' Mzmatrix_mean.'].');
fprintf(fid, '%.10e %.10f %.10f %.10f\n', [Tfixed.' Mxmatrix_mean2.' Mymatrix_mean2.' Mzmatrix_mean2.'].');
fclose(fid);

fid=fopen('fidb03.txt','w');
fprintf(fid, '%.10e %.10f\n', [Tfixed.' Mxmatrix_mean2.'].');
fclose(fid);


header1 = 'Tfixed';
header2 = 'df';
fid=fopen('dfb03.txt','w');
fprintf(fid, [ header1 ' ' header2 '\n']);
fprintf(fid, '%.10e %.10f %.10f %.10f %.10f\n', [Tfixed.' log(abs(Mxmatrix_mean2)).' log(abs(mpmatrix_mean2)).' log(abs(real(mpmatrix_mean2))).' log(sign(real(mpmatrix_mean2)).*(real(mpmatrix_mean2).^2)).'].');
fclose(fid);


function [T1, M1] = myfun(b, tstep, M0, dt)
stepsize = dt;
nElements = ceil((tstep(end) - tstep(1))/(stepsize)) + 1;
T1 = linspace(tstep(1),tstep(end),nElements);
M1 = zeros(2, length(T1));
M1(:, 1) = M0;
for k = 2 : numel(T1)
    M1(1, k) = M0(1)*cos(b*(T1(k) - T1(1))) - M0(2)*sin(b*(T1(k) - T1(1)));
    M1(2, k) = M0(2)*cos(b*(T1(k) - T1(1))) + M0(1)*sin(b*(T1(k) - T1(1)));
end

function Mdot = precession(b , M)
   Mdot = zeros(2,1);
   Mdot(1) = -b*M(2);
   Mdot(2) = b*M(1);
  