function onefdist
%uniform distribution of Y = log(X) in [a, b]
%means
%1/((b-a)*X) distribution of X in [exp(a), exp(b)]
gammamin = 1; %exp(a)
gammamax = 10; %exp(b)

gammalist = gammamin:1:gammamax;
p = @(gamma)1/gamma;
arrayfun(p,gammalist)
pgamma = arrayfun(p,gammalist);
%pgamma = pgamma./sum(pgamma)

%%%%%
la = log(1);
lb = log(10);
Y = exp(la + (lb-la)*rand(1,100000));


if 1 == 1
figure(1)
h = histogram(Y);
h.NumBins = 100;
hold on
h2 = plot(gammalist, pgamma, '-');
set(h2, 'color','[0 .5 0]', 'linewidth',1);

xlabel('$\gamma$','Interpreter','latex')
xlim([0 gammamax])
ylabel('$p(\gamma)$','Interpreter','latex')
title('Distribution function of relaxation rates $p(\gamma)$','Interpreter','latex')


figure(2)
h = histogram(log(Y));
h.NumBins = 100;
xlabel('$\ln(\gamma)$','Interpreter','latex')
xlim([log(0) log(gammamax)])
ylabel('frequency')
title('Histogram of $\ln(\gamma)$','Interpreter','latex')


figure(3)
h2 = plot(gammalist, pgamma, '-')
set(h2, 'color','[0 .5 0]', 'linewidth',1);

end

%https://stackoverflow.com/questions/19067014/simulation-of-custom-probability-distribution-in-matlab
% F = [0 0.3 0.4 0.8 1];
% f100000 = repmat(F,1000,1);
% r5 = repmat(rand(100000,1),1,5);
% diff(sum(f100000>r5));
% 
% largeNumber = 100;
% a=repmat( [0], 1, largeNumber*0.34 );
% e=repmat( [13], 1, largeNumber*0.11 );
% j = [a e];
% jshu = j(randperm(length(j)));
% result = jshu(1:10);


%https://math.stackexchange.com/questions/507795/simulation-of-custom-probability-distribution-in-matlab

nd = 100000;
Y = linspace(1,10,nd);
Yf = 1./Y;
R = randsample(Y,nd,true,Yf);
figure(5)
h = histogram(R);
h.BinEdges = 1:10;
h.NumBins = 100;
xlabel('$\gamma$','Interpreter','latex')
xlim([1 10])
ylabel('Occurrence')
title('Distribution function of relaxation rates $p(\gamma)$ according to stackoverflow','Interpreter','latex')


la = log(1);
lb = log(10);
Y = exp(la + (lb-la)*rand(1,100000));
figure(6)
h = histogram(Y)
h.BinEdges = 1:10;
h.NumBins = 100;
xlabel('$\gamma$','Interpreter','latex')
xlim([1 10])
ylabel('Occurrence')
title('Distribution function of relaxation rates $p(\gamma)$ according to stackoverflow','Interpreter','latex')

