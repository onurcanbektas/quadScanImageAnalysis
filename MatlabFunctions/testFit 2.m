% clear
% addpath('/Volumes/Junction/matlab/Onurcan/MatlabFunctions')
% B = [0.155314, 0.170532, 0.186198, 0.201417, 0.21709, 0.232301, 0.247967, 0.26363228, 0.27885044, 0.29451619, 0.31018195]';
% k = B ./ 0.43835;
% load data29.txt
% data = data29(:,2);
% %data = data .* 10.^(-3);
% l = 0.244;
% L = 1.88;
% %phi = (k.^(1./2)) .* l;
% syms MQ_sym(kvar) MQD_sym(kvar) sigma(s11, s12, s22)
% MQ_sym(kvar) = [cos(kvar.^(1./2)), sin(kvar.^(1./2)) ./ (kvar.^(1./2)); (-1).* (kvar.^(1./2) .* sin(kvar.^(1./2))), cos( l .* kvar.^(1/2))];
% MDrift = [1, L; 0, 1];
% MQD_sym(kvar) = MDrift * MQ_sym(kvar)

% %sigma(s11, s12, s22) = [s11, s12; s12, s22];
% %syms result_sym(s11, s12, s22, k)
% %result_sym(s11, s12, s22, k) = MQD_sym(k) * sigma(s11, s12, s22) * MQD_sym(k)';
% %result_sym[1]
% guess = [2.9;4.3;6.2];
% guess = guess(:);
% %model = @(s) (cos(conj(k.^(1./2))) - (47.*conj(k.^(1./2)).*sin(conj(k.^(1./2))))./25).*(s(2).*((47.*cos((61.*k.^(1./2))./250))./25 + sin(k.^(1./2))./k.^(1./2)) + s(1).*(cos(k.^(1./2)) - (47.*k.^(1./2).*sin(k.^(1./2)))./25)) + ((47.*cos((61.*conj(k.^(1./2)))./250))./25 + sin(conj(k.^(1./2)))./conj(k.^(1./2))).*(s(3).*((47.*cos((61.*k.^(1./2))./250))./25 + sin(k.^(1./2))./k.^(1/2)) + s(2).*(cos(k.^(1./2)) - (47.*k.^(1./2).*sin(k.^(1./2)))./25)) - data;

% model = @(s) (L.^2 .* l.^2 .* s(1)) .* k.^2 + 2.*(L .* l .* s(1) + L.^2 .* l .* s(2)).*k + (s(1) + 2.*L.*s(2) + L.^2 .* s(2)) - data;
% [soln ,ssq,cnt, res, XY] = LMFnlsq(model,guess,'Display',1,'XTol', 1e-12)

% asd_temp = matlabFunction(MQD_sym);
% asd = @(t) asd_temp(t);
% asd(0.4952267) * [soln(1), soln(2);soln(2), soln(3)]

% figure(3)
% plot(k, data, 'x')
% hold on
% plot(k, model(soln) + data)
% hold off

clear
addpath('/Volumes/Junction/matlab/Onurcan/MatlabFunctions')
B = [0.155314, 0.170532, 0.186198, 0.201417, 0.21709, 0.232301, 0.247967, 0.26363228, 0.27885044, 0.29451619, 0.31018195]';
k = B ./ 0.043835;
l = 0.244;%% .* 10.^3;
b = (1) .* k .* l;
%load data29.txt
load data29y.txt
data = data29y(:,2);
data = data .^2;
L = 1.88;%% .* 10^3;
k = B ./ 0.43835;

p = polyfit(b,data,3)     % fit to a 1st degree polynomial: linefit
xfit = polyval(p,b);     % evaluate that polynomial for those 'time's 
plot(b,xfit,'-b');

sigma11 = sqrt(p(1) ./ (L .^2)) 
sigma12 = ((p(2) - (2 .* L .* sigma11.^2)) ./ (2.* L.^2))
sigma22 = sqrt(( p(3) - (sigma11.^2) - (2.* L .* sigma12)) ./ (L.^2))
