clear
lib_dir = strcat(pwd, '/matlabFunctions')
addpath(pwd, '/matlabFunctions')
filename = "0308Y.txt";
%% (+1) for defocusing plane, (-1) for focusing plane
sign_focus_or_defocus = 1;
beam_check_ON = 1;
%% Needed only if sign_focus_or_defocus is equal to 1
rowofTheParameterToCheck = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% END OF LINE OF EDIT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen(filename);
C = textscan(fileID,'%s %f %f %d %d %d %d %d %d %d');
data = C(1,2);
data = cell2mat(data);
data = data * (10^(-3));
data = data.^2;
B = C(3);
B = cell2mat(B);
k = B ./ (0.043835);
l = 0.244;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (+1) for defocusing plane, (-1) for focusing plane
b = (sign_focus_or_defocus) .* k .* l;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 1.88;

if b(1) > 0
    w = b(1):0.001:(b(length(b) -1 ) + 0.2);
elseif b(1) < 0
    w = (b(length(b)-1 ) - 0.2):0.001:b(1);
end
p = polyfit(b,data,2)
xfit = polyval(p,w);  
plot(w,xfit,'r');
hold on
plot(b, data, '*')
hold off

sigma11 = sqrt(p(1) ./ (d .^2))
sigma12 = ((p(2) - (2 .* d .* sigma11.^2)) ./ (2.* d.^2))
sigma22 = sqrt(( p(3) - (sigma11.^2) - (2.* d .* sigma12)) ./ (d.^2))

emittance1 = sqrt(sigma11.^2 .* sigma22.^2 - sigma12.^2)


if(beam_check_ON == 1)

syms MQ_sym(kvar) MQD_sym(kvar) sigma(s11, s12, s22)
MQ_sym(kvar) = [cos((kvar.^(1./2)) .* l),  sin((kvar.^(1./2)) .* l) ./ (kvar.^(1./2)); (-1).* ((kvar.^(1./2)) .* sin((kvar.^(1./2)).* l)), cos((kvar.^(1./2)) .* l)];
MDrift = [1, d; 0, 1];
MQD_sym(kvar) = MDrift * MQ_sym(kvar);
tempMQD = matlabFunction(MQD_sym);
MQD  = @(kvar) tempMQD(kvar);

k_temp = k(rowofTheParameterToCheck);
result = MQD(k_temp) * [sigma11.^2, sigma12; sigma12, sigma22.^2] * MQD(k_temp)';

Endsigma11 = result(1,1).^(1./2)
Endsigma12 = result(1,2)
Endsigma22 = result(2,2).^(1./2)

emittance2 = sqrt(Endsigma11.^2 .* Endsigma22.^2 - Endsigma12.^2)
 
end
