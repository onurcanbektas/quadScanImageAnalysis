%% @author: Onurcan Bektaş
%%          onurcan.bektas@metu.edu.tr
clear

lib_dir = strcat(pwd, '/matlabFunctions');
output_dir = strcat(pwd, '/output.29072017');
addpath(lib_dir)
addpath(output_dir)

filename = '29072017Y.txt';
%% (+1) for defocusing plane, (-1) for focusing plane
sign_focus_or_defocus = -1;

beam_check_ON = 1;
%% Needed only if sign_focus_or_defocus is equal to 1
rowofTheParameterToCheck = 6;

fileID = fopen(filename);
C = textscan(fileID,'%s %f %f %d %d %d %d %d %d %d');
data = C(1,2);
data = cell2mat(data);
data = data * (10^(-3));
data = data.^2;
B = C(3);
B = cell2mat(B);
k = B ./ (0.043835);
l = 0.244; %radius 0,055m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (+1) for defocusing plane, (-1) for focusing plane
b = (sign_focus_or_defocus) .* k .* l;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 1.6;

if b(1) > 0
    w = b(1):0.001:(b(length(b)-1) + 0.05);
elseif b(1) < 0
    w = (b(length(b)-1) - 0.05):0.001:b(1);
end

p = polyfit(b,data,2)
xfit = polyval(p,w);  
plot(w,xfit,'r');
hold on
plot(b, data, '*')
title('The measurement of beam size versus 1/f');
legend('fit','data','Location', 'north');
xlabel(' \pm 1/f  (in 1/m) ');
ylabel('Beam size sigma^2 (in m^2)');
%text(-1.6, 7.5*10^(-5), 'y = 10^{-4} * (0.2835*x^2 + 0.4427*x + 0.6641)')
hold off

sigma11 = sqrt(p(1) ./ (d .^2))
sigma12 = ((p(2) - (2 .* d .* sigma11.^2)) ./ (2.* d.^2))
sigma22 = sqrt(( p(3) - (sigma11.^2) - (2.* d .* sigma12)) ./ (d.^2))

emittance1 = sqrt(sigma11.^2 .* sigma22.^2 - sigma12.^2)

phi = acos(emittance1 ./ (sigma11 .* sigma22));
kor = sin(phi)

syms thin(kvar)
thin(kvar) = [1 + kvar * l * d, d; kvar*l, 1];
temp_thin = matlabFunction(thin);
thin = @(kvar) temp_thin(kvar);

k_temp = (sign_focus_or_defocus) .* k(rowofTheParameterToCheck);
thin(k_temp);
result_thin = thin(k_temp) * [sigma11.^2, sigma12; sigma12, sigma22.^2] * (thin(k_temp)');

Endsigma11 = result_thin(1,1).^(1./2)
Endsigma12 = result_thin(1,2)
Endsigma22 = result_thin(2,2).^(1./2)

emittance2 = sqrt( Endsigma11.^2 * Endsigma22.^2 - Endsigma12.^2)