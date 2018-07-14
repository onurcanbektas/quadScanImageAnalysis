function [GoddnessofFit, A,mean, sigma] = newtonsMethodForGaussianFitEdit(data, index)

%%index = [1, 4, 7]';
%%data = [2 5 4]';
guess = [240, 200, 20]';
cut_edit1 = 5;
cut_edit2 = 15;
cut_edit3 = 35;
cut_edit4 = 55;
data = [data(cut_edit1:cut_edit2); data(cut_edit3:cut_edit4)];
index = [index(cut_edit1:cut_edit2);index(cut_edit3:cut_edit4)];

%% x's are the parameters that are to be fitted with the data, and t's are the data itself
%% A = x(1); mean = x(2); sigma = x(3); offset = x(4);
syms model_sym(x1, x2, x3, t)
model_sym(x1, x2, x3, t) = x1 .* exp( -(  (t - x2).^2  )./(2 .* x3.^2)   );
model = @(x, t) model_sym(x(1), x(2), x(3), t);
syms r_sym(x1, x2, x3)
r_sym(x1, x2, x3) = model_sym(x1, x2, x3, index) - data
r = @(x) r_sym(x(1), x(2), x(3));
r_temp = matlabFunction(r_sym);
R = @(x) r_temp(x(1), x(2), x(3));
syms f_sym(x1, x2, x3)
f_sym(x1, x2, x3) = (0.5) .* sum(r_sym(x1, x2, x3).^2);
tempf_sym = matlabFunction(f_sym);
f = @(x) tempf_sym(x(1), x(2), x(3));

syms f_A_sym(x1, x2, x3) f_u_sym(x1, x2, x3) f_sigma_sym(x1, x2, x3)
f_A_sym(x1, x2, x3) = diff(f_sym, x1);
f_u_sym(x1, x2, x3) = diff(f_sym, x2);
f_sigma_sym(x1, x2, x3) = diff(f_sym, x3);


syms r_A_sym(x1, x2, x3) r_u_sym(x1, x2, x3) r_sigma_sym(x1, x2, x3)
r_A_sym(x1, x2, x3) = diff(r_sym, x1);
r_u_sym(x1, x2, x3) = diff(r_sym, x2);
r_sigma_sym(x1, x2, x3) = diff(r_sym, x3);


syms J_sym(x1, x2, x3)
J_sym(x1, x2, x3) = [r_A_sym(x1, x2, x3), r_u_sym(x1, x2, x3), r_sigma_sym(x1, x2, x3)];
J_sym_temp = matlabFunction(J_sym);
J = @(x) J_sym_temp(x(1), x(2), x(3));

syms f_AA_sym(x1, x2, x3) f_Au_sym(x1, x2, x3) f_Asigma_sym(x1, x2, x3) f_Aoffset_sym(x1, x2, x3)
syms f_uA_sym(x1, x2, x3) f_uu_sym(x1, x2, x3) f_usigma_sym(x1, x2, x3) f_uoffset_sym(x1, x2, x3)
syms f_sigmaA_sym(x1, x2, x3) f_sigmau_sym(x1, x2, x3) f_sigmasigma_sym(x1, x2, x3) f_sigmaoffset_sym(x1, x2, x3)
f_AA_sym(x1, x2, x3) = diff(f_A_sym, x1);
f_Au_sym(x1, x2, x3) = diff(f_A_sym, x2);
f_Asigma_sym(x1, x2, x3) = diff(f_A_sym, x3);

f_uA_sym(x1, x2, x3) = diff(f_u_sym, x1);
f_uu_sym(x1, x2, x3) = diff(f_u_sym, x2);
f_usigma_sym(x1, x2, x3) = diff(f_u_sym, x3);

f_sigmaA_sym(x1, x2, x3) = diff(f_sigma_sym, x1);
f_sigmau_sym(x1, x2, x3) = diff(f_sigma_sym, x2);
f_sigmasigma_sym(x1, x2, x3) = diff(f_sigma_sym, x3);



syms D2f_sym(x1, x2, x3) invD2f_sym(x1, x2, x3)
D2f_sym(x1, x2, x3) = [f_AA_sym(x1, x2, x3), f_Au_sym(x1, x2, x3), f_Asigma_sym(x1, x2, x3);...
    f_uA_sym(x1, x2, x3), f_uu_sym(x1, x2, x3), f_usigma_sym(x1, x2, x3);...
    f_sigmaA_sym(x1, x2, x3), f_sigmau_sym(x1, x2, x3), f_sigmasigma_sym(x1, x2, x3)]; %%Don't look at this in the console, eyes might be hurt!.
D2f_sym_temp = matlabFunction(D2f_sym);
D2f = @(x) D2f_sym_temp(x(1), x(2), x(3));
%invD2f_sym(x1, x2, x3, x4) = inv(D2f_sym)


%s1 = matlabFunction(J_sym' * J_sym);
%tempA = @(x) s1(x(1), x(2), x(3));

%s2 = matlabFunction((-1 .* J_sym') * r_sym);
%tempB = @(x) s2(x(1), x(2), x(3));
%%syms DeltaX_sym(x1, x2, x3, x4)
%%DeltaX_sym(x1, x2, x3, x4) = linsolve(J_sym(x1, x2, x3, x4)' * J_sym(x1, x2, x3, x4),(-1 .* J_sym(x1, x2, x3, x4)') * r_sym(x1, x2, x3, x4))
%%tempDeltaX = matlabFunction(DeltaX_sym);
%%DeltaX = @(x) tempDeltaX(x(1), x(2), x(3), x(4));
it = 1;
while(1)
   [Q,Rf] = qr(D2f(guess));
   DeltaX = (-1) .* inv(Rf) * inv(Q) * J(guess)' * R(guess); %% (-1) .* pinv(D2f(guess)) * J(guess)' * R(guess);
   guess = guess + DeltaX
   if(f(guess) < 10)
       
       break;
   end
   it = it +1;
   if(it >= 100)
       break;
   end
end

GoddnessofFit = f(guess)
A = guess(1)
mean = guess(2)
sigma = guess(3)







end