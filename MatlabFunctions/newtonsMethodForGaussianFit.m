

function [GoddnessofFit, A,mean, sigma, offset] = newtonsMethodForGaussianFit(data, index)

%%index = [1, 4, 7]';
%%data = [2 5 4]';
guess = [250.0, 250.0, 40.0 ,20.0]';
cut_edit1 = 5;
cut_edit2 = 15;
cut_edit3 = 35;
cut_edit4 = 55;
data = [data(cut_edit1:cut_edit2); data(cut_edit3:cut_edit4)];
index = [index(cut_edit1:cut_edit2);index(cut_edit3:cut_edit4)];

%% x's are the parameters that are to be fitted with the data, and t's are the data itself
%% A = x(1); mean = x(2); sigma = x(3); offset = x(4);
syms model_sym(x1, x2, x3, x4, t)
model_sym(x1, x2, x3, x4, t) = x1 .* exp( -(  (t - x2).^2  )./(2 .* x3.^2)   )  + x4;
model = @(x, t) model_sym(x(1), x(2), x(3), x(4), t);
syms r_sym(x1, x2, x3, x4)
r_sym(x1, x2, x3, x4) = model_sym(x1, x2, x3, x4, index) - data;
r = @(x) r_sym(x(1), x(2), x(3), x(4));
r_temp = matlabFunction(r_sym);
R = @(x) r_temp(x(1), x(2), x(3), x(4));
syms f_sym(x1, x2, x3, x4)
f_sym(x1, x2, x3, x4) = (0.5) .* sum(r_sym(x1, x2, x3, x4).^2);
tempf_sym = matlabFunction(f_sym);
f = @(x) tempf_sym(x(1), x(2), x(3), x(4));

%{
syms f_A_sym(x1, x2, x3, x4) f_u_sym(x1, x2, x3, x4) f_sigma_sym(x1, x2, x3, x4) f_offset_sym(x1, x2, x3, x4)
f_A_sym(x1, x2, x3, x4) = diff(f_sym, x1);
f_u_sym(x1, x2, x3, x4) = diff(f_sym, x2);
f_sigma_sym(x1, x2, x3, x4) = diff(f_sym, x3);
f_offset_sym(x1, x2, x3, x4) = diff(f_sym, x4);


syms f_AA_sym(x1, x2, x3, x4) f_Au_sym(x1, x2, x3, x4) f_Asigma_sym(x1, x2, x3, x4) f_Aoffset_sym(x1, x2, x3, x4)
syms f_uA_sym(x1, x2, x3, x4) f_uu_sym(x1, x2, x3, x4) f_usigma_sym(x1, x2, x3, x4) f_uoffset_sym(x1, x2, x3, x4)
syms f_sigmaA_sym(x1, x2, x3, x4) f_sigmau_sym(x1, x2, x3, x4) f_sigmasigma_sym(x1, x2, x3, x4) f_sigmaoffset_sym(x1, x2, x3, x4)
f_offsetX_sym = [0 0 0 0];
f_AA_sym(x1, x2, x3, x4) = diff(f_A_sym, x1);
f_Au_sym(x1, x2, x3, x4) = diff(f_A_sym, x2);
f_Asigma_sym(x1, x2, x3, x4) = diff(f_A_sym, x3);
f_Aoffset_sym(x1, x2, x3, x4) = diff(f_A_sym, x4);
f_uA_sym(x1, x2, x3, x4) = diff(f_u_sym, x1);
f_uu_sym(x1, x2, x3, x4) = diff(f_u_sym, x2);
f_usigma_sym(x1, x2, x3, x4) = diff(f_u_sym, x3);
f_uoffset_sym(x1, x2, x3, x4) = diff(f_u_sym, x4);
f_sigmaA_sym(x1, x2, x3, x4) = diff(f_sigma_sym, x1);
f_sigmau_sym(x1, x2, x3, x4) = diff(f_sigma_sym, x2);
f_sigmasigma_sym(x1, x2, x3, x4) = diff(f_sigma_sym, x3);
f_sigmaoffset_sym(x1, x2, x3, x4) = diff(f_sigma_sym, x4);


syms D2f_sym(x1, x2, x3, x4) invD2f_sym(x1, x2, x3, x4)
D2f_sym(x1, x2, x3, x4) = [f_AA_sym(x1, x2, x3, x4), f_Au_sym(x1, x2, x3, x4), f_Asigma_sym(x1, x2, x3, x4), f_Aoffset_sym(x1, x2, x3, x4);...
    f_uA_sym(x1, x2, x3, x4), f_uu_sym(x1, x2, x3, x4), f_usigma_sym(x1, x2, x3, x4), f_uoffset_sym(x1, x2, x3, x4);...
    f_sigmaA_sym(x1, x2, x3, x4), f_sigmau_sym(x1, x2, x3, x4), f_sigmasigma_sym(x1, x2, x3, x4), f_sigmaoffset_sym(x1, x2, x3, x4);...
    f_offsetX_sym]; %%Don't look at this in the console, eyes might be hurt!.
D2f_sym_temp = matlabFunction(D2f_sym);
D2f = @(x) D2f_sym_temp(x(1), x(2), x(3), x(4));

%invD2f_sym(x1, x2, x3, x4) = inv(D2f_sym)
%}



syms r_A_sym(x1, x2, x3, x4) r_u_sym(x1, x2, x3, x4) r_sigma_sym(x1, x2, x3, x4) r_offset_sym(x1, x2, x3, x4)
r_A_sym(x1, x2, x3, x4) = diff(r_sym, x1);
r_u_sym(x1, x2, x3, x4) = diff(r_sym, x2);
r_sigma_sym(x1, x2, x3, x4) = diff(r_sym, x3);
r_offset_sym(x1, x2, x3, x4) = diff(r_sym, x4);

syms J_sym(x1, x2, x3, x4)
J_sym(x1, x2, x3, x4) = [r_A_sym(x1, x2, x3, x4), r_u_sym(x1, x2, x3, x4), r_sigma_sym(x1, x2, x3, x4), r_offset_sym(x1, x2, x3, x4)];
J_sym_temp = matlabFunction(J_sym);
J = @(x) J_sym_temp(x(1), x(2), x(3), x(4));

s1 = matlabFunction(J_sym' * J_sym);
tempA = @(x) s1(x(1), x(2), x(3), x(4));

s2 = matlabFunction((-1 .* J_sym') * r_sym);
tempB = @(x) s2(x(1), x(2), x(3), x(4));
%%syms DeltaX_sym(x1, x2, x3, x4)
%%DeltaX_sym(x1, x2, x3, x4) = linsolve(J_sym(x1, x2, x3, x4)' * J_sym(x1, x2, x3, x4),(-1 .* J_sym(x1, x2, x3, x4)') * r_sym(x1, x2, x3, x4))
%%tempDeltaX = matlabFunction(DeltaX_sym);
%%DeltaX = @(x) tempDeltaX(x(1), x(2), x(3), x(4));
it = 1;
while(1)
   DeltaX = tempA(guess)\tempB(guess);
   %d2 = D2f(guess);
   %if(any(isfinite(d2(:)))|any(isnan(d2(:))))
   %    d2 = d2 + eye(4)*10^(-2);
   %end
   %%DeltaX = (+1) .* pinv(D2f(guess)) * J(guess)' * R(guess);
   guess = guess - DeltaX
   if(f(guess) < 10)
       break;
   end
   it = it +1
   %%if(it == 723)
   %    continue;
   %end
   if(it >= 100)
       break;
   end
end

GoddnessofFit = f(guess)
A = guess(1)
mean = guess(2)
sigma = guess(3)
offset = guess(4)







%{
syms r_AA_sym(x1, x2, x3, x4) r_Au_sym(x1, x2, x3, x4) r_Asigma_sym(x1, x2, x3, x4) r_Aoffset_sym(x1, x2, x3, x4)
syms r_uA_sym(x1, x2, x3, x4) r_uu_sym(x1, x2, x3, x4) r_usigma_sym(x1, x2, x3, x4) r_uoffset_sym(x1, x2, x3, x4)
syms r_sigmaA_sym(x1, x2, x3, x4) r_sigmau_sym(x1, x2, x3, x4) r_sigmasigma_sym(x1, x2, x3, x4) r_sigmaoffset_sym(x1, x2, x3, x4)
r_AA_sym(x1, x2, x3, x4) = diff(r_A_sym, x2);
%}


%%Redundant
%{ 
%% x's are the parameters that are to be fitted with the data, and t's are the data itself
%% A = x(1); mean = x(2); sigma = x(3); offset = x(4);
model = @(x, t) x(1) .* exp( -(  (t - x(2)).^2  )./(2 .* x(3).^2)   )  + x(4);
r = @(x) model(x, index) - data;
R = @(x) r(x);
f = @(x) (0.5) .* sum(R(x).^2);


%%Partial derivatives of r_i's. Notice t_i = index(:,i);
r_A = @(x) exp( -(  (index - x(2)).^2  )./(2 .* x(3).^2)   );
r_u = @(x) ( x(1) ./ x(3).^2 ) .* (index - x(2)) .* exp( -(  (index - x(2)).^2  )./(2 .* x(3).^2)   );
r_sigma = @(x) ( ( x(1) .* (index - x(2)).^2 ) ./ x(3).^3 ) .* exp( -(  (index - x(2)).^2  )./(2 .* x(3).^2)   );
r_offset = ones(1,length(index))';
J = @(x) [r_A(x), r_u(x), r_sigma(x), r_offset];
%%Second partial derivative of r_i's.
r_AA = @(x) 0;
r_uA = @(x) r_u(x) ./ x(1);
r_sigmaA = @(x) r_sigma(x) ./ x(1);
r_offsetA = @(x) 0;

r_Au = @(x) r_u(x) ./ x(1);
r_uu = @(x) x(1) .*  exp( -(  (index - x(2)).^2  )./(2 .* x(3).^2)   ) .* (index.^2  - 2.* index .* x(2) + x(2).^2 - x(3).^2)./ x(3).^4;
%%r_sigmau = @(x) 

%%D2r = @(x) [0, ];

%%laplacianr = @(x) 0 + 
S = @(x) dot(r(x), laplacianr(x));
laplacianf = @(x) J(x)' * J(x) + S(x);

%%model([26, 150, 38 ,10], index)
%%r(guess)
%%R([1, 3, 1 ,0])
%%f(guess)
%%J(guess)
%%laplacianr([1, 3, 1 ,0])
%%S([1, 3, 1 ,0])
%5J([1, 3, 1 ,0])' * J([1, 3, 1 ,0])
%%laplacianf([1, 3, 1 ,0])
%%deltaf2 = @(x) (-1) .* ( J(x)' * J(x) +           )


%%syms x
%%matlabFunction( diff(r_A(x)) )

syms f(t, x1, x2, x3, x4)
f(t, x1, x2, x3, x4) = x1 .* exp( -(  (t - x2).^2  )./(2 .* x3.^2)   )  + x4;
df = diff(f,x2)
f_x2 = matlabFunction(df)
g = @(x) f_x2(index, x(1), x(2), x(3), x(4))
g(guess)
%%f_x2(index, guess(1), guess(2), guess(3), guess(4))
%}



end