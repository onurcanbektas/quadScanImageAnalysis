clear
lib_dir = strcat(pwd, '/MatlabFunctions')
scan_dir = strcat(pwd, '/Scan')
addpath(scan_dir)
addpath(lib_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STARTING LINE OF EDIT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strcat('c150.jpg');
image_inspection_ON = 0;
data_inspection_ON = 0;
fit_ON = 1;

% Do we want to save the result of the fitting parameters to the end of file ?
save_the_fit_parameters = 0;
result_filename = '0308X.txt';
magnetic_field_magnitude = 0.04654967;

%For x -> 0; for y-> 1
x_or_y =  0;
pixel_to_mm_X = 7.07 ./ 87;
pixel_to_mm_Y = 10 ./ 138; %% For 29.07 

%% Cleaning the erroneous data 
numberOfExtraCuts = 1;
z1 = 1;
z2 = 322;
% First extra little cut interval
fit_extra11 = 240; 
fit_extra12 = 270;
% Second extra little cut interval
fit_extra21 = 299;
fit_extra22 = 342;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% END OF LINE OF EDIT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RGP = imread(filename);
image = rgb2gray(RGP);
figure
imshow(image)

data_sum = sum(image, 2);
if(image_inspection_ON == 1)
	[x,y] = ginput(2)
end

i = 1:1:length(data_sum);
plot(i, data_sum, 'x')
if(data_inspection_ON == 1)
    [x,y] = ginput(4)
end

data = data_sum';
index = 1:1:length(data_sum);

%% Extracting the data that we want to use from the given cut intervals
if(numberOfExtraCuts == 0)
    wholerange = index(z1: z2)';
	wholedata = data(z1:z2)';
elseif (numberOfExtraCuts == 1)
		wholerange = [index(z1: fit_extra11), index(fit_extra12: z2)]';
		wholedata = [data(z1:fit_extra11), data(fit_extra12: z2)]';	
elseif (numberOfExtraCuts == 2)
	wholerange = [index(z1: fit_extra11), index(fit_extra12: fit_extra21), index(fit_extra22: z2)]';
	wholedata = [data(z1:fit_extra11), data(fit_extra12: fit_extra21), data(fit_extra22: z2)]';
end

if(x_or_y == 0)
    pixel_to_mm = pixel_to_mm_X;
elseif(x_or_y == 1)
    pixel_to_mm = pixel_to_mm_Y;
end
wholedata = double(wholedata);
wholerange = double(wholerange);



%% Fitting...
%% x = [A , mean, sigma, offset], where
%% our model is A .* exp( -(  (t - mean.^2  )./(2 .* sigma.^2)   )  + offset;
if(fit_ON == 1)
t = [max(wholedata), mean(wholerange), 20, min(wholedata)];
t = t(:);
model = @(x) x(1) .* exp( -(  (wholerange - x(2)).^2  )./(2 .* (x(3).^2))   )  + x(4) - wholedata;
[soln ,ssq,cnt, res, XY] = LMFnlsq(model,t,'Display',1,'XTol', 1e-7,'Lambda', 1);


%% Plotting the fit against data
r = z1:1:z2;
figure(3)
plot(wholerange, wholedata, 'x')
hold on
plot(r, soln(1).*exp( -(  (r - soln(2)).^2  )./(2 .* soln(3).^2)  ) + soln(4))
hold off

A = soln(1)
mean = soln(2) * pixel_to_mm
sigma = soln(3) * pixel_to_mm
offset = soln(4)

if(save_the_fit_parameters == 1)
	pFile = fopen (result_filename,'a+');
	fprintf(pFile, '\n%s %f %.5f %d %d %d %d %d %d %d',filename, sigma, magnetic_field_magnitude, numberOfExtraCuts, z1, z2, fit_extra11, fit_extra12, fit_extra21, fit_extra22)
end
end











