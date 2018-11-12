%% @author: Onurcan Bektaş
%%          onurcan.bektas@metu.edu.tr
clear
clear all

lib_dir = strcat(pwd, '/matlabFunctions');
scan_dir = strcat(pwd, '/scan.29072017/');
output_dir = strcat(pwd, '/output.29072017/');
addpath(scan_dir)
addpath(lib_dir)
addpath(output_dir)


i=1;
pixelToMMRatioInX = 10 ./ 86;


RGP = imread([scan_dir, imageFileName]);
image = rgb2gray(RGP);


saveFitResults = false;


if imageInspection_ON == true
	figure(1)
	imshow(image)
	[x,y] = ginput(3)
end

if isItXaxis == false
	data = sum(image, 1);
	sizeOfData = length(data);
	pxToMMRatio = pixelToMMRatioInX;
else
	data = sum(image, 2);
	data = data';
	sizeOfData = length(data);
	pxToMMRatio = pixelToMMRatioInY;
end


if (dataInspection_ON == true)
	figure(2)
	plot(xAxis, data, 'x')
	[x,y] = ginput(2)
end


numberOfCutsToMake = 0;
dataIndexBeg = 1;
dataIndexEnd = sizeOfData;

firstCutBeg = 250;
firstCutEnd = 308;

secondCutBeg = 421;
secondCutEnd = 557;

index = 1:1:sizeOfData;
switch numberOfCutsToMake
	case 0
		cleanIndex = dataIndexBeg:1:dataIndexEnd;
		cleanData = data(dataIndexBeg:dataIndexEnd);
	case 1
		cleanIndex = [dataIndexBeg:1:firstCutBeg, firstCutEnd:1:dataIndexEnd];
		cleanData = [data(dataIndexBeg: firstCutBeg), data(firstCutEnd:dataIndexEnd)];
	case 2
		cleanIndex = [dataIndexBeg:1:firstCutBeg, firstCutEnd:1:secondCutBeg, secondCutEnd:1:dataIndexEnd];
		cleanData = [data(dataIndexBeg:firstCutBeg), data(firstCutEnd:secondCutBeg), data(secondCutEnd:dataIndexEnd)];
	otherwise
		"Dude, I cannot make that much cut!"
end

cleanIndex = double(cleanIndex)';
cleanData = double(cleanData)';

if makeFit_ON == true
	t = [max(cleanData), mean(cleanIndex), 20, min(cleanData)];
	t = t(:);
	model = @(x) x(1) .* exp( -(  (cleanIndex - x(2)).^2  )./(2 .* (x(3).^2))   )  + x(4) - cleanData;
	[soln ,ssq,cnt, res, XY] = LMFnlsq(model,t,'Display',1,'XTol', 1e-7,'Lambda', 1);

	figure(3)
	plot(cleanIndex, cleanData, 'x')
	hold on
	plot(index, soln(1).*exp( -(  (index - soln(2)).^2  )./(2 .* soln(3).^2)  ) + soln(4))
	legend('Veri', 'fit')
	xlabel('Görüntü boyutu (pixel)')
	ylabel('Parlaklık değeri (birimsiz)')
	title('Parlaklık vs görüntü boyutu grafiği; B = 2.32 kG')
	hold off

	A = soln(1)
	mean = soln(2) * pxToMMRatio
	sigma = soln(3) * pxToMMRatio
	offset = soln(4)

	if(saveFitResults == true)
	    save_dir = strcat(output_dir, outputFileName);
		pFile = fopen (save_dir, 'a+');
		magneticFields = [0.01567, 0.03088, 0.04655, 0.06222, 0.07743, 0.09310, 0.10832, 0.12398, 0.13965, 0.15531, 0.17053, 0.18620, 0.20142, 0.21708, 0.23230, 0.24797, 0.26363, 0.27885, 0.29452, 0.31018, 0.32540, 0.34107];
	end

end


