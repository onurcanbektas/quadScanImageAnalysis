function [image] = loadImageFile (filename, image_inspection_ON)

	RGP = imread(filename);
	image = rgb2gray(RGP);
	figure(1)
	imshow(image)
	if (image_inspection_ON == 1)
	    [x,y] = ginput(4)
	    return
	end

end