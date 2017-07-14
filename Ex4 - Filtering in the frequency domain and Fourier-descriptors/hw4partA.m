
% to use the function
%   hw4partA(2.5, taskA)
function mat4 =  hw4partA(sigma, Image)
%converting to grayscale with double valuse between 0 & 1
I =  mat2gray(double(mean(Image, 3)));

%then adding gausian noise to image
Inoise = imnoise(I,'gaussian',0,0.01);
figure, imshow(Inoise);

%Custom 2-d 2d Gaussian filter
% Step 1
r = abs(ceil(3 * sigma));
r2 = (r * 2) + 1;
sizeA = ceil(r2);
counter = r * - 1;
Cx = zeros(sizeA, sizeA);

%create the Cx array

for n = 1:sizeA 
    for m = 1:sizeA
        Cx(n,m) = counter;
        counter = counter + 1;
    end
    counter = r * -1;
end
Cy = Cx.';
% Create the Gmask from Cx and Cy 
Gmask = zeros(sizeA,sizeA);
for i = 1:sizeA 
    for j = 1:sizeA
        Gmask(i,j) = (1*exp(-((Cx(i,j)^2 + (Cy(i,j))^2)/(2*sigma^2))))/(2*pi*sigma^2);
    end    
end


% add the gask to a matrix zero of image size
%image size 
[imageSx, imageSy] = size(I);
% mask size
[imageGx, imageGy] = size(Gmask);
% commpute the difference and add padding of zero after the image
DiffSGx = imageSx - imageGx;
DiffSGy = imageSy - imageGy;
PadG= padarray(Gmask, [DiffSGx DiffSGy],0, 'post');

%Center the filter using function circshift
X = ceil(abs(imageGx/2));
Y = ceil(abs(imageGy/2));
CntPadG = circshift(PadG, [-X, -Y]);

%Transform image and filter kernel to frequency-domain using the function fft2
ImageFFT = fft2(Inoise);
fitlerFFt = fft2(CntPadG);

convIF = ImageFFT.*fitlerFFt;
revIF = ifft2(convIF);
mat4 = revIF;
load trees
subplot(2,2,1), imshow(Inoise); title('Gaussian noise to image')
subplot(2,2,2), imshow(revIF); title('Final Image')
subplot(2,2,3), imshow(mat2gray(Gmask)); title('Gaussian Mask' )
subplot(2,2,4), imshow(mat2gray(CntPadG)); title('Gaussian centered Mask')

end
