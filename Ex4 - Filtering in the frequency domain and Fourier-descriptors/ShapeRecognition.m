%to use the function
%   ShapeRecognition(trainingB, test1B, test2B)

function A4 = ShapeRecognition(image, image1, image2)

%Number of descriptors
ndescriptors = 24;

%Part C: Process image for trainingB.png
I = mat2gray(double(mean(image, 3)));
I=histeq(I);
lvl = graythresh(I);
BW = imbinarize(I,lvl);
BOUND = bwboundaries(BW);
BOUND = BOUND{1};
if mod(size(BOUND, 1), 2) ~= 0
    BOUND = [BOUND; BOUND(end, :)];
end
res1 = partBC(BOUND, ndescriptors);
res1 = res1(2:24);


%Part D: Process image for test1B.jpg
Idoub2 = mat2gray(double(mean(image1, 3)));
lvl2=graythresh(Idoub2);
BW2 = im2bw(Idoub2,lvl2);
BOUND2 = bwboundaries(BW2);
BOUND2 = BOUND2{34};
res2=[50,1];
if mod(size(BOUND2, 1), 2) ~= 0
    BOUND2 = [BOUND2; BOUND2(end, :)];
    
end
res2 = partBC(BOUND2, ndescriptors);
res2 = res2(2:24);

A4 = res2;

%Part D: Process image for test2B.jpg
Idoub3=mat2gray(double(mean(image2, 3)));

lvl3=graythresh(Idoub3);
BW3 = im2bw(Idoub3,lvl3);
BOUND3 = bwboundaries(BW3);
BOUND3 = BOUND3{1};

if mod(size(BOUND3, 1), 2) ~= 0
    BOUND3 = [BOUND3; BOUND3(end, :)];
end

res3 = partBC(BOUND3, ndescriptors);
res3 = res3(2:24);


%Part E
d1 = norm(res1 - res2);
d2 = norm(res1 - res3);
A4= res1;



if(d1 < 0.06)
   figure, imshow(label2rgb(image1, @jet, [.5 .5 .5]));
   hold on
   for k = 1:length(BOUND)
   boundary = BOUND{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
   end
end


% img1 = res1 - min(res1(:));
% img1 = img1 ./ max(img1(:));
% 
% img2 = res2 - min(res2(:));
% img2 = img2 ./ max(img2(:));
% 
% img3 = res3 - min(res3(:));
% img3 = img3 ./ max(img3(:));

% partE = img1 - img2;
% 
% 
% if(partE < 0.06)
%    imshow(img1);
%    for k = 1 :length(partE)
%        currentB = partE{k};
%        x = currentB(:, 2);
%        y = currentB(:, 1);
%        hold on
%        figure, plot(x,y, 'red', 'LineWidth', 2);
%    end
% end


end

function DF1 = partBC(boundaries, n)
    D = boundaries(:,2) + 1i*boundaries(:,1);
    %Fourier descriptor Df
    desc = fft2(D);
    %Translation Invariance
    desc = desc(2:(n+1));
    %Scale Invariance
    desc = desc/(abs(desc(2)));
    %Orientation Invariance
    desc = abs(desc);
    %output
    DF1 = desc;
    
end