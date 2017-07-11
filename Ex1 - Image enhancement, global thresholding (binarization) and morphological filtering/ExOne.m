%Team: Flaming Corgis
%Members: Liselot Ramirez, Wael Al Atrash, Nazmus Saquib

%How to call our main function: import the image, and then specify the
%Strel, for example: ExOne(Unbenannt, strel('disk',4))

%EXECUTING TASK D To use all three tasks Together
function [resultA, resultB, resultC] = ExOne(inputImg, ourStrel)
    resultA = stepA(inputImg); %enhacing contrast
    assignin('base','resultA',resultA);
    resultB = stepB(resultA); %thresholding
    assignin('base','resultB',resultB);
    resultC = stepC(resultB, ourStrel);%morphological filtering
    assignin('base','resultC',resultC);
end

%EXECUTING TASK A Image Enhancement
function visualA = stepA(n)
%Visualize the img
    imshow(n);
    %converts image from RBG to GreyScale
    blackImg = uint8(mean(n, 3));
    imshow(blackImg);
    imhist(blackImg);
    %Enhance contrast using histogram equalization (contrast)
	%the intial histogram was spread to fill the full 255 bits
    contrastedImg=histeq(blackImg);
    imhist(contrastedImg);
    visualA = contrastedImg;
end

%EXECUTING TASK B thresholding
function visualB = stepB(m)
    level = graythresh(m); %Finds the threshold level
    BW = im2bw(m,level); %converts it to black and white
    BW = im2bw(m, 0.44); imshowpair(m,BW,'montage'); %use this to test different values and see the difference best value found.44
    inverseBW = imcomplement(BW); %done again to reverse the black and white
    imshowpair(BW,inverseBW,'montage');
    visualB = inverseBW;
end

%EXECUTING TASK C Binary mask and Morphological operators
function visualC = stepC(img, st)
    %Gets the neighboring strel to compare it to the pixels in our image
    %In our case, a disk size 4 on 2D (strel('disk',4))
    window=getnhood(st);

    %Given the neighborhood value, which is the comparison window,
    %we specify in which position we need our pixel matrix to grow and with
    %which value to fill it
    %vector position in which the matrix will grow
    vpos=floor(size(window,1)/2);
    %elements to add on the matrix
    elemAdd=floor(size(window,2)/2);

    %Pad array on sides
    %We take the image as the initial matrix (img), we setup the size and how many
    %elements to add ([vpos elemAdd]), and finally we fill it with ones (1)
    modMatrix=padarray(img,[vpos elemAdd],1);
    %Intialize an empty matrix with size as the image
    emptyMatrix=false(size(img));

    %Now we compare and replace the pixels if they match the window
    for i=1:size(modMatrix,1)-(2*vpos)
        for j=1:size(modMatrix,2)-(2*elemAdd)
           %Temp will create a vector that will take the current vector
           %position , or row, called (i) and it will increase based on our
           % vpos variable; then it will take the element added, or col, (j)
           % and it will increase it as well
            Temp=modMatrix(i:i+(2*vpos),j:j+(2*elemAdd));

            %We fill the empty matrix in case of a match
            emptyMatrix(i,j)=min(min(Temp-window));

        end
    end
    %Step C-D: Our erosion takes more area than the native function imopen from matlab
    visualC = imcomplement(emptyMatrix);
end

%EXECUTING TASK E
%Are the results satisfactory? As for humans go, they are relatively
%satisfactory, but a better job can be done in the future.
%What are the limitations of this approach for separating background and
%foreground? Limitations are the handling the smaller gaps in the high
%definition pictures; the more erosion and/or dilation the contouring that
%should exist gets distorted.
