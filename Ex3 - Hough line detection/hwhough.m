%to use the function
%   hwhough(0.5, input_ex3)

function matA =  hwhough(sigma, Image)
[Gimage,IX,IY,Img] = GoG(sigma, Image);
[imageSx,imageSy] = size(Img);
%create the image from the matrix
G1= mat2gray(Gimage);
level = graythresh(G1); %value 0.2156
% BW = edge(G,'canny');
BW = imbinarize(G1,level);

%Part2
Pmax= ceil(sqrt(imageSx.^2 + imageSy.^2));
Pind = -Pmax:Pmax;
Alpha = -90:89;
H = zeros(ceil(Pmax*2)+1,180);

% extended version
for z=1:imageSx
    for z2=1:imageSy 
        if BW(z,z2)==1
                alpha= atand(IY(z,z2)/IX(z,z2));
                p = z*cosd(alpha) + z2*sind(alpha);% to use the sin and cosine function in degree
                p1 = round(p);
                alpha1 = round(alpha);
                pis = find(Pind == p1);
                alphai = find(Alpha == alpha1);
                H(pis,alphai) = H(pis,alphai) + 1; 
        end
   end
end

% normal version
% for z=1:imageSx
%     for z2=1:imageSy 
%         if BW(z,z2)==1
%             for alpha = -90:89
%                 p = z*cosd(alpha) + z2*sind(alpha);
%                 p = round(p);
%                 pis = find(Pind == p);
%                 alphai = find(Alpha == alpha);
%                 H(pis,alphai) = H(pis,alphai) + 1;                
%             end
%         end
%    end
% end

%Plot Hough Voting Array H
colormap(autumn(5));
imshow(H,[],'XData',Alpha,'YData',Pind,...
            'InitialMagnification','fit');
        axis on, axis normal, hold on;
P  = houghpeaks(double(H),40,'threshold',ceil(0.05*max(H(:))));
kx = Alpha(P(:,2)); ky = Pind(P(:,1));
plot(kx,ky,'s','color','white');
lines = houghlines(BW,Alpha,Pind,P,'FillGap',255,'MinLength',7);
figure, imshow(BW), hold on

for kin = 1:length(lines)
   kxy = [lines(kin).point1; lines(kin).point2];
   plot(kxy(:,1),kxy(:,2),'LineWidth',2,'Color','green'); 
   plot(kxy(1,1),kxy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(kxy(2,1),kxy(2,2),'x','LineWidth',2,'Color','red');
end
hold off

%Use Houghpeaks to find maxima of H
% maxH = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
% imshow(H, [],'XData', Alpha, 'YData', Pind, 'InitialMagnification', 'fit');
% theX = Alpha(maxH(:,2));
% theY = Pind(maxH(:,1));
% axis on, axis normal, hold on;
% plot(theX, theY, 's', 'color', 'white');
% 
% %Use Houghlines to derive corresponding lines
% lines = houghlines(BW,Alpha,Pind,maxH);
% 
% figure, imshow(G1), hold on
% 
% for kx = 1:size(lines)
%     kyx = [lines(kx).point1; lines(kx).point2];
%     plot(kyx(:,1), kyx(:,2), 'LineWidth', 2, 'Color', 'green');
%     %start and end of lines
%     plot(kyx(1,1), kyx(1,2), 'x', 'LineWidth', 2, 'Color', 'red');
%     plot(kyx(2,1), kyx(2,2), 'x', 'LineWidth', 2, 'Color', 'yellow');
% end
% hold off

matA = lines;
end		





        
    
