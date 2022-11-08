 % Load the images and convert them to grayscale.
 
 
      I1 = rgb2gray(imread('scene_left.png'));
      I2 = rgb2gray(imread('scene_right.png'));
 
      figure; imshowpair(I1,I2,'ColorChannels','red-cyan');
      title('Red-cyan composite view of the stereo images');
 
      % Compute the disparity map.
      disparityMap = disparity(I1, I2, 'BlockSize', 15, ...
        'DisparityRange', [-6 10]);
 
      % For the purpose of visualizing the disparity, replace
      % the -realmax('single') marker with the minimum disparity value.
      marker_idx = (disparityMap == -realmax('single'));
      disparityMap(marker_idx) = min(disparityMap(~marker_idx));
 
      % Show the disparity map. Brighter pixels indicate objects which are
      % closer to the camera.
      figure; imshow(mat2gray(disparityMap));
      colormap jet; colorbar;
      
   %%%%%%dwt%%%%%     
 [A2L1,H2L1,V2L1,D2L1]=dwt2(disparityMap,'haar');
figure;
title('DWT Image');
dwtimg=[A2L1,H2L1;V2L1,D2L1];
imshow(dwtimg,[]);