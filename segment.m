function  F = segment(inp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nclusters=4;
for i=1:16
%     subplot(4,4,i); 
    blknam=strcat(num2str(i),'.jpg');
    cd datafin
    immg2(:,:,i)=imread(blknam);
    cd ..
end
k=0;
% cata=Input1/sum(sum(sum(Input1)));
cd datafin
% if isequal(cout,2) 
for  m =1:128:512
    for n=1:128:512
        k=k+1;
        blknam=strcat(num2str(k),'.jpg');
        inp11(m:m+127,n:n+127)=immg2(:,:,k);
    end
end
cd ..

[AA1, AA2, AA3, AA4] = spatialclust(inp11,nclusters);
   cd Clusim
   file = uigetfile('*.bmp','pick file');
   segout= imread(file);
   cd ..
   output=imread('1.png');
    % dilation
    se = strel('ball',5,5);
    F=imdilate(output,se);
end

