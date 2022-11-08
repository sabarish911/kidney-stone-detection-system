function [AA12, AA22,BinaryImage] = Lclusteringfin(Input,Clusterfin)
load file file
Input1 = imread('fina.jpg'); 
inp=imresize(Input1,[512 512]);
 
   if size(inp,3)>1
     inp = rgb2gray(inp);
   end
cout12 = str2num(file(1:2));
Clusters=2;
load finval finval
Input=inp;
Input=imcrop(Input,finval(cout12,:));
% Lclusteringfin(Input,Clusters,img2);
[r c] = size(Input);
Length  = r*c; 
wd1=r;
wd2=c;

a_mn=finval(cout12,1);
b_mn=finval(cout12,2);
a_mx=a_mn+finval(cout12,3);
b_mx=a_mn+finval(cout12,4);



Dataset = reshape(Input,[Length,1]);   %%%%%Reshape 2D Image to 1D Vectors 
    
Cluster1=zeros(Length,1);
Cluster2=zeros(Length,1);

miniv = min(min(Input));      
maxiv = max(max(Input));
range = maxiv - miniv;
stepv = range/Clusters;
incrval = stepv;

for i = 1:Clusters            %%%%Find the centroids to each Clusters
    K(i).centroid = incrval;
    incrval = incrval + stepv;
end

update1=0;
update2=0;


mean1=2;
mean2=2;

while  (mean1 ~= update1) && (mean2 ~= update2)

mean1=K(1).centroid;
mean2=K(2).centroid;

for i=1:Length                     %%%%%%Find the distance between Each Pixel and Centroids 
    for j = 1:Clusters
        temp = Dataset(i);
        difference(j) = abs(temp-K(j).centroid);
    end
    [y,ind]=min(difference);     %%%%Group Pixels to Each Cluster Based on Minimum Distance
    
	if ind==1
        Cluster1(i)   =temp;
	end
    if ind==2
        Cluster2(i)   =temp;
    end
end

%%%%%UPDATE CENTROIDS
cout1=0;
cout2=0;

for i=1:Length
    Load1=Cluster1(i);
    Load2=Cluster2(i);
    
    if Load1 ~= 0
        cout1=cout1+1;
    end
    
    if Load2 ~= 0
        cout2=cout2+1;
    end
end

Mean_Cluster(1)=sum(Cluster1)/cout1;
Mean_Cluster(2)=sum(Cluster2)/cout2;

for i = 1:Clusters
    K(i).centroid = Mean_Cluster(i);

end

update1=K(1).centroid;
update2=K(2).centroid;
end


AA12=reshape(Cluster1,[wd1 wd2]);
AA22=reshape(Cluster2,[wd1 wd2]);

% 
% figure('Name','Segmented Results');
% subplot(2,2,1); imshow(AA12,[]);
% subplot(2,2,2); imshow(AA22,[]);

[r112,c112]=size(AA22);
 %Convert to Binary Image
for i=1:r112
    for j=1:c112
    if AA22(i,j) >150
        BinaryImage(i,j)=255;
    else
        BinaryImage(i,j)=0;
    end
    end
end

a_mn=finval(cout12,2);
b_mn=finval(cout12,1);
a_mx=a_mn+finval(cout12,4);
b_mx=b_mn+finval(cout12,3);
finval(11,:)
imgfinn=zeros(512,512);

    for Cx =a_mn:a_mx-1                              
        for Cy11 =b_mn:b_mx-1 
            imgfinn(Cx,Cy11) = BinaryImage(Cx-a_mn+1,Cy11-b_mn+1);
        end
    end 
imwrite(imgfinn,'1.png');
return;
