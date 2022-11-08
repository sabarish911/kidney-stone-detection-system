function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cout=handles.cout;
inp=handles.inp;
save cout cout
a=inp;
figure;
cd datafin
% if isequal(cout,2) 
k=0;
for  m =1:128:size(a,1)
    for n=1:128:size(a,2)
        k=k+1;
        imnew(:,:,k)=a(m:m+127,n:n+127);
        blknam=strcat(num2str(k),'.jpg');
        imwrite(imnew(:,:,k),strcat(num2str(k),'.jpg'));
    end
end
cd ..

for i=1:16
    subplot(4,4,i); 
    blknam=strcat(num2str(i),'.jpg');
    cd datafin
    immg2(:,:,i)=imread(blknam);
    cd ..
    imshow(immg2(:,:,i));
end
% end
 output = ksegment(immg2);
 axes(handles.axes1);
 imshow(output,[]);
if isequal(cout,2) 
 output = ksegment(inp);
elseif isequal(cout,2) 
 output = ksegment(inp);
end
