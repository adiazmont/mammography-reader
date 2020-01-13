%myFolder = '/Users/adiaz/Documents/Research/Breast_Cancer_Image_Classification/Mammography/Cancer';
myFolder = '/Users/adiaz/Documents/MATLAB';
if ~isdir(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s',myFolder);
    uiwait(warndlg(errorMessage));
    return;
end
filePattern = fullfile(myFolder,'*.txt');
filePattern
size(filePattern)
overlayFiles = dir(filePattern);
number_file = length(overlayFiles);
 
%reportFile = fopen('result2.txt','w');   %In order to writer result, create the result2.txt
 
for FileNum = 1:number_file;
    baseFileName = overlayFiles(FileNum).name;
    fullFileName = fullfile(myFolder,baseFileName);
    fprintf(1,'Now reading %s\n',fullFileName);
    %fprintf(file_1,'Now reading %\n',fullFileName); %write to result file 
    
    fopen(fullFileName);
    fid = fopen(fullFileName);
% InputText = textscan(fid,'%s',8,'delimiter','\n');
% InputText
% Intro = InputText{1};
% Intro
%disp(Intro);                   %acquire the introduction in the every file
InputText=textscan(fid,'%d');
chaincode = InputText{1};      %obtain the Chaincode, here Chaincode is array format.
% for chaincode_index = 1:length(chaincode);
%     chaincode(chaincode_index)
%     str_chaincode_value = char(chaincode(chaincode_index));
%     class(str_chaincode_value)
%     double_chaincode_value = str2num(str_chaincode_value);
%     int32_chaincode_value = int32(double_chaincode_value);
%     chaincode(chaincode_index) = int32_chaincode_value;
% end
%disp(Chaincode);
[row column]=size(chaincode);  % the row of Chaincode denotes the length of Chaincode.
%for j = 1:row
%   Chaincode(j,1)    %if need to access some Chaincode, we can use Chaincode(i,j) to access that. 
%end
s2=sqrt(2);
perimeter=0;   %the perimeter of the breast tumor shape.              
area=0;        %the area of the breast tumor shape.    
y=0;  
%for j=3:row    %compute the area and the perimeter of breast tumor shape.
 %   if Chaincode(j,1)==0
 %      y=y-1;
 %      perimeter=perimeter+1;  %.......
 %   elseif Chaincode(j,1)==1
 %  end
 % end             %the end of computing the area and the perimeter of breast tumor shape.
%from chaincode to get the edge point
edge_point(1,1)=chaincode(1,1);
edge_point(1,2)=chaincode(2,1);
for j=3:row
    if chaincode(j,1)==0;
    edge_point(j-1,1)=edge_point(j-2,1);
    edge_point(j-1,2)=edge_point(j-2,2) +1;
    elseif chaincode(j,1)==1;
    edge_point(j-1,1)=edge_point(j-2,1)-1;
    edge_point(j-1,2)=edge_point(j-2,2)+1;
    elseif chaincode(j,1)==2;
    edge_point(j-1,1)=edge_point(j-2,1)-1;
    edge_point(j-1,2)=edge_point(j-2,2);
    elseif chaincode(j,1)==3;
    edge_point(j-1,1)=edge_point(j-2,1)-1;
    edge_point(j-1,2)=edge_point(j-2,2)-1;
    elseif chaincode(j,1)==4;
    edge_point(j-1,1)=edge_point(j-2,1);
    edge_point(j-1,2)=edge_point(j-2,2)-1;
    elseif chaincode(j,1)==5;
    edge_point(j-1,1)=edge_point(j-2,1)+1;
    edge_point(j-1,2)=edge_point(j-2,2)-1;
    elseif chaincode(j,1)==6;
    edge_point(j-1,1)=edge_point(j-2,1)+1;
    edge_point(j-1,2)=edge_point(j-2,2);
    else chaincode(j,1)==7;
    edge_point(j-1,1)=edge_point(j-2,1)+1;
    edge_point(j-1,2)=edge_point(j-2,2)+1;
    end
end

%edge_point(row-1,1)
%edge_point(row-1,2)
%The end of geting the edge point. Here, one thing that should be put
%attention is the last point row-1 is as same as the first point. That is
%to say, the Chaincode generates a complete shape. So if we want to use the
%edge point, we should be careful the length of edge point which should be
%row-2, not row-1. Also, we use edge_point(i,1) and edge_point(i,2) to get
%the edge point.
edge_point_min_x=edge_point(1,1);
edge_point_min_y=edge_point(1,2);
for i=1:row-2
    if edge_point(i,1)<edge_point_min_x;
        edge_point_min_x=edge_point(i,1);
    end
    if edge_point(i,2)<edge_point_min_y;
        edge_point_min_y=edge_point(i,2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%edge_point_min_x and edge_point_min_y can be considered as the new origin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%1.Fourier irregularity index
for i=1:row-2
    x_dist=edge_point(i,1)-edge_point_min_x; % the distance in x 
    y_dist=edge_point(i,2)-edge_point_min_y; % the distance in y
    x_y_dist=x_dist*x_dist+y_dist*y_dist;
    r(i)=sqrt(double(x_y_dist));    %r(i) is the distance between the edge_point to minimum point
end
FII_R =real(fft(r)); %%%%%%%%%%%%%% maybe need to change ifft to be fft for real input.be careful!
X = FII_R;


X = double(X);
X = X';
size(X);

% Centroid and Plolygeom functions
% xColumn = edge_point(:,1);
% yColumn = edge_point(:,2);
% [center_x,center_y,a] = centroid(xColumn,yColumn);
% round(center_x)
% round(center_y)
% 
% [ geom, iner, cpmo ] = polygeom(xColumn,yColumn);
% geom

k=2;


%figure(3)
%subplot(4,1,1)
%hist(X,100);
%figure(2)
% [W,M,V,L] = EM_GM(X,k,0.01,1000,1,[]);
% W
% M
% V
% L
% 
% 
% L = MyLikelihood(X,k,W,M,V)

%write the output into result.txt 
%fileID=fopen('result.txt','w');
%fprintf(fileID,'%f      ',FII_R);
%fclose(fileID);
%the end of the writing the output into result.txt

%plot(R)
%axis([0 1800 -0.3 0.3]);

%2.Original roughness index    here, we assume that the center point is
%constant. Of course, we need to change that so it is right. 
center_point_x=1500;
center_point_y=2000;
for i=1:row-2         %compute the distance from edge point to center point
    x_dist=abs(edge_point(i,1)-center_point_x);
    y_dist=abs(edge_point(i,2)-center_point_y);
    d(i) = sqrt(double(x_dist*x_dist+y_dist*y_dist));
end
d_max=d(1);           %obtain the maximum d so that we can normalize distance
for i=1:row-2
    if d_max<d(i);
        d_max=d(i);
    end
end                   %the end of the obtaining the max d

for i=1:row-2
    d(i)=d(i)/d_max;
end
for i=2:row-2            %compute the difference between two distance from edge
   RI_R(i)=d(i)-d(i-1);  %point to center point
end
RI_R(1)=d(1)-d(row-2);
RI_R;
X = RI_R;

size(X);
X = double(X);
X = X';


%k=2;

%figure(3)
%subplot(4,1,1)
%hist(X,100);
%figure(2)
% [W,M,V,L] = EM_GM(X,k,0.01,1000,1,[]);
% W
% M
% V
% L


%L = MyLikelihood(X,k,W,M,V)

%The end of the Original Roughness Index,Here RI_R(i) is original roughness index, we can
%use that for next step.

%3.Enhance Roughness Index
%here,we use some existed variables in the above Original roughness index computing
for i=3:row-2
    ERI_r(i)=d(i)-d(i-2);    %adjust the step length to determine which step length is the best.
end
ERI_r(1)=d(1)-d(row-3);
ERI_r(2)=d(2)-d(row-2);
ERI_R=real(ifft(ERI_r));     %apply ifft on the adjusted distance difference.
%The end of computing Enhance Roughness Index.
ERI_R;
X = ERI_R;

size(X);
X = double(X);
X = X';

%k=2;

% figure(3)
% subplot(4,1,1)
% hist(X,100);
% figure(2)
% [W,M,V,L] = EM_GM(X,k,0.01,1000,1,[]);
% W
% M
% V
% L


L = MyLikelihood(X,k,W,M,V)

disp('= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ');
fclose('all');    
end
% for fileNum = 1:1
%     baseFileName = overlayFiles(fileNum).name;
%     fullFileName = fullfile(myFolder,baseFileName);
%     fprintf(1,'Now reading %s\n',fullFileName);
%     fprintf(reportFile,'Now reading %s\n',fullFileName); %write to result file 
%     
%     fid = fopen(fullFileName);
%     inputText = textscan(fid,'%s',8,'delimiter','\n');
%     intro = inputText{1};          %acquire the introduction in the every file              
%     inputText=textscan(fid,'%d');
%     chainCode = inputText{1};      %obtain the Chaincode, here Chaincode is array format.
%     fclose(fid);
%     
%     chainCodeLength = size(chainCode);
%     
%     edge_point(1,1)=chainCode(1,1);
%     edge_point(1,2)=chainCode(2,1);
%     
%     % Obtaining the new origin
%     min_values = min(edge_point,[],1);
%     edge_point_min_x = min_values(1,1);
%     edge_point_min_y = min_values(1,2);
%     
%     %1.Fourier irregularity index
%     for i=1:chainCodeLength-2
%         x_dist=edge_point(i,1)-edge_point_min_x; % the distance in x 
%         y_dist=edge_point(i,2)-edge_point_min_y; % the distance in y
%         x_y_dist=x_dist*x_dist+y_dist*y_dist;
%         r(i)=sqrt(double(x_y_dist));    %roughness index
%                                         %r(i) is the distance between the edge_point to minimum point
%     end
%     
%     FII_R =real(ifft(r)); %%%%%%%%%%%%%% maybe need to change ifft to be fft for real input.be careful!
%     X = FII_R;
%     plot_(1) = subplot(2,1,1);
%     plot(r);
%     title(plot_(1),{'Roughness Index';'Distance between adjacent points'});
%     xlabel(plot_(1),'i''th edge point px');
%     ylabel(plot_(1),'y value px');
%     
%     X = double(X);
%     plot_(2) = subplot(2,1,2);
%     plot(X);
%     title(plot_(2),{'Enhanced Roughness Index';'Irregularity of breast tumor'})
%     xlabel(plot_(1),'i''th edge point');
% 
%     X = double(X);
%     X = X';
%     size(X);
%     
%     k=3;
% 
% 
%     figure(3)
%     subplot(4,1,1)
%     hist(X,100);
%     figure(2)
%     [W,M,V,L] = EM_GM(X,k,0.01,1000,1,[]);
%    
%      L = MyLikelihood(X,k,W,M,V);
%      
%      W
%      M
%      V
%      L
% %     
% %     %write the output into result.txt 
% %     fileID=fopen('result.txt','w');
% %     fprintf(fileID,'%f      ',FII_R);
% %     fclose(fileID);
%     %the end of the writing the output into result.txt
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %2.Original roughness index    here, we assume that the center point is
% %constant. Of course, we need to change that so it is right. 
% 
% % center_point_x=1500;
% % center_point_y=2000;
% % for i=1:chainCodeLength-2         %compute the distance from edge point to center point
% %     x_dist=edge_point(i,1)-center_point_x;
% %     y_dist=edge_point(i,2)-center_point_y;
% %     d(i) = sqrt(double(x_dist*x_dist+y_dist*y_dist));
% % end
% % d_max=d(1);           %obtain the maximum d so that we can normalize distance
% % for i=1:chainCodeLength-2
% %     if d_max<d(i);
% %         d_max=d(i);
% %     end
% % end                   %the end of the obtaining the max d
% % 
% % for i=1:chainCodeLength-2
% %     d(i)=d(i)/d_max;
% % end
% % for i=2:chainCodeLength-2            %compute the difference between two distance from edge
% %    RI_R(i)=d(i)-d(i-1);  %point to center point
% % end
% % RI_R(1)=d(1)-d(chainCodeLength(1)-2);
% % %RI_R;
% % X = RI_R;
% % 
% % size(X);
% % X = double(X);
% % X = X';
% 
% 
% %k=2;
% 
% %figure(3)
% %subplot(4,1,1)
% %hist(X,100);
% %figure(2)
% % [W,M,V,L] = EM_GM(X,k,0.01,1000,1,[]);
% % 
% % 
% % L = MyLikelihood(X,k,W,M,V)
% % 
% % %The end of the Original Roughness Index,Here RI_R(i) is original roughness index, we can
% % %use that for next step.
% % 
% % %3.Enhance Roughness Index
% % %here,we use some existed variables in the above Original roughness index computing
% % for i=3:chainCodeLength-2
% %     ERI_r(i)=d(i)-d(i-2);    %adjust the step length to determine which step length is the best.
% % end
% % ERI_r(1)=d(1)-d(chainCodeLength(1)-3);
% % ERI_r(2)=d(2)-d(chainCodeLength(1)-2);
% % ERI_R=real(ifft(ERI_r));     %apply ifft on the adjusted distance difference.
% % %The end of computing Enhance Roughness Index.
% % ERI_R;
% % X = ERI_R;
% % 
% % size(X);
% % X = double(X);
% % X = X';
% % 
% % %k=2;
% % 
% % %figure(3)
% % %subplot(4,1,1)
% % %hist(X,100);
% % %figure(2)
% % [W,M,V,L] = EM_GM(X,k,0.01,1000,1,[]);
% % 
% % 
% % L = MyLikelihood(X,k,W,M,V)
% %     
%  end
% 
% fclose(reportFile);