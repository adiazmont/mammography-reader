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
%X %%ADDED LINE 09/18/15 FOR GETTING RAW DATA
size(X);

k=2;

%% CIRCLE CLASS DATA
W = [0.58527 0.41473];
M = [4.34332E-05 3.91458E-05];
V(:,:,1) = 1.26E-10;
V(:,:,2) = 1.82E-09;

%% RECTANGLE CLASS DATA
% W = [0.44964 0.55036];
% M = [-0.014910245 -0.000774341];
% V(:,:,1) = 5.68E-10;
% V(:,:,2) = 1.95E-11;

L = Likelihood(X,k,W,M,V);
L
end
