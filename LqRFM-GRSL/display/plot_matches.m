function plot_matches(I1, I2, X, Y, VFCIndex, CorrectIndex)
%   PLOT_MATCHES(I1, I2, X, Y, VFCINDEX, CORRECTINDEX)
%   only plots the ture positives with blue lines, false positives with red
%   lines, and false negatives with green lines. For visibility, it plots at
%   most NUMPLOT (Default value is 50) matches proportionately.
%   
% Input:
%   I1, I2: Tow input images.
%
%   X, Y: Coordinates of intrest points in I1, I2 respectively.
%
%   VFCIndex: Indexes preserved by VFC.
%
%   CorrectIndex: Ground truth indexes.
%
%   See also:: VFC(), FastVFC(), SparseVFC().

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    04/17/2012

% Define the maximum number of matches to plot
NumPlot = 50000;

TruePos = intersect(VFCIndex, CorrectIndex);%Ture positive
FalsePos = setdiff(VFCIndex, CorrectIndex); %False positive
FalseNeg = setdiff(CorrectIndex, VFCIndex); %False negative

NumPos = length(TruePos)+length(FalsePos)+length(FalseNeg);
if NumPos > NumPlot
    t_p = length(TruePos)/NumPos;
    n1 = round(t_p*NumPlot);
    f_p = length(FalsePos)/NumPos;
    n2 = ceil(f_p*NumPlot);
    f_n = length(FalseNeg)/NumPos;
    n3 = ceil(f_n*NumPlot);
else
    n1 = length(TruePos);
    n2 = length(FalsePos);
    n3 = length(FalseNeg);
end

per = randperm(length(TruePos));
TruePos = TruePos(per(1:n1));
per = randperm(length(FalsePos));
FalsePos = FalsePos(per(1:n2));
per = randperm(length(FalseNeg));
FalseNeg = FalseNeg(per(1:n3));

im3 = appendimages(I1,I2);
figure
imshow(im3,[]);
% colormap('gray');
% imagesc(im3);
% interval = 20;
% WhiteInterval = 255*ones(size(I1,1), interval, 3);
% figure;imagesc(cat(2, I1, WhiteInterval, I2)) ;
hold on ;
line([X(FalsePos,1)'; Y(FalsePos,1)'+size(I1,2)], [X(FalsePos,2)' ;  Y(FalsePos,2)'],'linewidth', 1.5, 'color', 'r') ;
% line([X(FalseNeg,1)'; Y(FalseNeg,1)'+size(I1,2)+interval], [X(FalseNeg,2)' ;  Y(FalseNeg,2)'],'linewidth', 1, 'color', 'g') ;
line([X(TruePos,1)'; Y(TruePos,1)'+size(I1,2)], [X(TruePos,2)' ;  Y(TruePos,2)'],'linewidth', 1.5, 'color', 'b') ;
axis equal ;axis off  ; 
drawnow;