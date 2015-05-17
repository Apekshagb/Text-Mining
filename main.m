clear all
close all
clc

disp('started')

[M_fulltext,txt_fulltext,raw_fulltext] = xlsread('matrix.xlsx','full_text');
[M_abstract,txt_abstract,raw_abstract] = xlsread('matrix.xlsx','abstract');

nterms = max(max(M_fulltext(:,1)),max(M_abstract(:,1)));
fnames = load('common_filenames__uniq.txt');
ndocs = length(fnames(:,1));
% 
A_fulltext = zeros(nterms,ndocs);
A_abstract = zeros(nterms,ndocs);

% % Mf = M_fulltext;
% % Mf(:,4) = mod(Mfb (:,2),2);
% % Mf(:,2) = floor(Mf(:,2)/2) + 1;
% 

for i = 1:ndocs
    
    % fulltext
    %disp([num2str(fnames(i,1)),'-',num2str(fnames(i,2)),'-',num2str(fnames(i,3)),'-',num2str(fnames(i,4)),'.pdf'])
    filename = [num2str(fnames(i,1)),'-',num2str(fnames(i,2)),'-',num2str(fnames(i,3)),'-',num2str(fnames(i,4)),'.pdf'];
    idx = find(strcmp(raw_fulltext(:,4),filename)==1);
    
    A_fulltext(cell2mat(raw_fulltext(idx,1)),i) = cell2mat(raw_fulltext(idx,3));
end
     

for i = 1:ndocs
    
    % abstract
    %disp([num2str(fnames(i,1)),'-',num2str(fnames(i,2)),'-',num2str(fnames(i,3)),'-',num2str(fnames(i,4)),'ab.pdf'])
    filename = [num2str(fnames(i,1)),'-',num2str(fnames(i,2)),'-',num2str(fnames(i,3)),'-',num2str(fnames(i,4)),'ab.pdf'];
    idx = find(strcmp(raw_abstract(:,4),filename)==1);
    

    A_abstract(cell2mat(raw_abstract(idx,1)),i) = cell2mat(raw_abstract(idx,3));

end

%Fulltext code
[coeff,score,latent] = pca(A_fulltext);
reduced_fulltext = A_fulltext * coeff(:,1:2);
[id,C,sumd] = kmeans(reduced_fulltext,3);

ptsymb = {'bs','r^','md','go','c+'};

wordIndx = 1:114766;
wordIndx = wordIndx';
figure
scatter(reduced_fulltext(:,1),reduced_fulltext(:,2));

figure
for i = 1:3
    clust = find(id==i);
    plot(reduced_fulltext(clust,1),reduced_fulltext(clust,2),ptsymb{i});
%       scatter3(resMat(clust,1),resMat(clust,2),id(clust,1),ptsymb{i});
    hold on
end

figure
for i = 1:3
    clust = find(id==i);
    %plot(wordIndx(clust,1),resMat(clust,1),ptsymb{i});
    %scatter3(resMat(clust,1),resMat(clust,2),resMat(clust,1),ptsymb{i});
    scatter3(reduced_fulltext(clust,1),reduced_fulltext(clust,2),id(clust,1),ptsymb{i});
    hold on
end

%Abstarct Code
[coeff,score,latent] = pca(A_abstract);
reduced_abstract = A_abstract * coeff(:,1:2);
[id1,C1,sumd1] = kmeans(reduced_abstract,3);

ptsymb = {'bs','r^','md','go','c+'};

wordIndx = 1:114766;
wordIndx = wordIndx';
figure
scatter(reduced_abstract(:,1),reduced_abstract(:,2));

figure
for i = 1:3
    clust = find(id1==i);
    plot(reduced_abstract(clust,1),reduced_abstract(clust,2),ptsymb{i});
%       scatter3(resMat(clust,1),resMat(clust,2),id(clust,1),ptsymb{i});
    hold on
end

figure
for i = 1:3
    clust = find(id1==i);
    %plot(wordIndx(clust,1),resMat(clust,1),ptsymb{i});
    %scatter3(resMat(clust,1),resMat(clust,2),resMat(clust,1),ptsymb{i});
    scatter3(reduced_abstract(clust,1),reduced_abstract(clust,2),id1(clust,1),ptsymb{i});
    hold on
end

disp('Done')
