clear all, close all
warning off

load YaleBCrop025.mat

nSet = [2 3 5 8 10];
for i = 1:length(nSet)	
    n = nSet(i);
    idx = Ind{n};   
    for j = 1:size(idx,1)		
        X = [];
		for p = 1:n
			X = [X Y(:,:,idx(j,p))];
		end
		X = X/max(X(:));		
        [D,N] = size(X);                  
		
		K = max(s{n});
		lambda(1) = 6e-3;%5e-2;
		lambda(2) = 3e-3;%5e-4;			
		affine = false; outlier = true; r = 9*K; % in the original paper, we used r = 10*K
        post = true;				
		
		missrate = nullSpaceClustering(X,s{n},lambda,r,outlier,post,affine,9);
						
        missrateTot{n}(j) = missrate;		
		
		disp([num2str(n) ' subjects, ' 'sequence ' num2str(j) ': ' num2str(100*missrateTot{n}(j)) '%']);					
	end
	
    avgmissrate(n) = mean(missrateTot{n});
    medmissrate(n) = median(missrateTot{n});	
	disp([num2str(n) ' subjects: ']);
	disp(['Mean: ' num2str(100*avgmissrate(n)) '%, ' 'Median: ' num2str(100*medmissrate(n)) '%']);
end
save('NSC_Faces.mat', 'missrateTot', 'avgmissrate', 'medmissrate')