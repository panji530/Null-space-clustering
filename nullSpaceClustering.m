function [Missrate, C] = nullSpaceClustering(X,s,lambda,r,outlier,post,affine,Dim)
% Pan Ji, pan.ji@anu.edu.au
if(nargin<8)
	Dim = 4;
end
if(nargin<7)
	affine = false;
end
if(nargin<6)
	post = true;	
end
if(nargin<5)
	outlier = false;
end
if(nargin<4)
	r = 0;
end
if(nargin<3)
	if(~outlier)	
		lambda = 240;		
	end	
elseif(length(lambda) == 2)
	lambda1 = lambda(1);
	lambda2 = lambda(2);
else
	lambda1 = lambda;
	lambda2 = lambda/2;
end

K = max(s);
X = dataProjection(X,r);
[d,N] = size(X);
tic
if(~outlier)
	if(~affine)
		C = (eye(N)+lambda*(X'*X))\eye(N);
	else
		nu = null(ones(1,N));
		lhd = (nu'*nu)+lambda*nu'*(X'*X)*nu;
		rhd = nu';
		phi = lhd\rhd;
		C = nu*phi;
	end		
else
	%initialization of parameters
	rho = 1e-6;
	max_rho = 1e10;
	eta = 3.5;
	epsilon = 1e-8;
	maxIter = 1e4;

	Y = zeros(d,N);
	C = zeros(N,N);
	E = zeros(d,N);
	
	xtx = X'*X;
	%begin iterations
	iter = 0;
	while(true)
		iter = iter+1;
		%update C
		lhd = eye(N) + (lambda1+rho)*xtx;
		rhd = eye(N) - X'*Y + rho*X'*E;
		C = lhd\rhd;
		%update E
		temp = X*C + Y/rho;
		E = max(temp - lambda2/rho,0) + min(temp + lambda2/rho,0);
		
		leq = X*C - E;
		stpC = max(max(abs(leq)));
		if(iter == 1 || mod(iter,50)==0 || stpC<epsilon)
				%disp(['iter ' num2str(iter) ',rho=' num2str(rho,'%2.1e') ',stopALM=' num2str(stpC,'%2.3e')]);
		end
		if(stpC<epsilon || iter>maxIter)
			break;
		else
			Y = Y + rho*leq;
			rho = min(max_rho, rho*eta);
		end		
    end
end

C = C - diag(diag(C));
if(~post)	
	W = abs(C) + abs(C');
	grp = ncutW(W,K); % Install your ncutW function from https://www.cis.upenn.edu/~jshi/software/.
	tmp = zeros(N,1);
    for i = 1:K
        tmp = tmp+grp(:,i)*i;
    end
    grp = bestMap(s,tmp);
    Missrate = sum(s(:) ~= grp(:)) / length(s);    
else
	grp = postProC(C,K,Dim);
	grp = bestMap(s,grp);
    Missrate = sum(s(:) ~= grp(:)) / length(s);
end



















