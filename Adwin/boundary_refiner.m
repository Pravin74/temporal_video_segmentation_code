function idx = boundary_refiner(data, fi, p, win_start, idx_lr, idx_rl, win_end, threshold)
%set parameters
epsilon= 1;
% epsilon, 1 - sensitive segmentation; 0 - robust segmentation;
% epsilon = 0;
% INPUT:
%   data   -> n x k data stream (each data in [0, 1]).
%   fi -> segmentation parameter
%   p -> norm parameter
%   flag 0 or 1 -> bound type
%
% OUTPUT:
%   w3    -> n x k data streams with means. 
%   idx   -> idx inside a stream with detected change.

flag = 1;
[~,k]=size(data);
n_tot = win_end - win_start +1;
w=zeros(1, n_tot);
len = min(5,n_tot);
W=data(win_start:win_end,1:k);
pNorm = sum(abs(W).^p,2).^(1/p);
pNorm=pNorm./((k)^(1/p));
variance=var(pNorm);
% Drop elements for the tail of W
% while all splits of W (||mu_w0-mu_w1||)>=ecut)
cut=false;
r1=zeros(1,size(W,1)-1);
r2=zeros(1,size(W,1)-1);
z=cumsum(W);
n_frames=size(W,1)-1;
for i=1:n_frames
    n0=i;
    n1=size(W,1)-i;

    m=(n0*n1)/((sqrt(n0)+sqrt(n1))^2);
    fi_prime=fi/(k*(n0+n1));

    mu_w0=(z(i,:))/i;
    mu_w1=(z(end,:)-z(i,:))/(n_frames+1-i);
    r1(i)=norm(mu_w0 - mu_w1, p);
    r2(i)=threshold;
end
[aa bb]=max(r1-r2);
idx = win_start+bb;
end


