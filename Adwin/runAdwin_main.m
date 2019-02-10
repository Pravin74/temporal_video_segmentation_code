function [indexes, win_start, win_end, max_mean] = runAdwin_main(data, fi, p)
%set parameters
epsilon= 1; % epsilon, 1 - sensitive segmentation; 0 - robust segmentation;
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
% data = data(:,1);
% data = repmat(data, [1,84]) ;
%% code to make the video data IID.



%%
flag = 1;
[n_tot, ~]=size(data);
% Normalizing dataset
data = (data - min(data(:))) / ( max(data(:)) - min(data(:)));
% calulating norm to convert form x_t to L_t
%data = vecnorm(data,p,2);

% k =1 for original ADWIN
% k = 1;
%w3=zeros(n_tot, k);

w=zeros(1, n_tot);

len = min(5,n_tot);

if len == 1
    indexes = [1];
else
k=1;
W=data(1:len,:);
win_start = [];
win_end = [];
max_mean = [];
temp_start = [1];
cut_threshold = [];
mean_diff = [];
index_to_pick = 0;
for t=len+1:n_tot
    W=[W; data(t,:)] ;
    pNorm = sum(abs(W).^p,2).^(1/p);
    variance=var(pNorm);
    
    % Drop elements for the tail of W
    % while all splits of W (||mu_w0-mu_w1||)>=ecut)
    cut=false;
    
    while(cut==false)
        
        r1=zeros(1,size(W,1)-1);
        r2=zeros(1,size(W,1)-1);
        
        z=cumsum(W);
        n_frames=size(W,1)-1;
        
        for i=1:n_frames
            n0=i;
            n1=size(W,1)-i;
      
            m=(n0*n1)/((sqrt(n0)+sqrt(n1))^2);
            m=(n0*n1)/ (n0+n1);
            fi_prime=fi/(k*(n0+n1));
            fi_prime=fi/((n0+n1));
                                             
            if flag ==1
                ecut2=sqrt( (2/m) * variance * log(2/fi_prime) ) + (2/(3*m))*log(2/fi_prime);
            elseif flag==0             
                ecut2=((1/(2*m)) * log((4)/fi_prime) )^(1/2);
            end

            mu_w0=(z(i,:))/i;
            mu_w1=(z(end,:)-z(i,:))/(n_frames+1-i);
            r1(i)=norm(mu_w0 - mu_w1, p);
            r2(i)=ecut2;
        end
         [aa bb]=max(r1-r2);
         cut_threshold = [cut_threshold, r2(bb)];
         mean_diff = [mean_diff, r1(bb)];
        if(max(r1-r2)<0)
            cut=true;
        end
        if(cut==false)
            %disp(ecut2)
            %disp(max(r1))
            %disp(max(r1-r2))
            win_start = [win_start, sum(temp_start)];
            [aa bb]=max(r1-r2);
            win_end = [win_end, t];
            temp_start = [temp_start, bb];
            max_mean = [max_mean, max(r1)];
            W_tmp=W(bb+1:end,:);
            clear W;
            W=W_tmp;
            clear W_tmp;
           
        end
    end
    w(t)=size(W,1);
    
   
end

figure;
plot([cut_threshold;mean_diff; mean_diff-cut_threshold]');
legend({'\epsilon_{cut}','|W_0 - W_1|', '|W_0 - W_1|-\epsilon_{cut}'},'Location','southeast')

x=w(2:end)-w(1:end-1);
indexes=abs(x(find(x<0)))+1;
indexes=cumsum(indexes);
%indexes=[1 indexes length(w)];

end


