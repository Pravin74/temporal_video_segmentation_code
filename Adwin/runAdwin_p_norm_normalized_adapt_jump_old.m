function [fMeasureMerge] = runAdwin_p_norm_normalized_adapt_jump(data, fi, p, params, folder_name, jump, files)
%set parameters
% load('/home/hubble/Desktop/workspace_yair.mat')
% data= data(1:5000,:);
epsilon= 1; % epsilon, 1 - sensitive segmentation; 0 - robust segmentation;
%
[n_tot, ~]=size(data);
disp(['Total Frames : ' num2str(n_tot)]);
fi
p
flag =0;
flag_want_skip = 0 %1 if you want to use skip factor otherwise 0
if flag_want_skip == 1
    corr_coef_th = .94
else
    corr_coef_th = 0;
    disp ('No skip used')
end

tolerance = 600 ;
% Normalizing dataset
% data = (data - min(data(:))) / ( max(data(:)) - min(data(:)));
% Normalizing using pnd-norm 
p_norms = vecnorm(data',p);
data_new = data'./(p_norms);
data = data_new';

%% calculating the correlation coeff. 
i=1;
jump_adap = [];
jump_idx = [];
while (1)    
    for j=i:n_tot
        corr = corrcoef(data(i,:),data(j,:));
        if corr(1,2)<corr_coef_th
            break;
        end
    end
    jump_adap= [jump_adap ,j-i];
    jump_idx = [jump_idx, i];
%     tt=imread(strcat(path,num2str(i,'%06.f'), '.jpg'));
%     tt_2 =imread(strcat(path,num2str(i+(j-i),'%06.f'), '.jpg')) ;
%     imshowpair(tt, tt_2 , 'montage');
%     w = waitforbuttonpress;
%     axes;
    i=j;
    if i==n_tot
        break;
    end
end 
%%
% k =1 for original ADWIN
% k = 1;
%w3=zeros(n_tot, k);
w=zeros(1, n_tot);
len = min(50,n_tot);
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
confidence_n = [];
confidence_sqrt_n = [];
for t=len+1:n_tot-1
    W=[W; data(t,:)] ;
    % Drop elements for the tail of W
    % while all splits of W (||mu_w0-mu_w1||)>=ecut)
    cut=false;
    while(cut==false)
        r1 = zeros(1,size(W,1)-1);
        r2 = zeros(1,size(W,1)-1);
     
        z=cumsum(W);
        n_frames=size(W,1)-1;
        % code to calculate adaptive skip factor
        if corr_coef_th ~=0
            skipp = adaptive_independence(W, corr_coef_th);
            if isnan(skipp)
                skipp=  mean(jump_adap(find(jump_idx<=t)));
            end
        else
            skipp =1;
        end
        for i=1:n_frames
            n0 = i;
            n1 = size(W,1)-i;
            m = ((n0*n1)/ (n0+n1))/skipp;
            fi_prime = fi/((n0+n1)/skipp);
            %ecut2 = ((2/m) * log((4)/fi_prime))^(1/2);
            ecut2 =((2/m) * log((4)/fi_prime))^(1/2);
            %ecut2 = ecut2*1.5;
            mu_w0=(z(i,:))/i;
            mu_w1=(z(end,:)-z(i,:))/(n_frames+1-i);
            %% code for signed norm.
%             mu_w0_unit = mu_w0/sqrt(sum(mu_w0.^2));
%             mu_w1_unit = mu_w1/sqrt(sum(mu_w1.^2));
%             cos_theta = sum(mu_w0_unit.*mu_w1_unit);
            
%           cos_theta = sum(mu_w0.*mu_w1)/ (sqrt(sum(mu_w0.^2))* sqrt(sum(mu_w1.^2)));
%           sin_theta = sqrt(1-cos_theta^2);
%             theta = (acos(cos_theta));

            %r1(i) = (norm(mu_w0-mu_w1,p)) * (1/2);
            r1(i) = norm(mu_w0-mu_w1,p); 
            %r1_first_norm(i) = (norm(mu_w0-mu_w1,1)); 
            r2(i) = ecut2 ;
%             r3_mean_diff(i) = norm(mu_w0-mu_w1,p);
%             r4_theta(i) = theta;
        end
        [aa bb]=max(r1-r2);
        cut_threshold = [cut_threshold, r2(bb)];
        mean_diff = [mean_diff, r1(bb)];
%         theta_mean_diff = [theta_mean_diff, r1(bb)];
%         theta_only = [theta_only, r4_theta(bb)];
%         % calculate the variance of right window
%         W_right=W(bb+1:end,:);
%         pNorm = sum(abs(W_right).^p,2).^(1/p);
%         variance=var(pNorm);
%         variances_t = [variances_t variance];

        if(max(r1-r2)<0)
            cut=true;
        end
        if(cut==false)
            %disp(ecut2)
            %disp(max(r1))
            %disp(max(r1-r2))
            win_start = [win_start, sum(temp_start)]
            %confidence_n = [confidence_n, fi/(n0+n1)];
            %confidence_sqrt_n = [confidence_sqrt_n, fi/sqrt(n0+n1)];
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
     if mod(t,1000)==0
            disp(['Frames: ' num2str(t)]);
     end
end

x=w(2:end)-w(1:end-1);
indexes=abs(x(find(x<0)))+1;
indexes=cumsum(indexes);

num_frames = size(data,1) ;
indexes = [indexes, num_frames];

figure;
fig = plot([cut_threshold; mean_diff; cut_threshold-mean_diff]');
fig = legend(fig, {'\epsilon_{cut}','|W_0 - W_1|', '|W_0 - W_1|-\epsilon_{cut}'},'Location','southeast');
saveas(fig, strcat(params.RC_results_path, '/result_', folder_name, '_delta_', num2str(fi) , '_Corr_coef_threshold_' ,num2str(corr_coef_th), '_Number_of_segments_', num2str(length(indexes)), '.png'))

labels_modified = zeros(1,num_frames);
label = 1;
j=1;
for i=1:num_frames
    if i<=indexes(j)
      labels_modified(i)=label;
    else
      label = label + 1;
      j=j+1;
      labels_modified(i)=label;
    end
end

events = labels_modified;
num_clusters = max(events);
nFrames = length(events);
result_data = cell(1, num_clusters);
for i = 1:nFrames
    if(events(i) ~= 0)
        result_data{events(i)} = [result_data{events(i)} i];
    end
end
segmentation = {};
segmentation_file = [params.RC_results_path '/result_' folder_name '_delta_' num2str(fi) '_Corr_coef_threshold_' num2str(corr_coef_th) '_Number_of_segments_' num2str(length(indexes)) '.csv'];
f = fopen(segmentation_file, 'w');

[final_boundaries]=compute_boundaries(labels_modified, num_frames);
num_clusters = length(final_boundaries)+1;
save([segmentation_file(1:end-4) '_predicted.mat'], 'final_boundaries');

% final_boundaries_with_confidence = [final_boundaries(1:end-1); confidence_n; confidence_sqrt_n];
% save([segmentation_file(1:end-4) '_predicted_segments_with_confidence.mat'], 'final_boundaries_with_confidence');

for s = 1:num_clusters
    nImgs = length(result_data{s});
    line = ['Segment_' num2str(s)];
    for i = 1:nImgs
        img_name = files(result_data{s}(i)).name;
        segmentation{s}{i} = img_name;
%             segmentation{s} = {segmentation{s}, img_name};
        line = [line ',' img_name];
    end
    fprintf(f, [line '\n']);
end

fclose(f);
disp(['A .csv file with the result has been stored in ' params.RC_results_path]); 
 
doEvaluation = params.doEvaluation;
if(doEvaluation)
    GT = params.GT;
else
    GT = [];
end

disp('-------- Results small segments merging of adapt_jump--------');

if(doEvaluation)
    [recMerge,precMerge,accMerge,fMeasureMerge] = Rec_Pre_Acc_Evaluation(GT,final_boundaries,num_frames,tolerance);

    disp(['Precision: ' num2str(precMerge)]);
    disp(['Recall: ' num2str(recMerge)]);
    disp(['F-Measure: ' num2str(fMeasureMerge)]);
    GT_one_hot = zeros(num_frames,1);
    GT_one_hot(GT)=1;
    save([segmentation_file(1:end-4) '_GT.mat'], 'GT_one_hot');

else
    fMeasureMerge= 0;
end

disp(['Number of events: ' num2str(num_clusters)]);
disp(['Mean frames per event: ' num2str(num_frames/num_clusters)]);
disp(' ');


end

