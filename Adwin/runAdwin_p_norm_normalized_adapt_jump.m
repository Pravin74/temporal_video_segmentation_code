function [fMeasureMerge] = runAdwin_p_norm_normalized_adapt_jump(data, fi, p, params, folder_name, jump)
%set parameters
% load('/home/hubble/Desktop/workspace_yair.mat')
% data= data(1:5000,:);
%
% data Normalization
%data = normalize(data);
[n_tot, dimension]=size(data);
disp(['Total Frames : ' num2str(n_tot)]);
fi
p
flag_want_skip = 1  %1 if you want to use skip factor otherwise 0
if flag_want_skip == 1
    corr_coef_th = .95
else
    corr_coef_th = 0;
    disp ('No skip used')
end
tolerance = 600 ;

% Normalizing using pth-norm 
p_norms = vecnorm(data',p);
data_new = data'./(p_norms);
data = data_new';

%empirically testing variacne 
% load('/home/enigma/Desktop/result_HUJI_Yair_4fps_Sampling_1_Corr_coef_threshold_0.8_Number_of_segments_43_adapt_jump_mansi_GT.mat')
% GT = find(GT_one_hot==1);
% GT = [1;GT]
% all_seg_var = []
% for ii = 2: size(GT,1)
%     data_temp = data(GT(ii-1):GT(ii),:);
%     temp_var = trace(data_temp*data_temp')/size(data_temp, 1);
%     all_seg_var = [all_seg_var,temp_var]
% end
% empirically testing of variance
% all_seg_var = [];
% for ii = 1: n_tot
%     data_temp = data(ii,:);
%     temp_var = trace(data_temp*data_temp')/size(data_temp, 1);
%     all_seg_var = [all_seg_var,temp_var];
% end

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
for t=len+1:n_tot-1
    W=[W; data(t,:)] ;
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
            skipp = 1;
        end
        W1_sum = zeros(dimension+1, dimension+1);
        for i=1:n_frames
            n0 = i;
            n1 = size(W,1)-i;
            %% calculations for variance
            %temp = W(i,:);
            %pilation_i = [temp*temp', zeros(1,dimension); zeros(dimension,1), temp'*temp];
            %pilation_2 = [temp'*temp, zeros(dimension,1); zeros(1, dimension), temp*temp'];
            %pilation = [zeros(dimension,dimension), temp'; temp, 0];
            %sqrt(max(eig(pilation'*pilation)))
            %% calcuations for E_cut
            m_by_2  = ((n0*n1)/ (n0+n1))/skipp;
            fi_prime = fi/((n0+n1)/skipp);
            % Assuming worst case variance
            sigma = 1/4;
            outer_term = (8* sigma)/(m_by_2) ;
            ecut2 = (outer_term * log((2*(dimension+1))/(fi_prime)))^(1/2);
            r2(i) = ecut2 ;
            % calculation for mean
            mu_w0 = (z(i,:))/i;
            mu_w1 = (z(end,:)-z(i,:))/(n_frames+1-i);
            %diff_mean = abs(mu_w0- mu_w1);
            diff_mean = mu_w0- mu_w1;
            r1(i) = sqrt(eig(diff_mean* diff_mean'));
            
        end
        [aa bb]=max(r1-r2);
        cut_threshold = [cut_threshold, r2(bb)];
        mean_diff = [mean_diff, r1(bb)];
        if(max(r1-r2)<0)
            cut=true;
        end
        if(cut==false)
            win_start = [win_start, sum(temp_start)]
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
            disp(['Frames processed: ' num2str(t)]);
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
file_name = strcat(params.RC_results_path, '/result_', folder_name, '_delta_' ,num2str(fi) ,'_Corr_coef_threshold_' ,num2str(corr_coef_th), '_Number_of_segments_', num2str(length(indexes)), '.png');
saveas(fig, file_name)

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
% segmentation = {};
% segmentation_file = [params.RC_results_path '/result_' folder_name '_delta_' num2str(fi) '_Corr_coef_threshold_' num2str(corr_coef_th) '_Number_of_segments_' num2str(length(indexes)) '.csv'];
% f = fopen(segmentation_file, 'w');

[final_boundaries]=compute_boundaries(labels_modified, num_frames);
num_clusters = length(final_boundaries)+1;
save([file_name(1:end-4) '_predicted.mat'], 'final_boundaries');


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

