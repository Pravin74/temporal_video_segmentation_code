addpath('../Data_Loading;..;../Tests;../Features_Extraction;../Utils;../Concept_Detector;../Concept_Detector/fastsmooth;../LSDA;../Evaluation');

%%%%%%%%%%% Parameters %%%%%%%%%%%

%% Data parameters
%%% Features path
data_params.features_path = [pwd '/Features'];
data_params.RC_results_path = [pwd '/Results'];
data_params.RC_plot_results_path = [pwd '/Plot_Results'];

% 
%% CNN parameters (Global Features)
% Installation-dependent
CNN_params.use_gpu = 1;
CNN_params.batch_size = 30; % Depending on the deploy net structure!!
CNN_params.model_file ='/media/hubble/Drive1/temporal_video_segmentation/codes/SR-Clustering/bvlc_reference_caffenet.caffemodel';
CNN_params.size_features = 4096;
CNN_params.caffe_path = '/home/hubble/caffe/matlab/+caffe'; % installation path

% Model-dependent
CNN_params.parallel = false; % allow loading images in parallel
CNN_params.mean_file = '../Utils/ilsvrc_2012_mean.mat';

[structure_path, ~, ~] = fileparts(pwd);
CNN_params.model_def_file = [structure_path '/Utils/deploy_signed_features.prototxt'];

%% Create some folders for results
if(~exist(data_params.features_path, 'dir'))
    mkdir(data_params.features_path);
    mkdir([data_params.features_path '/CNNfeatures']);
    mkdir([data_params.features_path '/SemanticFeatures']);
end
if(~exist(data_params.RC_results_path, 'dir'))
    mkdir(data_params.RC_results_path);
end
if(~exist(data_params.RC_plot_results_path, 'dir'))
    mkdir(data_params.RC_plot_results_path);
end
