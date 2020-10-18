%Temporal segmentation of Day long egocentric videos
% start matlab by following command
%LD_PRELOAD="/usr/lib64/libstdc++.so.6" sh matlab

%
clc;
clear;
% Add caffe/matlab to your Matlab search PATH in order to use matcaffe
if exist('/home/hubble/caffe/matlab/+caffe', 'dir')
  addpath('/home/hubble/caffe/matlab/'); 
else
  error('Please run this demo from caffe/matlab/demo');
end

% Load parameters
loadParametersDemo;

% Folder with lifelogging images that we want to segment
%folder = [pwd '/test_data/Subject1'];

%folder = '/media/hubble/Drive1/temporal_video_segmentation/datasets/Disney/pngs/Alireza_Day_2';
% folder = '/media/hubble/Drive1/temporal_video_segmentation/datasets/Disney/pngs/Alireza_Day_2_4fps';
folder_path='/media/hubble/Drive1/temporal_video_segmentation/datasets/Disney/pngs/2fpm/';
folder_path='/media/hubble/Drive1/temporal_video_segmentation/datasets/Disney/pngs/5fps/';


% Format of the images in 'folder
images_format = '.png';

% .csv, .txt or .xls file with the GT segmentation of the data
% (OPTIONAL, set to empty string if you don't want an evaluation of the result) 
GT_path = '/media/hubble/Drive1/temporal_video_segmentation/datasets/EDUB-Seg/GT/';
GT_path = '/media/hubble/Drive1/temporal_video_segmentation/datasets/Disney/annotations/Michael_Day_2.xls';
GT_path = '';
tic


%folders = folders(3:end);
GT='';

videos_all = dir([folder_path '/*5fps']);
videos = videos_all(1:1);
%% The R-Clustering segmentation is applied
f1_score_all = [];
frame_sampling_rate=1/30;
for i =1:1
    tic
    disp(['Start processing: ', videos(i).name])
    [folder_path videos(i).name]

    fMeasureMerge = CNN_features_extraction([folder_path videos(i).name], images_format, data_params, CNN_params,  GT_path, frame_sampling_rate);
    cd '/media/hubble/Drive1/temporal_video_segmentation/codes/SR-Clustering_p_norm_normalized_Disney/Demo/'

    f1_score_all = [f1_score_all fMeasureMerge];
    close all;
    toc
    disp('Total time (Feature extraction and ADWIN)');
end
% disp(['Average F1- Score for Disney-Seg: ' num2str(mean(f1_score_all))]);
% disp('Total time (Feature extraction and ADWIN)');
disp('FINISH');

%%





