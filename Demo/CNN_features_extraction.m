function fMeasureMerge = RClustering(folder, format, data_params,  CNN_params,   GT_path, frame_sampling_rate)

    %% If the optional parameter 'GT_path' is used, the segmentation will be evaluated
    if(nargin < 9)
        eval_GT = false;
        disp('Not performing final evaluation.');
    elseif(~exist(GT_path, 'file'))
        eval_GT = false;
        disp('Not performing final evaluation: non existent GT_file');
    else
        eval_GT = true;
    end
    
    params.features_path = data_params.features_path;

    
    %% Start data processing

    [~, folder_name, ~] = fileparts(folder);
    files = dir([folder '/*' format]);
    files = files(arrayfun(@(x) x.name(1) ~= '.', files));

    %% Check if global features are computed
    path_features = [data_params.features_path '/CNNfeatures/CNNfeatures_' folder_name '.mat']
    
    if(~exist(path_features, 'file'))
        % Compute CNN features
        disp(['Extracting CNN global features of folder ' folder_name]);
        features = extractCNNFeaturesDemo(folder, files, CNN_params);
        this_feat_path = [data_params.features_path '/CNNfeatures'];
        save([this_feat_path '/CNNfeatures_' folder_name '.mat'], 'features','-v7.3');
        clear features;
	tag_matrix = [];
    end
    

    %% Calling ADWIN
    if(eval_GT)
        params.doEvaluation = true;
        %[~,~,cl_limGT, ~]=analizarExcel_Narrative(GT_path, files);
        [cl_limGT]=analizarExcel_Narrative_pravin(GT_path, frame_sampling_rate);

        delim=cl_limGT';
        if delim(1) == 1, delim=delim(2:end); end
        params.GT = delim;

    else
        params.doEvaluation = false;
    end
    params.files_path= folder;
    params.formats = format;
    %%% Get cluster indices applying R-Clustering
    path_here = pwd;
    cd ..
    fMeasureMerge=process_single_sequence_v2(folder, params);
    
      
end
