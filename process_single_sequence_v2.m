function [jump] = process_single_sequence_v2(folder, params)

    %% Load paths
    addpath('Adwin;Data_Loading;Evaluation;Features_Preprocessing');
    addpath('GCMex;GraphCuts;PCA;Tests;Utils;SpectralClust');
    load_features = true;

    %% Parameters loading
    fichero = params.files_path;               
    formats = params.formats;
    doEvaluation = params.doEvaluation;
        if(doEvaluation)
            GT = params.GT;
        else
            GT = [];
        end

   
    paramsfeatures.type = 'CNN'; %CNN ....
    paramsPCA.minVarPCA=0.95;
    % paramsPCA.minVarPCA=0.8;
    paramsPCA.standarizePCA=false;
    paramsPCA.usePCA_Clustering = true;
    
    plotFigResults = false;
    %% Adwin parameters
    pnorm = 2;
    confidence =  0.000001;
    paramsPCA.usePCA_Adwin = true;

 

    %% Build paths for images, features and results
    [~, folder_name, ~] = fileparts(folder);
    path_features = [params.features_path '/CNNfeatures/CNNfeatures_' folder_name '.mat'];
    path_features_PCA = [params.features_path '/CNNfeatures/CNNfeaturesPCA_' folder_name '.mat'];


     %% Images
    files_aux=dir([fichero '/*' formats]);
    count = 1;
	files = struct('name', []);
    for n_files = 1:length(files_aux)
        if(files_aux(n_files).name(1) ~= '.')
            files(count).name = files_aux(n_files).name;
            count = count+1;
        end
    end
    Nframes=length(files);
    



    %% Global Features
    if strcmp(paramsfeatures.type, 'CNN')
        if(load_features)
            load(path_features);
            % features = features(1:10:end,:);
            [features_norm] = signedRootNormalization(features);
        end

	if(size(features,1) ~= Nframes)
		error('The number of Global features does not match the number of images. TIP: remove the existent features file for re-calculation.');
	end

        %PCA FEATURES
%         if(exist(path_features_PCA) > 0)
%             load(path_features_PCA);
%         else
            [ featuresPCA, ~, ~ ] = applyPCA( features_norm, paramsPCA ) ; 
            if(load_features) % if we wanted to load the stored features, then we will also store PCA features
                save(path_features_PCA, 'featuresPCA');
            end
%         end
    end
    
   
    %% Check if we only have one sample
    if(paramsPCA.usePCA_Adwin && strcmp(paramsfeatures.type, 'CNN'))
        num_samp = size(featuresPCA,1);
    elseif( strcmp(paramsfeatures.type, 'CNN'))
        num_samp = size(features_norm,1);
    end
    if(num_samp == 1)
        events = [1];

    else

% picking a vidoe of 35 minutes. 
 %featuresPCA = featuresPCA(1:end,:);
 %% code for GT and jump
jump = 1;



featuresPCA = featuresPCA(1:jump:end,:);
 
%% ADWIN
 disp(['Start ADWIN ' folder_name]);
        [fMeasureMerge] = runAdwin_p_norm_normalized_adapt_jump(featuresPCA , confidence, pnorm, params, folder_name);

end
