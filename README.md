# Concept Drift Detection forMultivariate Data Streams and Temporal Segmentation ofDaylong Egocentric Videos

This is the implementation of our paper titled "Concept Drift Detection for Multivariate Data Streams and Temporal Segmentation of Daylong Egocentric Videos"  in ACMMM2020. We present a novel unsupervised temporal segmentation technique especially suited for day-long egocentric videos.



## Requirements
The code has been tested on:
 - Caffe Deep Learning Framework and matcaffe wrapper (for global features calculation)  
	       Caffe main page: http://caffe.berkeleyvision.org/  
		     Good Linux installation tutorial: https://github.com/tiangolo/caffe/blob/ubuntu-tutorial-b/docs/install_apt2.md  
		     CaffeNet model: http://dl.caffe.berkeleyvision.org/bvlc_googlenet.caffemodel  
- Compile files in GCMex for your system.  
- MATLAB 2016 (higher version also works)  
- Python 2.7  


## Get started
To reproduce the results:
- Download the features from project page in './Demo/Features/CNNfeatures' folder.
- Change the parameters '\delta' ('confidence' variable in 'process_single_sequence_v2.m' file) and '\rho_c' ('corr_coef_th' in 'runAdwin_p_norm_normalized_adapt_jump.m' file) as mentioned in supplementary material.
- Now run './Demo/demo.m'
- Find the output file in the './Demo/Results' folder.

To run on a new video sequence:
- Extract the images and give the folder path in 'demo.m'.
- Adjust '\delta' and '\rho_c' and execute the 'demo.m'.
- Find the output file in the './Demo/Results' folder.


## Citation
```
@inproceedings{nagar2020concept,
  title={Concept Drift Detection for Multivariate Data Streams and Temporal Segmentation of Daylong Egocentric Videos},
  author={Nagar, Pravin and Khemka, Mansi and Arora, Chetan},
  booktitle={Proceedings of the 28th ACM International Conference on Multimedia},
  pages={1065--1074},
  year={2020}
}
```
