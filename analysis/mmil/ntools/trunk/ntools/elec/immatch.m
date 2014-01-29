clear all; close all;
[input_img input_img_path] = uigetfile({'*.jpg';'*.png';'*.bmp';'*.*'},'Select the target image');
[base_img base_img_path] = uigetfile({'*.jpg';'*.png';'*.bmp';'*.*'},'Select the template image',input_img_path);
input = imread([input_img_path input_img]);
base = imread([base_img_path base_img]);
[input_points,base_points] = cpselect(input,base,'Wait',true);
input_points_corr = cpcorr(input_points, base_points,input(:,:,1),base(:,:,1));
t_concord = cp2tform(input_points_corr,base_points,'lwm');
input_trans = imtransform(input,t_concord,'bicubic');
imshow(input_trans,[])