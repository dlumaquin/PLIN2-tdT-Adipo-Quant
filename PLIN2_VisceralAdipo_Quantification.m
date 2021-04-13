%% Image Segmentation of the Visceral Adipose Tissue in 21 dpf ubb:plin2-tdtomato zebrafish
% 2020, Dianne Lumaquin, Richard M. White Lab at MSKCC
% for any questions, contact Dianne Lumaquin at dil2009@med.cornell.edu

%% Pipeline requires Bio-Formats
% download at https://docs.openmicroscopy.org/bio-formats/6.1.0/users/matlab/index.html

%% Set up

clear all
close all

analyze_flag = 1; % run analysis on images

% code demonstrates how to conduct analysis on 4 groups
% modify code as needed depending on number of experimental conditions

crop_wt_flag = 1; % run crop on WT images
load_wtcrop_flag = 0; % load crop WT positions
crop_mut_flag = 1; % run crop on mut images
load_mutcrop_flag = 0; % load crop mut positions
crop_mut_b_flag = 1; % run crop on mut_b
load_mutcrop_b_flag = 0; % load crop mut_b positions
crop_mut_c_flag = 1;% run crop on mut_c
load_mutcrop_c_flag = 0;% load crop mut_c positions

load_analysis_flag = 0; % load analysis .mat files

showIm = 1;
show_plots = 1;
split_channels = 1; % bioformats plugin: save data for each channel as .mat file
save_server = 1;
chName = {'BF', 'mCh1', 'tdT1', 'tdT2', 'GFP1', 'GFP2'};

% list channels
BF = find(strcmp(chName, 'BF'));
mCh1 = find(strcmp(chName, 'mCh1'));
tdT1 = find(strcmp(chName, 'tdT1'));
tdT2 = find(strcmp(chName, 'tdT2'));
GFP1 = find(strcmp(chName, 'GFP1'));
GFP2 = find(strcmp(chName, 'GFP2'));


folder_wt = ('INSERT FOLDER PATH'); 
folder_mut = ('INSERT FOLDER PATH');  
folder_mut_b = ('INSERT FOLDER PATH'); 
folder_mut_c = ('INSERT FOLDER PATH');

Analysisoutput = ('INSERT FOLDER PATH');

%% %% BIOFORMATS .CZI IMPORT - WT

if analysis_wt_flag
    
    fprintf('\nAnalyzing WT images...\n');
    
    filenames_wt = fullfile(folder_wt, '*.czi');
    CZI_files_wt = dir(filenames_wt);
    CZI_file_paths_wt = fullfile(folder_wt, {CZI_files_wt.name});
    
    
if split_channels
        
        for ii = 1:length(CZI_files_wt);
            reader = bfGetReader(CZI_file_paths_wt{1,ii});
            fprintf('Importing WT image %d/%d\n', ii, length(CZI_files_wt));
            
            BF_wt{1,ii} = bfGetPlane(reader,1);
            mCh1_wt{1,ii} = bfGetPlane(reader,2);
            tdT1_wt{1,ii} = bfGetPlane(reader,3);
            tdT2_wt{1,ii} = bfGetPlane(reader,4);
            GFP1_wt{1,ii} = bfGetPlane(reader,5);
            GFP2_wt{1,ii} = bfGetPlane(reader,6);
        end
        
    end
    
    % save to same folder on server where data is from:
    if save_server
        
        path_save_server1 = strcat(fullfile(folder_wt, 'wt_channel1.mat'));
        save([path_save_server1], 'BF_wt')
        path_save_server2 = strcat(fullfile(folder_wt, 'wt_channel2.mat'));
        save([path_save_server2], 'mCh1_wt')
        path_save_server3 = strcat(fullfile(folder_wt, 'wt_channel3.mat'));
        save([path_save_server3], 'tdT1_wt')
        path_save_server4 = strcat(fullfile(folder_wt, 'wt_channel4.mat'));
        save([path_save_server4], 'tdT2_wt')
        path_save_server5 = strcat(fullfile(folder_wt, 'wt_channel5.mat'));
        save([path_save_server5], 'GFP1_wt')
        path_save_server6 = strcat(fullfile(folder_wt, 'wt_channel6.mat'));
        save([path_save_server6], 'GFP2_wt')
    end
    
    % metadata
    reader = bfGetReader(CZI_file_paths_wt{1,1});
    omeMeta = reader.getMetadataStore();
    imageSizeX = omeMeta.getPixelsSizeX(0).getValue();
    imageSizeY = omeMeta.getPixelsSizeY(0).getValue();
    pixelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value();
    pixelSizeX_num = pixelSizeX.doubleValue();
    pixelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value();
    pixelSizeY_num = pixelSizeY.doubleValue();
    pixelSize_um2 = pixelSizeX_num * pixelSizeY_num;
    
    % save metadata
    metadata_wt = [pixelSize_um2 pixelSizeX_num pixelSizeY_num];
    path_save_server_wt = strcat(fullfile(folder_wt, 'metadata_wt.mat'));
    save([path_save_server_wt], 'metadata_wt');
    fprintf('wt metadata .mat saved')
    
    
end

%% BIOFORMATS .CZI IMPORT - MUTANT

if analysis_mut_flag
    
    % mutant
    fprintf('\nAnalyzing mutant images...\n');
    
    filenames_mut = fullfile(folder_mut, '*.czi');
    CZI_files_mut = dir(filenames_mut);
    CZI_file_paths_mut = fullfile(folder_mut, {CZI_files_mut.name});
    
    if split_channels
            
            % RUN TO SAVE ALL IMAGES/CHANNEL AS 1 .MAT:
            
            for ii = 1:length(CZI_files_mut);
                reader = bfGetReader(CZI_file_paths_mut{1,ii});
                fprintf('Importing Mut image %d/%d\n', ii, length(CZI_files_mut));
                
                BF_mut{1,ii} = bfGetPlane(reader,1);
                mCh1_mut{1,ii} = bfGetPlane(reader,2);
                tdT1_mut{1,ii} = bfGetPlane(reader,3);
                tdT2_mut{1,ii} = bfGetPlane(reader,4);
                GFP1_mut{1,ii} = bfGetPlane(reader,5);
                GFP2_mut{1,ii} = bfGetPlane(reader,6);

            end
            
    end
            
            % save to same folder on server where data is from:
            if save_server
                
                path_save_server1 = strcat(fullfile(folder_mut, 'mut_channel1.mat'));
                save([path_save_server1], 'BF_mut')
                path_save_server2 = strcat(fullfile(folder_mut, 'mut_channel2.mat'));
                save([path_save_server2], 'mCh1_mut')
                path_save_server3 = strcat(fullfile(folder_mut, 'mut_channel3.mat'));
                save([path_save_server3], 'tdT1_mut')
                path_save_server4 = strcat(fullfile(folder_mut, 'mut_channel4.mat'));
                save([path_save_server4], 'tdT2_mut')
                path_save_server5 = strcat(fullfile(folder_mut, 'mut_channel5.mat'));
                save([path_save_server5], 'GFP1_mut')
                path_save_server6 = strcat(fullfile(folder_mut, 'mut_channel6.mat'));
                save([path_save_server6], 'GFP2_mut')
            end
            
    
    % metadata
    reader = bfGetReader(CZI_file_paths_mut{1,1});
    omeMeta = reader.getMetadataStore();
    imageSizeX = omeMeta.getPixelsSizeX(0).getValue();
    imageSizeY = omeMeta.getPixelsSizeY(0).getValue();
    pixelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value();
    pixelSizeX_num = pixelSizeX.doubleValue();
    pixelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value();
    pixelSizeY_num = pixelSizeY.doubleValue();
    pixelSize_um2 = pixelSizeX_num * pixelSizeY_num;
    
    % save metadata
    metadata_mut = [pixelSize_um2 pixelSizeX_num pixelSizeY_num];
    path_save_server_mut = strcat(fullfile(folder_mut, 'metadata_mut.mat'));
    save([path_save_server_mut], 'metadata_mut');
    fprintf('mut metadata .mat saved')


end

 %% BIOFORMATS .CZI IMPORT - MUTANT_b

if analysis_mut_flag
    
    % mutant
    fprintf('\nAnalyzing mutant_b images...\n');
    
    filenames_mut_b = fullfile(folder_mut_b, '*.czi');
    CZI_files_mut_b = dir(filenames_mut_b);
    CZI_file_paths_mut_b = fullfile(folder_mut_b, {CZI_files_mut_b.name});
    

    if split_channels
            
            % RUN TO SAVE ALL IMAGES/CHANNEL AS 1 .MAT:
            
            for ii = 1:length(CZI_files_mut_b);
                reader = bfGetReader(CZI_file_paths_mut_b{1,ii});
                fprintf('Importing mut_b image %d/%d\n', ii, length(CZI_files_mut_b));
                
                BF_mut_b{1,ii} = bfGetPlane(reader,1);
                mCh1_mut_b{1,ii} = bfGetPlane(reader,2);
                tdT1_mut_b{1,ii} = bfGetPlane(reader,3);
                tdT2_mut_b{1,ii} = bfGetPlane(reader,4);
                GFP1_mut_b{1,ii} = bfGetPlane(reader,5);
                GFP2_mut_b{1,ii} = bfGetPlane(reader,6);
            
            end
    end
            
            % save to same folder on server where data is from:
            if save_server
                
                path_save_server1 = strcat(fullfile(folder_mut_b, 'mut_b_channel1.mat'));
                save([path_save_server1], 'BF_mut_b')
                path_save_server2 = strcat(fullfile(folder_mut_b, 'mut_b_channel2.mat'));
                save([path_save_server2], 'mCh1_mut_b')
                path_save_server3 = strcat(fullfile(folder_mut_b, 'mut_b_channel3.mat'));
                save([path_save_server3], 'tdT1_mut_b')
                path_save_server4 = strcat(fullfile(folder_mut_b, 'mut_b_channel4.mat'));
                save([path_save_server4], 'tdT2_mut_b')
                path_save_server5 = strcat(fullfile(folder_mut_b, 'mut_b_channel5.mat'));
                save([path_save_server5], 'GFP1_mut_b')
                path_save_server6 = strcat(fullfile(folder_mut_b, 'mut_b_channel6.mat'));
                save([path_save_server6], 'GFP2_mut_b')
            end
            
    
    % metadata
    reader = bfGetReader(CZI_file_paths_mut_b{1,1});
    omeMeta = reader.getMetadataStore();
    imageSizeX = omeMeta.getPixelsSizeX(0).getValue();
    imageSizeY = omeMeta.getPixelsSizeY(0).getValue();
    pixelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value();
    pixelSizeX_num = pixelSizeX.doubleValue();
    pixelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value();
    pixelSizeY_num = pixelSizeY.doubleValue();
    pixelSize_um2 = pixelSizeX_num * pixelSizeY_num;
    
    % save metadata
    metadata_mut_b = [pixelSize_um2 pixelSizeX_num pixelSizeY_num];
    path_save_server_mut_b = strcat(fullfile(folder_mut_b, 'metadata_mut_b.mat'));
    save([path_save_server_mut_b], 'metadata_mut_b');
    fprintf('mut_b metadata .mat saved')

end

 %% BIOFORMATS .CZI IMPORT - MUTANT_C

if analysis_mut_flag
    
    % mutant
    fprintf('\nAnalyzing mutant_c images...\n');
    
    filenames_mut_c = fullfile(folder_mut_c, '*.czi');
    CZI_files_mut_c = dir(filenames_mut_c);
    CZI_file_paths_mut_c = fullfile(folder_mut_c, {CZI_files_mut_c.name});

    if split_channels

            
            for ii = 1:length(CZI_files_mut_c);
                reader = bfGetReader(CZI_file_paths_mut_c{1,ii});
                fprintf('Importing mut_c image %d/%d\n', ii, length(CZI_files_mut_c));
                
                BF_mut_c{1,ii} = bfGetPlane(reader,1);
                mCh1_mut_c{1,ii} = bfGetPlane(reader,2);
                tdT1_mut_c{1,ii} = bfGetPlane(reader,3);
                tdT2_mut_c{1,ii} = bfGetPlane(reader,4);
                GFP1_mut_c{1,ii} = bfGetPlane(reader,5);
                GFP2_mut_c{1,ii} = bfGetPlane(reader,6);
            
            end
    end
           
            % save to same folder on server where data is from:
            if save_server
                
                path_save_server1 = strcat(fullfile(folder_mut_c, 'mut_c_channel1.mat'));
                save([path_save_server1], 'BF_mut_c')
                path_save_server2 = strcat(fullfile(folder_mut_c, 'mut_c_channel2.mat'));
                save([path_save_server2], 'mCh1_mut_c')
                path_save_server3 = strcat(fullfile(folder_mut_c, 'mut_c_channel3.mat'));
                save([path_save_server3], 'tdT1_mut_c')
                path_save_server4 = strcat(fullfile(folder_mut_c, 'mut_c_channel4.mat'));
                save([path_save_server4], 'tdT2_mut_c')
                path_save_server5 = strcat(fullfile(folder_mut_c, 'mut_c_channel5.mat'));
                save([path_save_server5], 'GFP1_mut_c')
                path_save_server6 = strcat(fullfile(folder_mut_c, 'mut_c_channel6.mat'));
                save([path_save_server6], 'GFP2_mut_c')
            end

    
    % metadata
    reader = bfGetReader(CZI_file_paths_mut_c{1,1});
    omeMeta = reader.getMetadataStore();
    imageSizeX = omeMeta.getPixelsSizeX(0).getValue();
    imageSizeY = omeMeta.getPixelsSizeY(0).getValue();
    pixelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value();
    pixelSizeX_num = pixelSizeX.doubleValue();
    pixelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value();
    pixelSizeY_num = pixelSizeY.doubleValue();
    pixelSize_um2 = pixelSizeX_num * pixelSizeY_num;
    
    % save metadata
    metadata_mut_c = [pixelSize_um2 pixelSizeX_num pixelSizeY_num];
    path_save_server_mut_c = strcat(fullfile(folder_mut_c, 'metadata_mut_c.mat'));
    save([path_save_server_mut_c], 'metadata_mut_c');
    fprintf('mut_c metadata .mat saved')

end


%% BACKGROUND CORRECTION VIA GFP: Group 1 (wt)
%goal is to threshold GFP and use this to background subtract from the tdTomato channel

if analyze_flag

GFPthresh = 400; %change to threshold to fish gut

for m = 1:length(GFP2_wt)
GFP2_thresh_wt = GFP2_wt{1,m};
GFP2_thresh_wt_bg{1,m} = GFP2_thresh_wt > GFPthresh;
figure; imshow(GFP2_thresh_wt_bg{1,m})
end

tdT2thresh = 350; %change threshold to capture adipocytes

for n = 1:length(tdT2_wt)
tdT2_thresh_wt = tdT2_wt{1,n};
tdT2_thresh_wt_bg{1,n} = tdT2_thresh_wt > tdT2thresh;
figure; imshow(tdT2_thresh_wt_bg{1,n})
tdT2_thresh_wt_corr{1,n} = tdT2_thresh_wt_bg{1,n} - GFP2_thresh_wt_bg{1,n};
figure; imshow(tdT2_thresh_mut_corr{1,n})
tdT2_thresh_wt_filt{1,n} = bwareafilt(im2bw(tdT2_thresh_wt_corr{1,n}),[20 500000]);
figure; imshow(imfuse(BF_wt{1,n},tdT2_thresh_wt_filt{1,n},'montage'))
end

end

%% crop images: Group 1 (wt)
if analyze_flag

rectpos_xy_wt = NaN(36,4);
cropwindow = [0,0,1200,700];
vehfile = fullfile(folder_wt,'croppositions_veh.csv');

if crop_wt_flag
%crop tdTomato image
for cc = 1:length(tdT2_wt)
I = tdT2_thresh_wt_filt{1,cc};
figure; imshow(I)
rect = drawrectangle('Position',cropwindow,'Color','r');
pause; %make sure red box is around area of interest and then press return key
currkey=get(gcf,'CurrentKey'); 
        if currkey=='return' 
            currkey==1
        else
            currkey==0
        end
rectpos = get(rect,'Position');
rectpos_xy_wt(cc,1:4) = rectpos;
crop = imcrop(I,rectpos_xy_wt(cc,1:4));
crop_tom{1,cc} = crop;
end
close all
writematrix(rectpos_xy_wt,vehfile);
end

%load_wtcrop_flag = 1;
if load_wtcrop_flag
for cc = 1:length(tdT2_wt)
    I = tdT2_thresh_wt_filt{1,cc};
    rectpos_xy_wt = readmatrix(vehfile);
    crop = imcrop(I,rectpos_xy_wt(cc,1:4));
    crop_tom{1,cc} = crop;
end
end

end

%% Mask, second thresholding and final dilation: Group 1 (wt)

if analyze_flag
    
for p = 1:length(tdT2_wt)
crop_GFP{1,p} = imcrop(GFP2_thresh_wt_bg{1,p},rectpos_xy_wt(p,1:4));
corr_tdT2{1,p} = crop_tom{1,p} - crop_GFP{1,p}; 
figure; imshow(corr_tdT2{1,p})
tdT2_filt{1,p} = bwareafilt(im2bw(corr_tdT2{1,p}),[20 500000]);
figure; imshow(tdT2_filt{1,p})
se = strel('disk',2);
tdT2_clean{1,p} = imopen(tdT2_filt{1,p},se);
figure; imshow(tdT2_clean{1,p})
se_dilate = strel('disk',8);
tdT2_dilate_wt{1,p} = imdilate(tdT2_clean{1,p},se_dilate);
crop_BF_wt{1,p} = imcrop(BF_wt{1,p},rectpos_xy_wt(p,1:4));
figure; imshow(imfuse(crop_BF_wt{1,p},tdT2_dilate_wt{1,p},'montage'))
end


tdT2thresh2 = 430; % change to higher threshold that tdT2thresh to only segment adipocytes

for q = 1:length(tdT2_wt)
crop_tdT2_wt{1,q} = imcrop(tdT2_wt{1,q},rectpos_xy_wt(q,1:4));
tdT2_bg_corr{1,q} = uint16(double(tdT2_dilate_wt{1,q}) .* double(crop_tdT2_wt{1,q}));

%second threshold
tdT2_new = tdT2_bg_corr{1,q};
tdT2_thresh_2{1,q} = tdT2_new > tdT2thresh2;
figure; imshow(tdT2_thresh_2{1,q})
tdT2_clean_3{1,q} = imclearborder(tdT2_thresh_2{1,q});
figure; imshow(imfuse(tdT2_thresh_2{1,q},tdT2_clean_3{1,q},'montage'));
se2 = strel('disk',15);
tdT2_clean_4{1,q} = imdilate(tdT2_clean_3{1,q},se2);
figure;imshow(imfuse(tdT2_clean_3{1,q},tdT2_clean_4{1,q},'montage')); 

tdT2_bg_corr_2{1,q} = uint16(double(tdT2_clean_4{1,q}) .* double(crop_tdT2_wt{1,q}));
figure; image_wt = imshow(imfuse(tdT2_bg_corr_2{1,q},crop_BF_wt{1,q},'montage'));

baseFileName = sprintf('wt_corr%d.tif',q);
fullFileName = fullfile(folder_wt, baseFileName);
saveas(image_wt,fullFileName,'tif')
end
close all

end

%% BACKGROUND CORRECTION VIA GFP_mutant: Group 2 (mut)

if analyze_flag

GFPthresh; %change threshold depending on image


for a = 1:length(GFP2_mut)
GFP2_thresh_mut = GFP2_mut{1,a};
GFP2_thresh_mut_bg{1,a} = GFP2_thresh_mut > GFPthresh;
figure; imshow(GFP2_thresh_mut_bg{1,a})
end

tdT2thresh;

for b = 1:length(tdT2_mut)
tdT2_thresh_mut = tdT2_mut{1,b};
tdT2_thresh_mut_bg{1,b} = tdT2_thresh_mut > tdT2thresh;
figure; imshow(tdT2_thresh_mut_bg{1,b})
tdT2_thresh_mut_corr{1,b} = tdT2_thresh_mut_bg{1,b} - GFP2_thresh_mut_bg{1,b};
figure; imshow(tdT2_thresh_mut_corr{1,b})
tdT2_thresh_mut_filt{1,b} = bwareafilt(im2bw(tdT2_thresh_mut_corr{1,b}),[20 500000]);
figure; imshow(imfuse(BF_mut{1,b},tdT2_thresh_mut_filt{1,b},'montage'))
end

end
%% Crop the images: Group 2 (mut)

if analyze_flag
rectpos_xy_mut = NaN(36,4);
atglfile = fullfile(folder_mut,'croppositions_mut.csv');

%crop tdTomato image
if crop_mut_flag
for dd = 1:length(tdT2_mut)
I2 = tdT2_thresh_mut_filt{1,dd};
figure; imshow(I2)
rect = drawrectangle('Position',cropwindow,'Color','r');
pause; %make sure red box is around area of interest and then press return key
currkey=get(gcf,'CurrentKey'); 
        if currkey=='return'
            
            currkey==1
        else
            currkey==0
        end
rectpos_mut = get(rect,'Position');
rectpos_xy_mut(dd,1:4) = rectpos_mut;
crop_mut = imcrop(I2,rectpos_xy_mut(dd,1:4));
crop_tom_mut{1,dd} = crop_mut;
end
close all
writematrix(rectpos_xy_mut,atglfile);
end

if load_mutcrop_flag
for dd = 1:length(tdT2_mut)% modify
    I2 = tdT2_thresh_mut_filt{1,dd};
    rectpos_xy_mut = readmatrix(atglfile);
    crop_mut = imcrop(I2,rectpos_xy_mut(dd,1:4));
    crop_tom_mut{1,dd} = crop_mut;
end
end
end
%% Mask, second level thresholding and final dilation: Group 2 (mut)
if analyze_flag
    
for e = 1:length(tdT2_mut)
crop_GFP_mut{1,e} = imcrop(GFP2_thresh_mut_bg{1,e},rectpos_xy_mut(e,1:4));
corr_tdT2_mut{1,e} = crop_tom_mut{1,e} - crop_GFP_mut{1,e}; 
figure; imshow(corr_tdT2_mut{1,e})
tdT2_filt_mut{1,e} = bwareafilt(im2bw(corr_tdT2_mut{1,e}),[20 500000]);
figure; imshow(tdT2_filt{1,e})
tdT2_clean_mut{1,e} = imopen(tdT2_filt_mut{1,e},se);
figure; imshow(imfuse(tdT2_filt_mut{1,e},tdT2_clean_mut{1,e},'montage'))
se_dilate = strel('disk',8);
tdT2_dilate_mut{1,e} = imdilate(tdT2_clean_mut{1,e},se_dilate);
crop_BF_mut{1,e} = imcrop(BF_mut{1,e},rectpos_xy_mut(e,1:4));
figure; imshow(imfuse(crop_BF_mut{1,e},tdT2_dilate_mut{1,e},'montage'))
end

tdT2thresh2;

for f = 1:length(tdT2_mut)
crop_tdT2_mut{1,f} = imcrop(tdT2_mut{1,f},rectpos_xy_mut(f,1:4));
tdT2_bg_corr_mut{1,f} = uint16(double(tdT2_dilate_mut{1,f}) .* double(crop_tdT2_mut{1,f}));
figure;imshow(imfuse(tdT2_bg_corr_mut{1,f},crop_BF_mut{1,f},'montage')); %check if the segmentation looks okay!


%second threshold
tdT2_new_mut = tdT2_bg_corr_mut{1,f};
tdT2_thresh_2_mut{1,f} = tdT2_new_mut > tdT2thresh2;
figure; imshow(tdT2_thresh_2_mut{1,f})
tdT2_clean_3_mut{1,f} = imclearborder(tdT2_thresh_2_mut{1,f});
figure; imshow(imfuse(tdT2_clean_2_mut{1,f},tdT2_clean_3_mut{1,f},'montage'));
tdT2_clean_4_mut{1,f} = imdilate(tdT2_clean_3_mut{1,f},se2);
figure;imshow(imfuse(tdT2_clean_3_mut{1,f},tdT2_clean_4_mut{1,f},'montage')); 

tdT2_bg_corr_2_mut{1,f} = uint16(double(tdT2_clean_4_mut{1,f}) .* double(crop_tdT2_mut{1,f}));
figure;image_mut = imshow(imfuse(tdT2_bg_corr_2_mut{1,f},crop_BF_mut{1,f},'montage'));

baseFileName_mut = sprintf('mut_corr%d.tif',f);
fullFileName_mut = fullfile(folder_mut, baseFileName_mut);
saveas(image_mut,fullFileName_mut,'tif')
end
close all

end
%% BACKGROUND CORRECTION VIA GFP_mutant_b: Group 3 (mut_b)

if analyze_flag

GFPthresh; %change threshold depending on image

for a = 1:length(GFP2_mut_b)
GFP2_thresh_mut_b = GFP2_mut_b{1,a};
GFP2_thresh_mut_b_bg{1,a} = GFP2_thresh_mut_b > GFPthresh;
figure; imshow(GFP2_thresh_mut_b_bg{1,a})
end

tdT2thresh;

for b = 1:length(tdT2_mut_b)
tdT2_thresh_mut_b = tdT2_mut_b{1,b};
tdT2_thresh_mut_b_bg{1,b} = tdT2_thresh_mut_b > tdT2thresh;
figure; imshow(tdT2_thresh_mut_b_bg{1,b})
tdT2_thresh_mut_b_corr{1,b} = tdT2_thresh_mut_b_bg{1,b} - GFP2_thresh_mut_b_bg{1,b};
figure; imshow(tdT2_thresh_mut_b_corr{1,b})
tdT2_thresh_mut_b_filt{1,b} = bwareafilt(im2bw(tdT2_thresh_mut_b_corr{1,b}),[20 500000]);
figure; imshow(imfuse(BF_mut_b{1,b},tdT2_thresh_mut_b_filt{1,b},'montage'))
end

end

%% Crop the images: Group 3 (mut_b)

if analyze_flag

rectpos_xy_mut_b = NaN(36,4);
aurafile = fullfile(folder_mut_b,'croppositions_mut_b.csv');

%crop tdTomato image
if crop_mut_b_flag
for dd = 1:length(tdT2_mut_b)% modify
I3 = tdT2_thresh_mut_b_filt{1,dd};
figure; imshow(I3)
rect = drawrectangle('Position',cropwindow,'Color','r');
pause; %make sure red box is around area of interest and then press return key
currkey=get(gcf,'CurrentKey'); 
        if currkey=='return'
            currkey==1
        else
            currkey==0
        end
rectpos_mut_b = get(rect,'Position');
rectpos_xy_mut_b(dd,1:4) = rectpos_mut_b;
crop_mut_b = imcrop(I3,rectpos_xy_mut_b(dd,1:4));
crop_tom_mut_b{1,dd} = crop_mut_b;
end
close all
writematrix(rectpos_xy_mut_b,aurafile);
end

if load_mutcrop_b_flag
for dd = 1:length(tdT2_mut_b)% modify
    I3 = tdT2_thresh_mut_b_filt{1,dd};
    rectpos_xy_mut_b = readmatrix(aurafile);
    crop_mut_b = imcrop(I3,rectpos_xy_mut_b(dd,1:4));
    crop_tom_mut_b{1,dd} = crop_mut_b;
end
end
end

%% Mask, second level thresholding and final dilation: Group 3 (mut_b)

if analyze_flag

for e = 1:length(tdT2_mut_b)
crop_GFP_mut_b{1,e} = imcrop(GFP2_thresh_mut_b_bg{1,e},rectpos_xy_mut_b(e,1:4));
corr_tdT2_mut_b{1,e} = crop_tom_mut_b{1,e} - crop_GFP_mut_b{1,e}; 
figure; imshow(corr_tdT2_mut_b{1,e})
tdT2_filt_mut_b{1,e} = bwareafilt(im2bw(corr_tdT2_mut_b{1,e}),[20 500000]);
figure; imshow(tdT2_filt{1,e})
tdT2_clean_mut_b{1,e} = imopen(tdT2_filt_mut_b{1,e},se);
figure; imshow(imfuse(tdT2_filt_mut_b{1,e},tdT2_clean_mut_b{1,e},'montage'))
se_dilate = strel('disk',8);
tdT2_dilate_mut_b{1,e} = imdilate(tdT2_clean_mut_b{1,e},se_dilate);
figure; imshow(imfuse(tdT2_clean_mut_b{1,e},tdT2_dilate_mut_b{1,e},'montage'))
end


for f = 1:length(tdT2_mut_b)
crop_tdT2_mut_b{1,f} = imcrop(tdT2_mut_b{1,f},rectpos_xy_mut_b(f,1:4));
crop_BF_mut_b{1,f} = imcrop(BF_mut_b{1,f},rectpos_xy_mut_b(f,1:4));
tdT2_bg_corr_mut_b{1,f} = uint16(double(tdT2_dilate_mut_b{1,f}) .* double(crop_tdT2_mut_b{1,f}));
figure;imshow(imfuse(tdT2_bg_corr_mut_b{1,f},crop_BF_mut_b{1,f},'montage')); %check if the segmentation looks okay!


%second threshold
tdT2_new_mut_b = tdT2_bg_corr_mut_b{1,f};
tdT2_thresh_2_mut_b{1,f} = tdT2_new_mut_b > tdT2thresh2;
figure; imshow(tdT2_thresh_2_mut_b{1,f})
tdT2_clean_3_mut_b{1,f} = imclearborder(tdT2_thresh_2_mut_b{1,f});
figure; imshow(imfuse(tdT2_clean_2_mut_b{1,f},tdT2_clean_3_mut_b{1,f},'montage'));
tdT2_clean_4_mut_b{1,f} = imdilate(tdT2_clean_3_mut_b{1,f},se2);
figure;imshow(imfuse(tdT2_clean_3_mut_b{1,f},tdT2_clean_4_mut_b{1,f},'montage')); 

tdT2_bg_corr_2_mut_b{1,f} = uint16(double(tdT2_clean_4_mut_b{1,f}) .* double(crop_tdT2_mut_b{1,f}));
figure;image_mut_b = imshow(imfuse(tdT2_bg_corr_2_mut_b{1,f},crop_BF_mut_b{1,f},'montage'));


baseFileName_mut_b = sprintf('mut_b_corr%d.tif',f);
fullFileName_mut_b = fullfile(folder_mut_b, baseFileName_mut_b);
saveas(image_mut_b,fullFileName_mut_b,'tif')
end
close all

end

%% BACKGROUND CORRECTION VIA GFP_mutant_c: Group 4 (mut_c)

if analyze_flag

GFPthresh; %change threshold depending on image


for a = 1:length(GFP2_mut_c)
GFP2_thresh_mut_c = GFP2_mut_c{1,a};
GFP2_thresh_mut_c_bg{1,a} = GFP2_thresh_mut_c > GFPthresh;
figure; imshow(GFP2_thresh_mut_c_bg{1,a})
end

tdT2thresh;

for b = 1:length(tdT2_mut_c)
tdT2_thresh_mut_c = tdT2_mut_c{1,b};
tdT2_thresh_mut_c_bg{1,b} = tdT2_thresh_mut_c > tdT2thresh;
figure; imshow(tdT2_thresh_mut_c_bg{1,b})
tdT2_thresh_mut_c_corr{1,b} = tdT2_thresh_mut_c_bg{1,b} - GFP2_thresh_mut_c_bg{1,b};
figure; imshow(tdT2_thresh_mut_c_corr{1,b})
tdT2_thresh_mut_c_filt{1,b} = bwareafilt(im2bw(tdT2_thresh_mut_c_corr{1,b}),[20 500000]);
figure; imshow(imfuse(BF_mut_c{1,b},tdT2_thresh_mut_c_filt{1,b},'montage'))
end

end

%% Crop the images: Group 4 (mut_c)
if analyze_flag

rectpos_xy_mut_c = NaN(36,4);
jskfile = fullfile(folder_mut_c,'croppositions_mut_c.csv');

%crop tdTomato image
if crop_mut_c_flag
for dd = 1:length(tdT2_mut_c)% modify
I3 = tdT2_thresh_mut_c_filt{1,dd};
figure; imshow(I3)
rect = drawrectangle('Position',cropwindow,'Color','r');
pause; %make sure circle is around area of interest and then press return key
currkey=get(gcf,'CurrentKey'); 
        if currkey=='return'
            currkey==1
        else
            currkey==0
        end
rectpos_mut_c = get(rect,'Position');
rectpos_xy_mut_c(dd,1:4) = rectpos_mut_c;
crop_mut_c = imcrop(I3,rectpos_xy_mut_c(dd,1:4));
crop_tom_mut_c{1,dd} = crop_mut_c;
end
close all
writematrix(rectpos_xy_mut_c,jskfile);
end

if load_mutcrop_c_flag
for dd = 1:length(tdT2_mut_c)% modify
    I3 = tdT2_thresh_mut_c_filt{1,dd};
    rectpos_xy_mut_c = readmatrix(jskfile);
    crop_mut_c = imcrop(I3,rectpos_xy_mut_c(dd,1:4));
    crop_tom_mut_c{1,dd} = crop_mut_c;
end
end
end

%% Mask, second level thresholding and final dilation: Group 4 (mut_c)

if analyze_flag

for e = 1:length(tdT2_mut_c)
crop_GFP_mut_c{1,e} = imcrop(GFP2_thresh_mut_c_bg{1,e},rectpos_xy_mut_c(e,1:4));
corr_tdT2_mut_c{1,e} = crop_tom_mut_c{1,e} - crop_GFP_mut_c{1,e}; 
figure; imshow(corr_tdT2_mut_c{1,e})
tdT2_filt_mut_c{1,e} = bwareafilt(im2bw(corr_tdT2_mut_c{1,e}),[20 500000]);
figure; imshow(tdT2_filt{1,e})
tdT2_clean_mut_c{1,e} = imopen(tdT2_filt_mut_c{1,e},se);
figure; imshow(imfuse(tdT2_filt_mut_c{1,e},tdT2_clean_mut_c{1,e},'montage'))
se_dilate = strel('disk',8);
tdT2_dilate_mut_c{1,e} = imdilate(tdT2_clean_mut_c{1,e},se_dilate);
figure; imshow(imfuse(tdT2_clean_mut_c{1,e},tdT2_dilate_mut_c{1,e},'montage'))
end


for f = 1:length(tdT2_mut_c)
crop_tdT2_mut_c{1,f} = imcrop(tdT2_mut_c{1,f},rectpos_xy_mut_c(f,1:4));
crop_BF_mut_c{1,f} = imcrop(BF_mut_c{1,f},rectpos_xy_mut_c(f,1:4));
tdT2_bg_corr_mut_c{1,f} = uint16(double(tdT2_dilate_mut_c{1,f}) .* double(crop_tdT2_mut_c{1,f}));
figure;imshow(imfuse(tdT2_bg_corr_mut_c{1,f},crop_BF_mut_c{1,f},'montage')); %check if the segmentation looks okay!


%second threshold
tdT2_new_mut_c = tdT2_bg_corr_mut_c{1,f};
tdT2_thresh_2_mut_c{1,f} = tdT2_new_mut_c > tdT2thresh2;
figure; imshow(tdT2_thresh_2_mut_c{1,f})
tdT2_clean_3_mut_c{1,f} = imclearborder(tdT2_thresh_2_mut_c{1,f});
figure; imshow(imfuse(tdT2_clean_2_mut_c{1,f},tdT2_clean_3_mut_c{1,f},'montage'));
tdT2_clean_4_mut_c{1,f} = imdilate(tdT2_clean_3_mut_c{1,f},se2);
figure;imshow(imfuse(tdT2_clean_3_mut_c{1,f},tdT2_clean_4_mut_c{1,f},'montage')); 

tdT2_bg_corr_2_mut_c{1,f} = uint16(double(tdT2_clean_4_mut_c{1,f}) .* double(crop_tdT2_mut_c{1,f}));
figure;image_mut_c = imshow(imfuse(tdT2_bg_corr_2_mut_c{1,f},crop_BF_mut_c{1,f},'montage'));


baseFileName_mut_c = sprintf('mut_c_corr%d.tif',f);
fullFileName_mut_c = fullfile(folder_mut_c, baseFileName_mut_c);
saveas(image_mut_c,fullFileName_mut_c,'tif')
end
close all

end
%% Quantification

pxCount_tdT2 = NaN(length(tdT2_mut),4);
for aa = 1:length(tdT2_wt)
    pxCount_tdT2(aa,1) = sum(tdT2_clean_4{1,aa}(:) == 1);
end

for bb = 1:length(tdT2_mut)
    pxCount_tdT2(bb,2) = sum(tdT2_clean_4_mut{1,bb}(:) == 1);
end

for c = 1:length(tdT2_mut_b)
    pxCount_tdT2(c,3) = sum(tdT2_clean_4_mut_b{1,c}(:) == 1);
end

for d = 1:length(tdT2_mut_c)
    pxCount_tdT2(d,4) = sum(tdT2_clean_4_mut_c{1,d}(:) == 1);
end

Filename = fullfile(Analysisoutput, 'Plin2pixelcount.csv');
writematrix(pxCount_tdT2,Filename);

%convert into area
area_tdT2 = NaN(length(tdT2_mut),4);
for ff = 1:length(tdT2_wt)
    area_tdT2(ff,1) = pxCount_tdT2(ff,1) .* metadata_wt(1,1);
end

for gg = 1:length(tdT2_mut)
    area_tdT2(gg,2) = pxCount_tdT2(gg,2) .* metadata_mut(1,1);
end

for hh = 1:length(tdT2_mut_b)
    area_tdT2(hh,3) = pxCount_tdT2(hh,3) .* metadata_mut_b(1,1);
end

for ii = 1:length(tdT2_mut_c)
    area_tdT2(ii,4) = pxCount_tdT2(ii,4) .* metadata_mut_c(1,1);
end

Areafilename = fullfile(Analysisoutput, 'Plin2area.csv');
writematrix(area_tdT2,Areafilename);


fprintf('Analysis Complete')
