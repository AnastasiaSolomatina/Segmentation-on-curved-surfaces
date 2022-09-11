%--------------------------------------------------------------------------------
% Implementation of preprocessing step before segmentation of the basal
% cell-matrix sites on laminin structure
%
% Input: a tif file with two channels: for cells and ECM
%           a masked image of the membrane
% Author: Anastasia Solomatina
% the code fits the surface to the membrane channel
clear all;
close all;

%% INPUT from the user:
numChann = 2;
numStacks = 35;
numTimeFr = 30;

%% parameters:
Npix = 2000; %number of pixels that are sampled for the surface fitting
numPL = 5;   %numPL - number of levelsets  with particles for pix->part interpolation
numLS = 7;    %numLS - number of levelsets for particle->pixel interpolation

%%
delta = [0.22 0.22 0.7]; %x y z
deltaP = min(delta);  %distance between the particles in the normal direction
num0st = 5;
Next = 3;       %after interp and max proj, I extend the values to the band for floor(numLS/2)+Next levels

%read the initial tif file, create a folder
list = dir('*.tif');
name = list(1).name;
fname = strcat(name(1:end-4),'/');
if ~exist(fname,'dir')
    mkdir(fname)
end
%read the file with the mask for membrane
if size(list,1)==2
    WShmask_name = list(2).name;
end

info = imfinfo(name);
num_images = numel(info);

stPerTime = numChann*numStacks;

sizeInputIm = [info(1).Height, info(1).Width, numStacks];
sizeX = sizeInputIm(1)+2*num0st;
sizeY = sizeInputIm(2)+2*num0st;
sizeZ = sizeInputIm(3)+2*num0st;
Xparticles(1:sizeX,1,1) = delta(1):delta(1):sizeX*delta(1);
Yparticles(1,1:sizeY,1) = delta(2):delta(2):sizeY*delta(2);
Zparticles(1,1,1:sizeZ) = delta(3):delta(3):sizeZ*delta(3);

for timeSt = 1:numTimeFr

    disp([num2str(timeSt), ' time step out of ', num2str(numTimeFr)])
    
    %split channels
    z1 = 0; z2 = 0;

    %allocate arrays for different channels
    sizeInputIm = [info(1).Height, info(1).Width, numStacks];
    channFA_orig = zeros(sizeInputIm);
    channECM_orig = zeros(sizeInputIm);
    
    for stack = stPerTime*(timeSt-1)+1:stPerTime*timeSt
        if mod(stack,numChann) == 0 %the second channel
            z1 = z1+1;
            channECM_orig(:,:,z1) = imread(name,stack);
        end
        if mod(stack,numChann) == 1 %the first channel
            z2 = z2+1;
            channFA_orig(:,:,z2) = imread(name,stack);
            if size(list,1)==2
                chann_WShmask(:,:,z2) = imread(WShmask_name,stack);
            end
        end
    end
    channECM_orig = medfilt3(channECM_orig,[3 3 3]);
    channFA_orig = medfilt3(channFA_orig,[3 3 3]);
    if size(list,1)==2
        chann_WShmask = double(chann_WShmask/max(max(max(chann_WShmask))));
    end
    
    %compute correlation between two channels
    corrM = corrMap(channECM_orig,channFA_orig);
    if timeSt==1
        maxCorr = max(corrM,[],'all');
    end
    corrM = uint16(corrM/maxCorr*double(intmax('uint16')));
    %save the correlation map
    for it_m = 1:size(corrM,3)
        imwrite(uint16(corrM(:,:,it_m)),strcat(fname,'corr','.tif'),'WriteMode','append');
    end

    channECM_init = channECM_orig;
    if size(list,1)==2
        channECM_orig = channECM_orig.*chann_WShmask;
    end
    
    %introduce zero-flux BC
    channFA = add0fluxBC(channFA_orig,num0st);
    channECM = add0Stacks(channECM_orig,num0st);
    
    %create x and y mesh
    X = Xparticles(:,:,1);
    Y = Yparticles(:,:,1);
    
    intensity_t = reshape(channECM,[size(channECM,1)*size(channECM,2)*size(channECM,3),1]);
    quant = quantile(intensity_t,0.95); 
    clear intensity_t
    
    mem_pix = channECM;
    mem_pix(channECM <= quant) = 0.0;
    
    %% pick pixels for the surface fitting
    lin_ind = find(mem_pix);
    %renormalize uniform r.n. by the total number of non-zero elements
    gen_num = uint32(rand(Npix,1)*(size(lin_ind,1)-2)+1);
    %convert generates linear indices into subscripts
    [gen_sub1,gen_sub2,gen_sub3] = ind2sub(size(mem_pix),lin_ind(gen_num));

    % save the ECM with bright pixels on the positions of the picked pixels for fitting
    max_val = 5*max(max(max(channECM)));
    partPos = zeros(Npix,3);
    ECM_px = channECM;
    for it_ind = 1:size(gen_num,1)
        ECM_px(gen_sub1(it_ind,1), gen_sub2(it_ind,1),gen_sub3(it_ind,1)) = max_val;
        partPos(it_ind,1:3) = [gen_sub1(it_ind,1)*delta(1) gen_sub2(it_ind,1)*delta(2) gen_sub3(it_ind,1)*delta(3)];
    end
    
    clear gen_sub1 gen_sub2 gen_sub3 gen_num
    
    ECM_px_size = ECM_px(1+num0st:end-num0st, 1+num0st:end-num0st, 1+num0st:end-num0st);
    
    for it = 1 : size(ECM_px_size,3)
        imwrite(uint16(ECM_px_size(:,:,it)), strcat(fname,'pxls4fit.tif'),'WriteMode','append');
    end
    
     %% fit a surface to the data
     [sf,~,output] = fit([partPos(:,2), partPos(:,3)], partPos(:,1),'poly33');
     
     resid = output.residuals;
     resid2 = resid;
     resid2(resid2>6) = NaN;
     resid2(resid2<-6) = NaN;
     num = sum(isnan(resid2));
     resid2 = repmat(resid2,[1 3]);
     partPos2 = partPos(~isnan(resid2));
     partPos2 = reshape(partPos2,[size(partPos,1)-num 3]);
     
     [sf,~,output] = fit([partPos2(:,2), partPos2(:,3)], partPos2(:,1),'poly33');
    
     clear partPos num resid resid2 partPos2
    
    %% set particles on the surface on the grid Y-Z
    meshY(1:sizeY-2,1) = 2:1:sizeY-1;
    meshY = repmat(meshY,1,sizeZ-2);
    meshZ(1,1:sizeZ-2) = 2:1:sizeZ-1;
    meshZ = repmat(meshZ,sizeY-2,1);
    meshX = sf(meshY*delta(2), meshZ*delta(3));
    sz = size(meshX,1)*size(meshX,2);
    
    %length of each normal vector is 1
    [ny,nz,nx] = surfnorm(sf(meshY*delta(2),meshZ*delta(3)));  %length = sqrt(nx(1,1)*nx(1,1) + ny(1,1)*ny(1,1) + nz(1,1)*nz(1,1)) = 1
    
    %set particles on the surface and extended band
    countP = 0;
    countP_ext = 0;
    half_numLS = floor(numLS/2);
    meshXP = zeros(size(meshX,1), size(meshY,2), numLS);
    meshYP = zeros(size(meshX,1), size(meshY,2), numLS);
    meshZP = zeros(size(meshX,1), size(meshY,2), numLS);
    meshXP_ext = zeros(size(meshX,1), size(meshY,2), numLS+2*Next);
    meshYP_ext = zeros(size(meshX,1), size(meshY,2), numLS+2*Next);
    meshZP_ext = zeros(size(meshX,1), size(meshY,2), numLS+2*Next);
    
    for it_L = -half_numLS-Next : half_numLS+Next
        if ((it_L >=-half_numLS) && (it_L <=half_numLS))
            countP = countP + 1;
            meshXP(:,:,countP) = meshX+deltaP*it_L*nx;
            meshYP(:,:,countP) = meshY*delta(2)+deltaP*it_L*ny;
            meshZP(:,:,countP) = meshZ*delta(3) + deltaP*it_L*nz;
        end
            countP_ext = countP_ext + 1;
            meshXP_ext(:,:,countP_ext) = meshX+deltaP*it_L*nx;
            meshYP_ext(:,:,countP_ext) = meshY*delta(2)+deltaP*it_L*ny;
            meshZP_ext(:,:,countP_ext) = meshZ*delta(3) + deltaP*it_L*nz;
    end
    clear countP countP_ext
    
    surfPart(:,1) = reshape(meshXP,[sz*numLS,1]);
    surfPart(:,2) = reshape(meshYP,[sz*numLS,1]);
    surfPart(:,3) = reshape(meshZP,[sz*numLS,1]);
    
     %% intepolate intensities from the pixels to the particles for the chann_1
     im_size = [sizeX sizeY sizeZ].*delta;
     partVal = pixPartInterpolation(channFA,surfPart,delta,im_size);
      
     %take a max projection and assign the max values to all particles
     %extension to numPL+2*Next level sets
     partVal2 = reshape(partVal,[sz,numLS]);
     partMax = max(partVal2,[],2);
     partMaxAll = repmat(partMax,[1 numLS+2*Next]);
     partMaxAll = reshape(partMaxAll,[size(partMaxAll,1)*size(partMaxAll,2),1]);
     partPos(:,1) = reshape(meshXP_ext,[sz*size(meshXP_ext,3) 1]);
     partPos(:,2) = reshape(meshYP_ext,[sz*size(meshXP_ext,3) 1]);
     partPos(:,3) = reshape(meshZP_ext,[sz*size(meshXP_ext,3) 1]);

     %% create the mask for the particle-to-mesh interpolation correction
     vec_dn2(:,1) = reshape(meshXP(:,:,1),[sz 1]);
     vec_dn2(:,2) = reshape(meshYP(:,:,1),[sz 1]);
     vec_dn2(:,3) = reshape(meshZP(:,:,1),[sz 1]);
     vec_up2(:,1) = reshape(meshXP(:,:,end),[sz 1]);
     vec_up2(:,2) = reshape(meshYP(:,:,end),[sz 1]);
     vec_up2(:,3) = reshape(meshZP(:,:,end),[sz 1]);
     
     sf_up2 = fit([vec_dn2(:,2), vec_dn2(:,3)], vec_dn2(:,1),'poly33');
     sf_dn2 = fit([vec_up2(:,2), vec_up2(:,3)], vec_up2(:,1),'poly33');
    
    pix_mask2 = NaN(size(channFA));
    pix_X_coor(:,1,1) = delta(1):delta(1):size(channFA,1)*delta(1);
    pix_X_coor = repmat(pix_X_coor,[1 size(channFA,2) size(channFA,3)]);
    pix_Y_coor(1,:,1) = delta(2):delta(2):size(channFA,2)*delta(2);
    pix_Y_coor = repmat(pix_Y_coor,[size(channFA,1) 1 size(channFA,3)]);
    pix_Z_coor(1,1,:) = delta(3):delta(3):size(channFA,3)*delta(3);
    pix_Z_coor = repmat(pix_Z_coor,[size(channFA,1) size(channFA,2) 1]);
    
    bound_up2 = sf_up2(pix_Y_coor(1,:,:), pix_Z_coor(1,:,:));
    bound_up2 = repmat(bound_up2,[size(channFA,1) 1 1]);
    bound_dn2 = sf_dn2(pix_Y_coor(1,:,:), pix_Z_coor(1,:,:));
    bound_dn2 = repmat(bound_dn2,[size(channFA,1) 1 1]);
    
    pix_mask2(pix_X_coor > bound_dn2) = 0;
    pix_mask2(pix_X_coor < bound_up2) = 0;
    pix_mask2(isnan(pix_mask2)) = 1;
    
    pix_mask_size2 = pix_mask2(1+num0st:end-num0st, 1+num0st:end-num0st, 1+num0st:end-num0st);
    
    for it = 1 : size(pix_mask_size2,3)
        imwrite(uint16(pix_mask_size2(:,:,it)), strcat(fname,'mask_ECM.tif'),'WriteMode','append');
        imwrite(uint16(channECM_init(:,:,it)),strcat(fname,'mask_ECM.tif'),'WriteMode','append');
    end
    
    %% particle to mesh interpolation
    %remove min intensity background from the particles
    tempA = partMaxAll;
    tempA(tempA==0) = NaN;
    tempA = tempA-min(tempA);
    tempA(isnan(tempA))=0;
    partMaxAll = tempA;
    image_FA = partPixInterpolation(partMaxAll, partPos, [delta(1) delta(2) delta(3)/2], [size(mem_pix,1) size(mem_pix,2) size(mem_pix,3)*2]);
    
    %create extended mask for the membrane
    pix_mask_ext = zeros(size(pix_mask2,1), size(pix_mask2,2), size(pix_mask2,3)*2);
    for it = 1:size(pix_mask2,3)
        pix_mask_ext(:,:,2*it-1) = pix_mask2(:,:,it);
        pix_mask_ext(:,:,2*it) = pix_mask2(:,:,it);
    end
    
    image_FA2 = image_FA.*pix_mask_ext;
    image_FA_size = image_FA2(1+num0st:end-num0st, 1+num0st:end-num0st, 1+2*num0st:end-2*num0st);
    pix_mask_ext_size = pix_mask_ext(1+num0st:end-num0st, 1+num0st:end-num0st, 1+2*num0st:end-2*num0st);
    
    image_FA_init_size = image_FA_size(:,:,2:2:end);
        
    %save 
    for it = 1 : size(image_FA_init_size,3)
        imwrite(uint16(image_FA_init_size(:,:,it)), strcat(fname,'clean_FA.tif'),'WriteMode','append');
    end
    
    clear bound_up bound_dn meshX meshXP meshY meshYP meshZ meshZP
    clear nx ny nz partMax partPos partVal partVal2 pix_X_coor pix_Y_coor pix_Z_coor
  
end

