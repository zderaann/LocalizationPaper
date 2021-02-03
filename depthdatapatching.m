folder = './Vicom_2020_08_20/HoloLensRecording__2020_08_20__08_36_27/';
pvhololens = readtable([folder, 'pv.csv']);
files=dir([folder, 'long_throw_depth/']);
depthposes = readtable([folder, 'long_throw_depth.csv']);
figure();
hold on;
grid on;
xlabel("x");
ylabel("y");
zlabel("z");
title('HoloLens Depth Data pointclouds');
pvtimestamps = pvhololens.Timestamp;

ind = 0;
cols = getColors(length(files));
for i = 1:length(files)
   if ~startsWith(files(i).name, 'world_') && endsWith(files(i).name, '.ply')
%         if ind == 47
       name = files(i).name;
       parsed = split(name, '.');
       
       
       poseID = ['long_throw_depth\' , parsed{1}, '.pgm'];
       row = depthposes(strcmp(depthposes.ImageFileName, poseID), :);

       timestamp = row.Timestamp;
       FrameToOrigin = [row.FrameToOrigin_m11 row.FrameToOrigin_m12 row.FrameToOrigin_m13 row.FrameToOrigin_m14;
                        row.FrameToOrigin_m21 row.FrameToOrigin_m22 row.FrameToOrigin_m23 row.FrameToOrigin_m24;
                        row.FrameToOrigin_m31 row.FrameToOrigin_m32 row.FrameToOrigin_m33 row.FrameToOrigin_m34;
                        row.FrameToOrigin_m41 row.FrameToOrigin_m42 row.FrameToOrigin_m43 row.FrameToOrigin_m44];
       
       CameraViewTransform = [row.CameraViewTransform_m11 row.CameraViewTransform_m12 row.CameraViewTransform_m13 row.CameraViewTransform_m14;
                              row.CameraViewTransform_m21 row.CameraViewTransform_m22 row.CameraViewTransform_m23 row.CameraViewTransform_m24;
                              row.CameraViewTransform_m31 row.CameraViewTransform_m32 row.CameraViewTransform_m33 row.CameraViewTransform_m34;
                              row.CameraViewTransform_m41 row.CameraViewTransform_m42 row.CameraViewTransform_m43 row.CameraViewTransform_m44];

        C2D = inv(CameraViewTransform)';
        D2C = [C2D(1:3,1:3)' -C2D(1:3,1:3)' * C2D(1:3,4); 0 0 0 1];
        D2O = FrameToOrigin'; 
        O2D = inv(D2O);


       pcdepth = pcread([folder, 'long_throw_depth/', name]);
       depthdata = pcdepth.Location / 1000;
       depthdata = depthdata * diag([-1 -1 -1]);
       
       [M, closestIndex] = min(abs(pvtimestamps - timestamp));
       closestImageName = pvhololens.ImageFileName(closestIndex);
       imgname = split(closestImageName{1}, '.');
       closestImage = imread([folder, imgname{1}, '.jpg']);
       

       pvFrameToOrigin = [pvhololens.FrameToOrigin_m11(closestIndex) pvhololens.FrameToOrigin_m12(closestIndex) pvhololens.FrameToOrigin_m13(closestIndex) pvhololens.FrameToOrigin_m14(closestIndex);
                          pvhololens.FrameToOrigin_m21(closestIndex) pvhololens.FrameToOrigin_m22(closestIndex) pvhololens.FrameToOrigin_m23(closestIndex) pvhololens.FrameToOrigin_m24(closestIndex);
                          pvhololens.FrameToOrigin_m31(closestIndex) pvhololens.FrameToOrigin_m32(closestIndex) pvhololens.FrameToOrigin_m33(closestIndex) pvhololens.FrameToOrigin_m34(closestIndex);
                          pvhololens.FrameToOrigin_m41(closestIndex) pvhololens.FrameToOrigin_m42(closestIndex) pvhololens.FrameToOrigin_m43(closestIndex) pvhololens.FrameToOrigin_m44(closestIndex)];

       pvCameraViewTransform = [pvhololens.CameraViewTransform_m11(closestIndex) pvhololens.CameraViewTransform_m12(closestIndex) pvhololens.CameraViewTransform_m13(closestIndex) pvhololens.CameraViewTransform_m14(closestIndex);
                              pvhololens.CameraViewTransform_m21(closestIndex) pvhololens.CameraViewTransform_m22(closestIndex) pvhololens.CameraViewTransform_m23(closestIndex) pvhololens.CameraViewTransform_m24(closestIndex);
                              pvhololens.CameraViewTransform_m31(closestIndex) pvhololens.CameraViewTransform_m32(closestIndex) pvhololens.CameraViewTransform_m33(closestIndex) pvhololens.CameraViewTransform_m34(closestIndex);
                              pvhololens.CameraViewTransform_m41(closestIndex) pvhololens.CameraViewTransform_m42(closestIndex) pvhololens.CameraViewTransform_m43(closestIndex) pvhololens.CameraViewTransform_m44(closestIndex)];
        
       projectedDepthData = ones(size(depthdata, 1), 2);
       pvC2D = inv(pvCameraViewTransform)';
       pvD2O = pvFrameToOrigin'; 
       pvO2D = inv(pvD2O);
       pvD2C = inv(pvC2D);
       cmatrix = ones(pcdepth.Count, 3);
       
       f = [1038.135254 1036.468140];
       pp = [664.387146 396.142090];
       K = [f(1) 0 pp(1); 0 f(2) pp(2); 0 0 1];
       visible = true(size(depthdata, 1), 1);
       
       for k = 1:pcdepth.Count
           tmp = ((D2O) * (C2D) * [depthdata(k,:) 1]')';
           depthdata(k, :) =  [tmp(1) tmp(2) tmp(3)];
           
           tmp2 = (pvD2C * pvO2D * [depthdata(k,:) 1]');
           %tmp3 = (([-tmp2(1) * 1000 tmp2(2) * 1000 tmp2(3)])/tmp2(3) + [1344/2 756/2 0]);
           %projectedDepthData(k, :) = (h2a(tmp3'))';
           projectedDepthData(k, :) = (h2a(K * tmp2(1:3)))';
           
           
           cmatrix(k,:) = closestImage(max(1, min(round(projectedDepthData(k, 2)), 756)), 1345 - max(1, min(round(projectedDepthData(k, 1)), 1344)), :);
           if round(projectedDepthData(k, 1)) < 1 || round(projectedDepthData(k, 1)) > 1344 || round(projectedDepthData(k, 2)) < 1 || round(projectedDepthData(k, 2)) > 756
            visible(k) = false;
           end
       end
       
       %cmatrix = ones(pcdepth.Count, 1) * cols{i}; 
       colors = cmatrix(visible, :)/255;
       pcdepthtransformed = pointCloud(depthdata(visible, :), 'Color', colors);
       %pcdepthprojected = pointCloud(projectedDepthData);
       pcshow(pcdepthtransformed, 'MarkerSize', 50);
       hold on;
       
%        figure()
%        hold on;
       %pcshow(pcdepthprojected, 'MarkerSize', 50);
%        for m = 1:
%        plot()
%        imshow(0.5 * closestImage);
%        hold on;
%        
%       visibleprojected = projectedDepthData(visible, :);
%        for m = 1:size(visibleprojected,1)
%            plot(visibleprojected(m,1), visibleprojected(m,2), '.', 'Color', colors(m, :));
%        end
       
        
%    break;
%         
%    
%        end
   ind = ind + 1;
   end
end
title("(D2O) * C2D * X_C")