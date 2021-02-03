%% ------------               -------- Depth data
folder = './B-640/';
files=dir([folder, 'long_throw_depth/']);
depthposes = readtable([folder, 'long_throw_depth.csv']);
figure();
hold on;
grid on;
xlabel("x");
ylabel("y");
zlabel("z");
title('HoloLens Depth Data pointclouds');


ind = 0;
cols = getColors(length(files));
for i = 1:length(files)
   if ~startsWith(files(i).name, 'world_') && endsWith(files(i).name, '.ply')
%        if ind == 10
       name = files(i).name;
       parsed = split(name, '.');
       
       
       poseID = ['long_throw_depth\' , parsed{1}, '.pgm'];
       row = depthposes(strcmp(depthposes.ImageFileName, poseID), :);

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
       depthdata = depthdata*diag([-1 -1 -1]);

       for k = 1:pcdepth.Count
           tmp = ((D2O) * (C2D) * [depthdata(k,:) 1]')';
           depthdata(k, :) =  [tmp(1) tmp(2) tmp(3)];
       end
       cmatrix = ones(pcdepth.Count, 1) * cols{i};   
       pcdepthtransformed = pointCloud(depthdata, 'Color', cmatrix);
       pcshow(pcdepthtransformed, 'MarkerSize', 50);
       hold on;

%        break;

%        end
%       ind = ind + 1;
   end
end

title("(D2O) * C2D * X_C")