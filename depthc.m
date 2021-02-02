folder = './2020_12_02_Vicon_HoloLens/HoloLensRecording__2020_12_02__12_57_18/';
depthtimestamps = {'00132513873764110925', '00132513873774364976'};
numOfPoses = 2;
depthposes = readtable([folder, 'long_throw_depth.csv']);
cols = getColors(numOfPoses);


depthdata = cell(numOfPoses);
figure();
hold on;
for i=1:numOfPoses
    depthpc = pcread([folder, 'long_throw_depth/', depthtimestamps{i}, '.ply']);
    depthPoints = depthpc.Location / 1000;
    

    poseID = ['long_throw_depth\' , depthtimestamps{i}, '.pgm'];
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
    
%     q =r2q(C2D);
%     acos(q(1))*2 /pi *180
    
    for k = 1:depthpc.Count
        tmp = ((D2O) * (C2D) * [-depthPoints(k,:) 1]')';
        depthPoints(k,:) = tmp(1:3);
    end
    cmatrix = ones(depthpc.Count, 1) * cols{i};             
    depthdata{i} = pointCloud(depthPoints, 'Color', cmatrix);
    pcshow(depthdata{i}, 'MarkerSize', 50);
end
title("inv(D2O) * (C2D) * X_C")
grid on;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');