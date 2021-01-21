% Hololens camera coordinate systems

folder = './2020_12_02_Vicon_HoloLens/HoloLensRecording__2020_12_02__12_57_18/';
pvnames = {'00132513873758872990', '00132513873872443301', '00132513874232633638'};
depthnames = {'00132513873753984088', '00132513873873738359', '00132513874181976975'};
pvhololens = readtable([folder, 'pv.csv']);
files=dir([folder, 'long_throw_depth/']);
depthhololens = readtable([folder, 'long_throw_depth.csv']);
vlc_ll = readtable([folder, 'vlc_ll.csv']);
vlc_lf = readtable([folder, 'vlc_lf.csv']);
vlc_rf = readtable([folder, 'vlc_rf.csv']);
vlc_rr = readtable([folder, 'vlc_rr.csv']);

figure(); view(3);
xlabel("x");
ylabel("y");
zlabel("z");
title('HoloLens Depth Data without transform');
hold;

for i = 2
%DEPTH
    loadData('long_throw_depth\', depthnames{i}, '.pgm', depthhololens);
%PV
    loadData('pv\' , pvnames{i}, '.ppm', pvhololens);
%VLC    
    loadData('vlc_ll\' , pvnames{i}, '.pgm', vlc_ll);
    loadData('vlc_lf\' , pvnames{i}, '.pgm', vlc_lf);
    loadData('vlc_rf\' , pvnames{i}, '.pgm', vlc_rf);
    loadData('vlc_rr\' , pvnames{i}, '.pgm', vlc_rr);
end

function [O2D,D2C] = loadData(prefix, name, sufix, data)
%
% x_C = D2C * O2D * x_O
%
poseID = [prefix , name, sufix];
row = data(strcmp(data.ImageFileName, poseID), :);

% T_{Device^{Origin}}
FrameToOrigin = [row.FrameToOrigin_m11 row.FrameToOrigin_m12 row.FrameToOrigin_m13 row.FrameToOrigin_m14;
                 row.FrameToOrigin_m21 row.FrameToOrigin_m22 row.FrameToOrigin_m23 row.FrameToOrigin_m24;
                 row.FrameToOrigin_m31 row.FrameToOrigin_m32 row.FrameToOrigin_m33 row.FrameToOrigin_m34;
                 row.FrameToOrigin_m41 row.FrameToOrigin_m42 row.FrameToOrigin_m43 row.FrameToOrigin_m44];

% T_{Camwera}^{Device}
CameraViewTransform = [row.CameraViewTransform_m11 row.CameraViewTransform_m12 row.CameraViewTransform_m13 row.CameraViewTransform_m14;
                       row.CameraViewTransform_m21 row.CameraViewTransform_m22 row.CameraViewTransform_m23 row.CameraViewTransform_m24;
                       row.CameraViewTransform_m31 row.CameraViewTransform_m32 row.CameraViewTransform_m33 row.CameraViewTransform_m34;
                       row.CameraViewTransform_m41 row.CameraViewTransform_m42 row.CameraViewTransform_m43 row.CameraViewTransform_m44];
      
D2C = inv(CameraViewTransform');
O2D = FrameToOrigin';
    
drawCamera(cameraPointsInWorld);
text(cameraPointsInWorld(1,1), cameraPointsInWorld(1,2), cameraPointsInWorld(1,3), prefix);
end


function drawCamera(p)
plot3([p(4, 1), p(1, 1)], [p(4, 2), p(1, 2)], [p(4, 3), p(1, 3)], 'r');
plot3([p(4, 1), p(2, 1)], [p(4, 2), p(2, 2)], [p(4, 3), p(2, 3)], 'g');
plot3([p(4, 1), p(3, 1)], [p(4, 2), p(3, 2)], [p(4, 3), p(3, 3)], 'b');
end