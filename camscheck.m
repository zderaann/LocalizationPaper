folder = './Vicon_session_2020_12_02/HoloLensRecording__2020_12_02__12_57_18/';
pvnames = {'00132513873758872990', '00132513873872443301', '00132513874232633638'};
depthnames = {'00132513873753984088', '00132513873873738359', '00132513874181976975'};
pvhololens = readtable([folder, 'pv.csv']);
figure();
hold on;
% matdata = pcMatterport.Location(1:10:end, :);
% matcol = pcMatterport.Color(1:10:end, :);
% pcMinMatterport = pointCloud(matdata, 'Color', matcol);
% pcshow( pcMinMatterport, 'MarkerSize', 50); %kazdy desaty cca
grid on;
xlabel("x");
ylabel("y");
zlabel("z");
axis equal;
title('HoloLens Depth Data without transform');
files=dir([folder, 'long_throw_depth/']);
depthhololens = readtable([folder, 'long_throw_depth.csv']);
vlc_ll = readtable([folder, 'vlc_ll.csv']);
vlc_lf = readtable([folder, 'vlc_lf.csv']);
vlc_rf = readtable([folder, 'vlc_rf.csv']);
vlc_rr = readtable([folder, 'vlc_rr.csv']);





for i = 1:3
%DEPTH
    [d_FrameToOrigin, d_CameraViewTransform] = loadData('long_throw_depth\', depthnames{i}, '.pgm', depthhololens);
      
%PV
    [pv_FrameToOrigin, pv_CameraViewTransform] = loadData('pv\' , pvnames{i}, '.ppm', pvhololens);
%VLC    
    [vlcll_FrameToOrigin, vlcll_CameraViewTransform] = loadData('vlc_ll\' , pvnames{i}, '.pgm', vlc_ll);
    [vlclf_FrameToOrigin, vlclf_CameraViewTransform] = loadData('vlc_lf\' , pvnames{i}, '.pgm', vlc_lf);
    [vlcrf_FrameToOrigin, vlcrf_CameraViewTransform] = loadData('vlc_rf\' , pvnames{i}, '.pgm', vlc_rf);
    [vlcrr_FrameToOrigin, vlcrr_CameraViewTransform] = loadData('vlc_rr\' , pvnames{i}, '.pgm', vlc_rr);

end

function [FrameToOrigin, CameraViewTransform] = loadData(prefix, name, sufix, data)
%rho = 1.0260;
%St = [0.999769729972730 0.021430304952185 -0.001108629744351;-0.021424487117784 0.999757878839966 0.005017474379792;0.001215887327680 -0.004992567182154,0.999986797858321];


    poseID = [prefix , name, sufix];
    row = data(strcmp(data.ImageFileName, poseID), :);
    
    FrameToOrigin = [row.FrameToOrigin_m11 row.FrameToOrigin_m12 row.FrameToOrigin_m13 row.FrameToOrigin_m14;
                        row.FrameToOrigin_m21 row.FrameToOrigin_m22 row.FrameToOrigin_m23 row.FrameToOrigin_m24;
                        row.FrameToOrigin_m31 row.FrameToOrigin_m32 row.FrameToOrigin_m33 row.FrameToOrigin_m34;
                        row.FrameToOrigin_m41 row.FrameToOrigin_m42 row.FrameToOrigin_m43 row.FrameToOrigin_m44];
       
     CameraViewTransform = [row.CameraViewTransform_m11 row.CameraViewTransform_m12 row.CameraViewTransform_m13 row.CameraViewTransform_m14;
                              row.CameraViewTransform_m21 row.CameraViewTransform_m22 row.CameraViewTransform_m23 row.CameraViewTransform_m24;
                              row.CameraViewTransform_m31 row.CameraViewTransform_m32 row.CameraViewTransform_m33 row.CameraViewTransform_m34;
                              row.CameraViewTransform_m41 row.CameraViewTransform_m42 row.CameraViewTransform_m43 row.CameraViewTransform_m44];
      
    cameraPointsInWorld = [[diag([1 1 -1]) * 0.05; [0 0 0]] ones(4,1)];            
    transformedCamPoints = zeros(5,3);
    %CamTran = [[inv(CameraViewTransform(1:3, 1:3));0 0 0] [0; 0; 0; 1]];
    CamTran = [[CameraViewTransform(1:3, 1:3);0 0 0] [0; 0; 0; 1]];
    F2O = [[FrameToOrigin(1:3, 1:3);0 0 0] [0; 0; 0; 1]];
    %translation = FrameToOrigin(4,:) + CameraViewTransform(4,:);
    %translation = CameraViewTransform(4,:);
    translation = FrameToOrigin(4,:);
    D2C = inv(CameraViewTransform)'
    
    O2D = FrameToOrigin';
    for m = 1:4
      %cameraPoints(m, :) =  cameraPoints(m, :)  * (FrameToOrigin * CamTran);
      cameraPointsInWorld(m, :) =  (D2C * cameraPointsInWorld(m, :)')' + translation;
      %transformedCamPoints(m, :) = 1/rho * St *  (cameraPoints(m, 1:3)* Rii)' + T;
    end
      drawCamera(cameraPointsInWorld);
      text(cameraPointsInWorld(1,1), cameraPointsInWorld(1,2), cameraPointsInWorld(1,3), prefix);
      hold on;
      axis equal;
      ax = gca;               
      ax.Clipping = 'off';
end


function drawCamera(p)
plot3([p(4, 1), p(1, 1)], [p(4, 2), p(1, 2)], [p(4, 3), p(1, 3)], 'r');
plot3([p(4, 1), p(2, 1)], [p(4, 2), p(2, 2)], [p(4, 3), p(2, 3)], 'g');
plot3([p(4, 1), p(3, 1)], [p(4, 2), p(3, 2)], [p(4, 3), p(3, 3)], 'b');
end