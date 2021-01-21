croppingStart = 0;
croppingEnd = 100;
Binit = 0;
BMax = 100;
%folder = './Vicon_session_2020_12_02/HoloLensRecording__2020_12_02__12_57_18/';
folder = './HoloLensRecording__2020_09_14__11_54_02/';
vicom = readtable('./Vicon_session_2020_12_02/hololens_seq04.txt');
pvhololens = readtable([folder, 'pv.csv']);
Xmin = 1.316;
Xmax = 1.445;
Ymin = -1.2105;
Ymax = -1.145;
Zmin = 1.212;
Zmax = 1.275;
%pri zmene sekvence se musi zmenit i folder dole!

% 
% croppingStart = 100;
% croppingEnd = 0;
% Binit = 120;
% BMax = 160;
% vicom = readtable('./Vicom_2020_08_20/2020_08_20_Vicon_HoloLens_Session/seq1.txt');
% pvhololens = readtable('./Vicom_2020_08_20/HoloLensRecording__2020_08_20__08_36_27/pv.csv');


figure();

indexes = vicom.Var4(:) ~= 1;

%% --------------- Odstraneni chyb
for i = 1:size(vicom,1)
    if vicom.Var5(i) > Xmin && vicom.Var5(i) < Xmax && vicom.Var6(i) > Ymin && vicom.Var6(i) < Ymax && vicom.Var7(i) > Zmin && vicom.Var7(i) < Zmax 
        indexes(i) = 0;        
    end
    
    
end
%% --------------------------------

pcVicom = pointCloud([vicom.Var5(indexes), vicom.Var6(indexes), vicom.Var7(indexes)]);
ViconRot = [vicom.Var5(indexes), vicom.Var6(indexes), vicom.Var7(indexes)];
pcHoloLens = pointCloud([pvhololens.Position_X(croppingStart+1:end-croppingEnd), pvhololens.Position_Y(croppingStart+1:end-croppingEnd), pvhololens.Position_Z(croppingStart+1:end-croppingEnd)]);
pcshowpair(pcHoloLens, pcVicom, 'MarkerSize', 50);

axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title('Hololens and Vicom tracking before registration');


ax = 269;
ay = 69;
az = 248;


Rx = [1 0 0 0; 0 cos(ax) -sin(ax) 0; 0 sin(ax) cos(ax) 0; 0 0 0 1];
Ry = [cos(ay) 0 sin(ay) 0; 0 1 0 0; -sin(ay) 0 cos(ay) 0; 0 0 0 1];
Rz = [cos(az) -sin(az) 0 0; sin(az) cos(az) 0 0; 0 0 1 0; 0 0 0 1];

R = Rx * Ry * Rz;

tform_rotate = affine3d(R);
ptrotHolol = pctransform(pcHoloLens,tform_rotate);

[tform,hololensReg] = pcregrigid(ptrotHolol, pcVicom, 'MaxIterations', 1000);

figure();
pcshowpair(hololensReg, pcVicom, 'MarkerSize', 50);
axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title(['Registered Hololens tracking and Vicom tracking, ax = ', num2str(ax), ', ay = ', num2str(ay), ', az = ', num2str(az)]);

%% --------------------------------------------------------------

numHolol = size(pvhololens,1);
numVicom = size(vicom,1);

params = [];
%inds = [];
vals = [];
initvals = [];


% holFPS = pvhololens.Timestamp(end - croppingEnd) - pvhololens.Timestamp(1 + croppingStart);
% holFPS = holFPS / 10^7;
% holFPS = size(pvhololens,1) / holFPS;
% a = 100/holFPS;
%k = 1 : hololensReg.Count;

hol = hololensReg.Location;
minB = 0;

timediff = pvhololens.Timestamp(2 + croppingStart:end - croppingEnd) - pvhololens.Timestamp(1 + croppingStart:end - croppingEnd -1);
timediff = timediff / 10^7;
cs = cumsum([0; timediff]);
cs  = round(cs * 100) + 1;


for B = Binit:BMax 
    i = cs + B; 
    vic = pcVicom.Location(i, :);
    Rvic = ViconRot(i, :);
    
    fun = @(par)closerTransformation([par(1) par(2) par(3)]', par(4), euler2mat([par(5) par(6) par(7)]), [par(8) par(9) par(10)]', vic, Rvic, hol);
    initpar = [0 0 0 1 0 0 0 0 0 0];
    options = optimset( 'MaxFunEvals', 1000 * 16);
    [bestpar, bestval] = fminsearch(fun, initpar, options);
    
    %t, rho, S, d, C, Rv, D
    initvals = [initvals closerTransformation([initpar(1) initpar(2) initpar(3)]', initpar(4), euler2mat([initpar(5) initpar(6) initpar(7)]), [initpar(8) initpar(9) initpar(10)]', vic, Rvic, hol)];
    
    params = [params bestpar];
    vals = [vals bestval];
    
    if B == Binit
      minpar = bestpar;
      minval = bestval;
      minB = B;
    end
    
    if bestval < minval
      minpar = bestpar;
      minval = bestval; 
      minB = B;
    end
    fprintf(['For B = ', num2str(B), ' the optimized value is ',  num2str(bestval),' parameters are:', '\n' ]);
    bestpar
    
    %t, rho, S, d, C, R, D
end

fprintf(['Best B is ', num2str(minB), ' the optimized value is ',  num2str(minval),' parameters are:', '\n' ]);
minpar

%% -------------------------------------------------------------------------------------------
%pcVicom = pointCloud(vic);
holtrans = zeros(hololensReg.Count, 3);
hol2vicon = zeros(hololensReg.Count, 3);
j = uint64(cs + minB); 
vic = pcVicom.Location(j, :);

Rvic = ViconRot(j, :);

for i = 1:hololensReg.Count
    holtrans(i,:) = (euler2mat([minpar(5) minpar(6) minpar(7)])' * ((1/minpar(4)) * hol(i,:)') - [minpar(8) minpar(9) minpar(10)]')'; %+ Rvic * t cary jinou barvou;
    hol2vicon(i, :) = (holtrans(i, :)' + euler2mat(Rvic(i, :)) * [minpar(1) minpar(2) minpar(3)]')';
end

pcHol = pointCloud(holtrans);
%% ---------------------------------------------------------------------------




figure();
%%

pcshowpair(pcHol, pcVicom, 'MarkerSize', 50);
axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title('Registered Hololens tracking and Vicom tracking after tuning');

%%
hold on;
%kresleni car

cmatrix = ones(size(hol2vicon)).*[0 1 1];
pcHol2Vicon = pointCloud(hol2vicon, 'Color', cmatrix);
pcshow(pcHol2Vicon, 'MarkerSize', 50);
hold on;

for i = 1:hololensReg.Count
    plot3([vic(i, 1), holtrans(i,1)], [vic(i, 2), holtrans(i,2)], [vic(i, 3), holtrans(i,3)], 'w');
    plot3([vic(i, 1), hol2vicon(i,1)], [vic(i, 2), hol2vicon(i,2)], [vic(i, 3), hol2vicon(i,3)], 'r');
end


legend('\color{white} Transformed HoloLens (without Rvic * t)','\color{white} Vicon', '\color{white} Transformed Hololens (with Rvic * t)', '\color{white} Error between HoloLens without Rvic * t and Vicon', '\color{white} Error between HoloLens with Rvic * t and Vicon');


%%
figure();
plot(Binit:BMax , vals, 'b', Binit:BMax , initvals, 'r');

holCheck = pcHoloLens.Location;

rho = minpar(4);
St = euler2mat([minpar(5) minpar(6) minpar(7)])';
Rii = R(1:3, 1:3) * tform.Rotation;
T = St * tform.Translation' - [minpar(8) minpar(9) minpar(10)]';

for i = 1:hololensReg.Count
    holCheck(i,:) = 1/rho * St *  (holCheck(i,:) * Rii)' + T;
end

pcHolCheck = pointCloud(holCheck);

figure();
pcshowpair(pcHolCheck, pcVicom, 'MarkerSize', 50);
axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title('Registered Hololens tracking and Vicom tracking after tuning');


%%
hold on;

plot3(vic(:, 1), vic(:, 2), vic(:, 3), 'ro');

%% ------------------- IN MATTERPORT
pcMatterport = pcread('D:/Documents/Work/B-315/matterport2vicon.ply');
figure();

pcshowpair(pcHolCheck, pcVicom, 'MarkerSize', 50);
hold on;
pcshow(pcMatterport, 'MarkerSize', 50);
axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
title('Transformed HoloLens and Vicon data in transformed Matterport data');
legend('\color{white} Transformed HoloLens','\color{white} Vicon', '\color{white} Matterport');

%% ------------------------ graphs/histograms

errorHolTrans = sum((holtrans - vic).^2, 2);
errorHol2Vicon = sum((hol2vicon - vic).^2, 2);

figure()
hold on;
grid on;
xlabel("pose id");
ylabel("error (m)");

plot(errorHolTrans);
plot(errorHol2Vicon);
title('Graph of errors');
legend('Transformed HoloLens (without Rvic * t)','Transformed HoloLens (with Rvic * t)');


figure()
hold on;
grid on;
xlabel("error (m)");
ylabel("number of occurence");
histogram(errorHolTrans);
histogram(errorHol2Vicon);

title('Histogram of errors');
legend('Transformed HoloLens (without Rvic * t)','Transformed HoloLens (with Rvic * t)');



%% ------------               -------- Depth dat
% CSV DATA:
% Timestamp,ImageFileName,
% Position.X,Position.Y,Position.Z,
% Orientation.W,Orientation.X,Orientation.Y,Orientation.Z,
% FrameToOrigin.m11,FrameToOrigin.m12,FrameToOrigin.m13,FrameToOrigin.m14,
% FrameToOrigin.m21,FrameToOrigin.m22,FrameToOrigin.m23,FrameToOrigin.m24,
% FrameToOrigin.m31,FrameToOrigin.m32,FrameToOrigin.m33,FrameToOrigin.m34,
% FrameToOrigin.m41,FrameToOrigin.m42,FrameToOrigin.m43,FrameToOrigin.m44,
% CameraViewTransform.m11,CameraViewTransform.m12,CameraViewTransform.m13,CameraViewTransform.m14,
% CameraViewTransform.m21,CameraViewTransform.m22,CameraViewTransform.m23,CameraViewTransform.m24,
% CameraViewTransform.m31,CameraViewTransform.m32,CameraViewTransform.m33,CameraViewTransform.m34,
% CameraViewTransform.m41,CameraViewTransform.m42,CameraViewTransform.m43,CameraViewTransform.m44,
% CameraProjectionTransform.m11,CameraProjectionTransform.m12,CameraProjectionTransform.m13,CameraProjectionTransform.m14,
% CameraProjectionTransform.m21,CameraProjectionTransform.m22,CameraProjectionTransform.m23,CameraProjectionTransform.m24,
% CameraProjectionTransform.m31,CameraProjectionTransform.m32,CameraProjectionTransform.m33,CameraProjectionTransform.m34,
% CameraProjectionTransform.m41,CameraProjectionTransform.m42,CameraProjectionTransform.m43,CameraProjectionTransform.m44



figure();
hold on;
matdata = pcMatterport.Location(1:10:end, :);
matcol = pcMatterport.Color(1:10:end, :);
pcMinMatterport = pointCloud(matdata, 'Color', matcol);
pcshow( pcMinMatterport, 'MarkerSize', 50); %kazdy desaty cca
grid on;
xlabel("x");
ylabel("y");
zlabel("z");
title('Matterport');
%legend('\color{white} HoloLens Depth Data', '\color{white} Matterport');

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
       name = files(i).name;
       parsed = split(name, '.');
       
       
       poseID = ['long_throw_depth\' , parsed{1}, '.pgm'];
       %poseID = 'long_throw_depth\00132513873764110925.pgm';
       row = depthposes(strcmp(depthposes.ImageFileName, poseID), :);
       cameraT = [row.Position_X row.Position_Y row.Position_Z];
       cameraQ = [row.Orientation_W row.Orientation_X row.Orientation_Y row.Orientation_Z];
       FrameToOrigin = [row.FrameToOrigin_m11 row.FrameToOrigin_m12 row.FrameToOrigin_m13 row.FrameToOrigin_m14;
                        row.FrameToOrigin_m21 row.FrameToOrigin_m22 row.FrameToOrigin_m23 row.FrameToOrigin_m24;
                        row.FrameToOrigin_m31 row.FrameToOrigin_m32 row.FrameToOrigin_m33 row.FrameToOrigin_m34;
                        row.FrameToOrigin_m41 row.FrameToOrigin_m42 row.FrameToOrigin_m43 row.FrameToOrigin_m44];
       
       CameraViewTransform = [row.CameraViewTransform_m11 row.CameraViewTransform_m12 row.CameraViewTransform_m13 row.CameraViewTransform_m14;
                              row.CameraViewTransform_m21 row.CameraViewTransform_m22 row.CameraViewTransform_m23 row.CameraViewTransform_m24;
                              row.CameraViewTransform_m31 row.CameraViewTransform_m32 row.CameraViewTransform_m33 row.CameraViewTransform_m34;
                              row.CameraViewTransform_m41 row.CameraViewTransform_m42 row.CameraViewTransform_m43 row.CameraViewTransform_m44];
      
%        cameraPoints = [0 0 0 1; 
%                        -0.3 -0.3 -0.3 1;
%                        0.3 -0.3 -0.3 1;
%                        0.3 0.3 -0.3 1;
%                        -0.3 0.3 -0.3 1];
       cameraPoints = [[diag([1 1 -1]) * 0.3; [0 0 0]] ones(4,1)];              
       transformedCamPoints = zeros(5,3);
       Rit = [[Rii'; 0 0 0] [0; 0; 0; 1]];
       CamTran = [[inv(CameraViewTransform(1:3, 1:3));0 0 0] [0; 0; 0; 1]];
       F2O = [[FrameToOrigin(1:3, 1:3);0 0 0] [0; 0; 0; 1]];
       
       D2C = inv(CameraViewTransform)';
       C2D = inv(D2C);
       translation = [FrameToOrigin(4,1) -FrameToOrigin(4,2) -FrameToOrigin(4,3) FrameToOrigin(4,4)]; %+ CameraViewTransform(4,:);
       for m = 1:4
          %cameraPoints(m, :) =  cameraPoints(m, :)  * (FrameToOrigin * CamTran);
          cameraPoints(m, :) =  (D2C * cameraPoints(m, :)')' + translation;
          %transformedCamPoints(m, :) = 1/rho * St *  (cameraPoints(m, 1:3)* Rii)' + T;
       end
       drawCamera(transformedCamPoints);
       cameraR = q2r(cameraQ);
       %cameraR = cameraR;
       %cameraT = (inv(cameraR) * cameraT')';
       cameraTT = (1/rho * St *  (cameraT * Rii)' + T)';    
       %text(transformedCamPoints(1,1), transformedCamPoints(1,2), transformedCamPoints(1,3), parsed{1});
       pose = rigid3d(cameraR(1:3, 1:3),cameraTT);
       %plotCamera('AbsolutePose',pose,'Opacity',0, 'Size', 0.1);
       %pcdepth = pcread([folder, 'long_throw_depth/00132513873764110925.ply']);
       pcdepth = pcread([folder, 'long_throw_depth/', name]);
       depthdata = pcdepth.Location / 1000;
       %cmatrix = ones(size(depthdata)) .* cols{i};
       %cmatrix = ones(size(depthdata)) .* [0 1 0];
       Rx = [1 0 0 0; 0 cos(pi) -sin(pi) 0; 0 sin(pi) cos(pi) 0; 0 0 0 1];
       Rx3 = [1 0 0; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)];
       for k = 1:pcdepth.Count
           %tmp = (Rx * [depthdata(k, 1) depthdata(k, 2) depthdata(k, 3) 1]')' * (F2O * CamTran) + translation;
           %tmp = (((D2C * [depthdata(k, 1) depthdata(k, 2) depthdata(k, 3) 1]')' + translation)')'; 
           tmp = (D2C * ([depthdata(k, 1) depthdata(k, 2) depthdata(k, 3) 1])')' * inv(FrameToOrigin); 
           depthdata(k, :) =  [tmp(1) tmp(2) tmp(3)];
           %depthdata(k, :) = (cameraR(1:3, 1:3)  * depthdata(k, :)')' + cameraT(1:3);
           %depthdata(k, :) = ((1/rho * St *  (depthdata(k, :) * Rii)' + T))';
       end
       pcdepthtransformed = pointCloud(depthdata);
%        pcdepthtransformed = pointCloud(depthdata);
        pcshow(pcdepthtransformed, 'MarkerSize', 50);
       hold on;
       ind = ind + 1;
%        if ind == 10
%            break
%        end
   end
end


%% -------------------------------------------------------------------------------------------
function res = closerTransformation(t, rho, S, d, C, Rv, D)
%     C-known
%     R(e)-known
%     t-unknown
%     rho-unknown
%     S(e) - unknown
%     D - known
%     d - unknown
    
    num = size(Rv,1);
    res = 0;
    
    for i = 1:num
        R = euler2mat(Rv(i,:));
        f = C(i,:)' + R * t - (1/rho) * S' * D(i,:)' - d;
%         res = res + norm(f)^2;
        res = res + log(norm(f) + 1);
    end
    
    res = res + 10 * log(norm(t) + 1) + 10 * log(norm(d) + 1);
    %sum(norm(C + R*t - 1/ro * S' * D - d))
end

function R = euler2mat(e)
    x = e(1);
    y = e(2);
    z = e(3);
    R = [(cos(y) * cos(z)) (-1 * cos(y) * sin(z)) sin(y);
         (cos(x) * sin(z) + sin(x) * sin(y) * cos(z)) (cos(x) * cos(z) - sin(x) * sin(y) * sin(z)) (-1 * sin(x) * cos(y));
         (sin(x) * sin(z) - cos(x) * sin(y) * cos(z)) (sin(x) * cos(z) + cos(x) * sin(y) * sin(z)) (cos(x) * cos(y))];

end

function drawCamera(p)
% plot3([p(1, 1), p(2, 1)], [p(1, 2), p(2, 2)], [p(1, 3), p(2, 3)], 'r');
% plot3([p(1, 1), p(3, 1)], [p(1, 2), p(3, 2)], [p(1, 3), p(3, 3)], 'r');
% plot3([p(1, 1), p(4, 1)], [p(1, 2), p(4, 2)], [p(1, 3), p(4, 3)], 'r');
% plot3([p(1, 1), p(5, 1)], [p(1, 2), p(5, 2)], [p(1, 3), p(5, 3)], 'r');
% plot3([p(2, 1), p(3, 1)], [p(2, 2), p(3, 2)], [p(2, 3), p(3, 3)], 'r');
% plot3([p(3, 1), p(4, 1)], [p(3, 2), p(4, 2)], [p(3, 3), p(4, 3)], 'r');
% plot3([p(4, 1), p(5, 1)], [p(4, 2), p(5, 2)], [p(4, 3), p(5, 3)], 'r');
% plot3([p(5, 1), p(2, 1)], [p(5, 2), p(2, 2)], [p(5, 3), p(2, 3)], 'r');


plot3([p(4, 1), p(1, 1)], [p(4, 2), p(1, 2)], [p(4, 3), p(1, 3)], 'r', 'LineWidth', 1);
plot3([p(4, 1), p(2, 1)], [p(4, 2), p(2, 2)], [p(4, 3), p(2, 3)], 'g', 'LineWidth', 1);
plot3([p(4, 1), p(3, 1)], [p(4, 2), p(3, 2)], [p(4, 3), p(3, 3)], 'b', 'LineWidth', 1);
end


