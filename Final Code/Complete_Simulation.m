% This DEMO calculates a Voronoi diagram with arbitrary points in arbitrary
% polytope/polyheron in 2D/3D
clear all;close all;clc
%% Size of our Simulation - number of robots,...
n = 60;        % number of robots
T = 50;       % Number of Time Steps - one step per one time step
speed = 0.1;   % Speed of each Robots = The max distance a robot can cover in each time step
m = 4;         % number of boundary point-candidates
d = 2;         % dimension of the space
tol = 1e-07;            % tolerance value used in "inhull.m" (larger value high precision, possible numerical error)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Points can be choosen uniformly or randomly %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please choose one - Uncomment the one you want and comment the other one
%% Generating Random Values
Xv=rand([1 n]);           % Random x-coordinates inside given boundary
Yv=rand([1 n]);           % Random y-coordinates inside given boundary
pos = [Xv; Yv]';          % Random Initial Coordinates of Robots
%% Generating evenly spaced initial Positions
% Xn=8;       % Robots along x Direction
% Ym=8;       % Robots along Y direction
% n=Xn*Ym;    % Total Number of robots
% nx = linspace(0.001,0.999,Xn);
% ny = linspace(0.001,0.999,Ym);
% [y_mesh,x_mesh] = meshgrid(ny,nx);
% pos = [x_mesh(:),y_mesh(:)];
% Xv = pos(:,1)';
% Yv = pos(:,2)';
%% Defining variable to store distance data
movement = zeros(n,T);    % Distance moved by each point in each time step 
%% Boundary Assignmenmt
bnd0 = [0 1 1 0; 0 0 1 1]';  % Regular Square Boundary
K = convhull(bnd0);
bnd_pnts = bnd0(K,:);   % take boundary points from vertices of convex polytope formed with the boundary point-candidates
%% Initial Plot
[vornb,vorvx] = polybnd_voronoi(pos,bnd_pnts);  % This Function generates voronoi vertices with the boundary
for i = 1:size(vorvx,2)
    col(i,:)= rand(1,3);
end
switch d
    case 2
        figure('position',[0 0 600 600],'Color',[1 1 1]);
        for i = 1:size(pos,1)
        plot(vorvx{i}(:,1),vorvx{i}(:,2),'-r')
        hold on;
        end
        plot(bnd_pnts(:,1),bnd_pnts(:,2),'-');
        hold on;
        plot(pos(:,1),pos(:,2),'Marker','o','MarkerFaceColor',[0 .75 .75],'MarkerEdgeColor','k','LineStyle','none');
        axis('equal')
        axis([0 1 0 1]);
        set(gca,'xtick',[0 1]);
        set(gca,'ytick',[0 1]);        
    case 3
        figure('position',[0 0 600 600],'Color',[1 1 1]);
        for i = 1:size(pos,1)
        K = convhulln(vorvx{i});
        trisurf(K,vorvx{i}(:,1),vorvx{i}(:,2),vorvx{i}(:,3),'FaceColor',col(i,:),'FaceAlpha',0.5,'EdgeAlpha',1)
        hold on;
        end
        scatter3(pos(:,1),pos(:,2),pos(:,3),'Marker','o','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');
        axis('equal')
        axis([0 1 0 1 0 1]);
        set(gca,'xtick',[0 1]);
        set(gca,'ytick',[0 1]);
        set(gca,'ztick',[0 1]);
        set(gca,'FontSize',16);
        xlabel('X');ylabel('Y');zlabel('Z');    
end
%% Making a Video file for Simulation
v=VideoWriter('Video.avi');     % Creating Video File
videoLength = 10;               % Length of the Video
v.FrameRate = T/videoLength;    % Images per second
open(v);
%% Assigning Our Density Function/ Sensory Function
syms Px Py
% sensoryFunc = ((Px-0.5)^4+(Py-0.5)^4)^(1/9);
% sensoryFunc = exp(-9*((Px-0.2)^2+(Py-0.2)^2)) + exp(-9*((Px-0.8)^2+(Py-0.8)^2));
sensoryFunc = exp(-9*((Px-0.5)^2+(Py-0.5)^2));
DelsensoryFunc = [diff(sensoryFunc,Px) diff(sensoryFunc,Py)];
%% Defining Variables for Effecient Calculation of integrals
syms xi yi xi1 yi1
EqIv00 = (1/2)*(yi1-yi)*(xi1+xi);
EqIv01 = (1/6)*( (xi-xi1)*(yi1^2 + yi1*yi + yi^2) + (3*(xi1*yi1^2 - xi*yi^2)) );
EqIv10 = (1/6)*(yi1-yi)*(xi1^2 + xi1*xi +xi^2);
EqIv11 = (1/24)*( (yi1-yi)*(2*xi1*xi*(yi1+yi) + (xi1^2)*(3*yi1+yi) + (xi^2)*(yi1+3*yi)) );
EqIv02 = (1/12)*( (xi-xi1)*(yi1^3 + yi*yi1^2 + yi1*yi^2 + yi^3) + 4*(xi1*yi1^3 - xi*yi^3) );
EqIv20 = (1/12)*(yi1-yi)*(xi1^3 + xi*xi1^2 + xi1*xi^2 + xi^3);
%% Shutting off Figure Display
set(0,'DefaultFigureVisible','off')
%% Initiating Time Step
for t=1:1:T
    for i=1:1:n
        Px=Xv(i);  % X vertices of all Voronoi regions
        Py=Yv(i);  % Y vertices of all voronoi regions
        c1 = (eval(sensoryFunc) - (eval(DelsensoryFunc) * [Px Py].'));
        c2c3 = eval(DelsensoryFunc).';
        vert = cell2mat(vorvx(i));
        Iv00 = 0;
        Iv01 = 0;
        Iv10 = 0;
        Iv11 = 0;
        Iv02 = 0;
        Iv20 = 0;
        for j=1:1:length(vert)-1
            xi = vert(j,1);
            yi = vert(j,2);
            xi1 = vert(j+1,1);
            yi1 = vert(j+1,2);
            % Evaluating the needed values
            Iv00 = Iv00 + eval(EqIv00);
            Iv01 = Iv01 + eval(EqIv01);
            Iv10 = Iv10 + eval(EqIv10);
            Iv11 = Iv11 + eval(EqIv11);
            Iv02 = Iv02 + eval(EqIv02);
            Iv20 = Iv20 + eval(EqIv20);
        end
        % Calculating The Centroid Estimation
        Mvi = Iv00*c1 + [Iv10 Iv01]*c2c3;
        Cv(:,i) = (c1/Mvi)*[Iv10; Iv01] + (1/Mvi)*[Iv20 Iv11; Iv11 Iv02]*c2c3; 
    end
    Cv=Cv';
    for i=1:1:n
        dist=sqrt((pos(i,1)-Cv(i,1))^2 + (pos(i,2)-Cv(i,2))^2);
        movement(i,t)=dist;
        if dist <= speed
            pos(i,:) = Cv(i,:);
        else
            pos(i,1) = ((1-(speed/dist))*pos(i,1)+(speed/dist)*Cv(i,1));
            pos(i,2) = ((1-(speed/dist))*pos(i,2)+(speed/dist)*Cv(i,2));
        end
    end
    Cv=zeros(2,n); % Clearing The values of Centroid for the next time step
    %% Plotting the image
    [~,vorvx] = polybnd_voronoi(pos,bnd_pnts);
    Xv=pos(:,1)';  % Assigning New values to X coordinates
    Yv=pos(:,2)';  % Assigning New values to Y coordinates
    for i = 1:size(vorvx,2)
        col(i,:)= rand(1,3);
    end
    switch d
        case 2
            figure('position',[0 0 600 600],'Color',[1 1 1]);
            for i = 1:size(pos,1)
            plot(vorvx{i}(:,1),vorvx{i}(:,2),'-r')
            hold on;
            end
            plot(bnd_pnts(:,1),bnd_pnts(:,2),'-');
            hold on;
            plot(pos(:,1),pos(:,2),'Marker','o','MarkerFaceColor',[0 .75 .75],'MarkerEdgeColor','k','LineStyle','none');
            axis('equal')
            axis([0 1 0 1]);
            set(gca,'xtick',[0 1]);
            set(gca,'ytick',[0 1]);        
        case 3
            figure('position',[0 0 600 600],'Color',[1 1 1]);
            for i = 1:size(pos,1)
            K = convhulln(vorvx{i});
            trisurf(K,vorvx{i}(:,1),vorvx{i}(:,2),vorvx{i}(:,3),'FaceColor',col(i,:),'FaceAlpha',0.5,'EdgeAlpha',1)
            hold on;
            end
            scatter3(pos(:,1),pos(:,2),pos(:,3),'Marker','o','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');
            axis('equal')
            axis([0 1 0 1 0 1]);
            set(gca,'xtick',[0 1]);
            set(gca,'ytick',[0 1]);
            set(gca,'ztick',[0 1]);
            set(gca,'FontSize',16);
            xlabel('X');ylabel('Y');zlabel('Z');
    end
    %% Taking the image into a Frame
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
h.Visible = 'on';
figure(1)
figure(T+1)
time=1:1:T;
figure(T+2)
r=3;         % Number of robots you want the graph for - MAXIMUM of "n"
for i=1:1:3
    plot(time,movement(i,:));
    hold on;
end
figure(T+2)

