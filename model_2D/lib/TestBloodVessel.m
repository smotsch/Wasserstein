% function [ output_args ] = untitled( input_args )

% This function takes in a phsyical grid of dimension M x N
% with "blood vessel" points, then computes the dispersion of the 
% initial condition according to some equation (heat equation with variable 
% diffusion and birth rates which depend on the point's proximity to the 
% "blood vessel" points. 

clc; clear all; close all;

% Load Data
rawData = imread('RawData.jpg');
rhoData = double(rawData(:,:,1)) ;

M = 20;
N = 10;
O = 1000;
BC = 0;

% Sample random points since file is large.
DataSample = rhoData(floor(200*rand(M,1)+1),floor(200*rand(N,1)+1));


x = linspace(0,1,M);
y = linspace(0,1,N);
t = linspace(0,1,O);

dx = x(2) - x(1);
dy = y(2) - y(1);
dt = t(2) - t(1);

%g(distance), measures diffusion in relation to "blood vessels"

%g = @(d) 1./(d+1); %Closer to blood vessel implies slower diffusion

g = @(d) d; %Closer to blood vessel implies faster diffusion

% Create solution vector
rho = zeros(M,N,numel(t));

% Initial/Boundary Condition
[X, Y] = meshgrid(x,y);
rho(:,:,1) = reshape(mvnpdf([X(:) Y(:)], [.5 .5], [1 0; 0 1]), M, N);
rho(1,:,:) = BC; rho(:,1,:) = BC; rho(end,:,:) = BC; rho(:,end,:) = BC; %hold the boundaries at zero. 

% Solve using explicit method, Dirichlet BC
    tic
    for k = 1:numel(t) %k for time
        for i = 2:length(x)-1 %i for x
            for j = 2:length(y)-1 % j for y              
                % compute the distance to closest "blood vessel"
                rho(i,j,k+1) = rho(i,j,k) + (dt/dx^2)*(rho(i+1,j,k) -2*rho(i,j,k) + rho(i-1,j,k)) ...
                + (dt/dy^2)*(rho(i,j+1,k) -2*rho(i,j,k) + rho(i,j-1,k)) + g(DataSample(i,j))*dt*rho(i,j,k)*(1-rho(i,j,k));     
            end
        end
    end
    toc

% end

