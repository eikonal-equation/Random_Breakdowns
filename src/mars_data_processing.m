%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
location = "../data/";
filename = "jezero_100km_elevation_128ppd.csv";

file = sprintf("%s%s", location, filename);

% Extract data as vector
data = csvread(file);
latitude = data(:,2);
longitude = data(:,3);
elevation = data(:,4);
slope = data(:,5);

% Determine length of square array and reshape
lat_change = latitude(1:end-1) - latitude(2:end);
logical = lat_change < 0;
v = 1:length(latitude);
row_widths = v(logical);
row_width = row_widths(1);

latitude = reshape(latitude, row_width, []);
longitude = reshape(longitude, row_width, []);
elevation = reshape(elevation, row_width, []);
slope = reshape(slope, row_width, []);

% Input average radius of Mars (in kilometers)
r_bar = 3396; % https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use equirectangular projection
ave_lat = mean(latitude(1,:));
ave_long = mean(longitude(:,1));
x_mesh = r_bar*(longitude - ave_long)*cos(ave_lat)*(2*pi/360);
y_mesh = r_bar*(latitude - ave_lat)*(2*pi/360);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE ROUGHNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,m] = size(x_mesh);
roughness = zeros(n,m);
ws = 3; % set window size
for i = 1+ws:n-ws
    for j = 1+ws:m-ws
        current_x_data = latitude(i-ws:i+ws, j-ws:j+ws);
        current_y_data = longitude(i-ws:i+ws, j-ws:j+ws);
        current_data = elevation(i-ws:i+ws, j-ws:j+ws);
        
        X = [reshape(current_x_data, 1, []); reshape(current_y_data, 1, [])]';
        Y = reshape(current_data, [], 1);
        A = (X'*X)\(X'*Y);
        
        roughness(i,j) = sqrt(sum((Y - X*A).^2)/(length(Y)-1)); % equivalent to detrended RMS height
    end
end

% Visualize results
figure;
contourf(x_mesh, y_mesh, elevation, 20); axis equal;
title('Contour plot of elevation')
figure;
contourf(x_mesh, y_mesh, slope, 20, 'LineStyle', 'none'); axis equal;
title('Contour plot of slope')
figure;
contourf(x_mesh, y_mesh, roughness, 20, 'LineStyle', 'none'); axis equal;
title('Contour plot of roughness')


%% Find interesting square locations: Example 4 Valley
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% EXAMPLE 4: VALLEY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

% Determine coordinates of region of interest
i_min = 35;
i_max = 80;
j_min = 30;
j_max = 80;

% Compute equirectangular projection for new area
lat_ex4 = latitude(j_min:j_max, i_min:i_max);
long_ex4 = longitude(j_min:j_max, i_min:i_max);
ave_lat = mean(lat_ex4(1,:));
ave_long = mean(long_ex4(:,1));
x_ex4 = r_bar*(long_ex4 - ave_long)*cos(ave_lat)*(2*pi/360);
y_ex4 = r_bar*(lat_ex4 - ave_lat)*(2*pi/360);

% Extract elevation, slope, and roughness for region of interest
ele_ex4 = elevation(j_min:j_max, i_min:i_max);
slope_ex4 = slope(j_min:j_max, i_min:i_max);
rough_ex4 = roughness(j_min:j_max, i_min:i_max);

% Interpolate on new grid
nx = 500;
ny = 500;
y_min = min(y_ex4(1,:));
y_max = max(y_ex4(1,:));
x_min = min(x_ex4(:,1));
x_max = max(x_ex4(:,1));

scale = (x_max-x_min + y_max - y_min)/2;

[x_mesh_q_ex4, y_mesh_q_ex4] = meshgrid(x_min:((x_max-x_min)/nx):x_max, y_min:((y_max-y_min)/ny):y_max);
ele_interp_ex4 = interpn(x_ex4, y_ex4, ele_ex4, x_mesh_q_ex4, y_mesh_q_ex4, 'linear');
slope_interp_ex4 = interpn(x_ex4, y_ex4, slope_ex4, x_mesh_q_ex4, y_mesh_q_ex4, 'linear');
rough_interp_ex4 = interpn(x_ex4, y_ex4, rough_ex4, x_mesh_q_ex4, y_mesh_q_ex4, 'linear');

% Visualize results
figure;
contourf(x_mesh_q_ex4, y_mesh_q_ex4, ele_interp_ex4, 20); axis equal;
title('Contour plot of interpolated elevation')
figure;
contourf(x_mesh_q_ex4, y_mesh_q_ex4, slope_interp_ex4); axis equal;
title('Contour plot of interpolated slope')
figure;
contourf(x_mesh_q_ex4, y_mesh_q_ex4, rough_interp_ex4); axis equal;
title('Contour plot of interpolated roughness')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% GENERATE SPEED AND BREAKDOWN FUNCTIONS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1_ex4 = zeros(nx,ny);
f2_ex4 = zeros(nx,ny);
f1_line = @(s) (-199/20)*(s) + 200;
f2_line = @(s) (-99/10)*(s) + 100;
for i=1:nx+1
    for j=1:ny+1
        slope_current = slope_interp_ex4(i,j);
        if slope_current > 20
            f1_ex4(i,j) = 1;
        else
            f1_ex4(i,j) = f1_line(slope_current);
        end
        if slope_current > 10
            f2_ex4(i,j) = 1;
        else
            f2_ex4(i,j) = f2_line(slope_current);
        end
    end
end

% rescale to kilometers, and to match grid on [0,1]x[0,1]
f1_ex4 = f1_ex4*scale/1000;
f2_ex4 = f2_ex4*scale/1000;

% Compute breadown rate
phi_ex4 = (rough_interp_ex4).^2/5e3;

% Visualize results
figure;
contourf(x_mesh_q_ex4, y_mesh_q_ex4, f1_ex4, 'LineStyle', 'none');
title('Contour plot of mode 1 speed')
figure;
contourf(x_mesh_q_ex4, y_mesh_q_ex4, f2_ex4, 'LineStyle', 'none');
title('Contour plot of mode 2 speed')
figure;
contourf(x_mesh_q_ex4, y_mesh_q_ex4, phi_ex4, 'LineStyle', 'none');
title('Contour plot of partial breakdown rate')

% Write results to file
filename = "Ex4_f_1.csv";
file = sprintf("%s%s", location, filename);
writematrix(f1_ex4, file);

filename = "Ex4_f_2.csv";
file = sprintf("%s%s", location, filename);
writematrix(f2_ex4, file);

filename = "Ex4_phi.csv";
file = sprintf("%s%s", location, filename);
writematrix(phi_ex4, file);

filename = "Ex4_elevation.csv";
file = sprintf("%s%s", location, filename);
writematrix(ele_interp_ex4, file);

