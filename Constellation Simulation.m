%% ASEN 3200 Orbits Projects
% 
%
%%
% Main script (part D) that does the following
%
%  i) read in a JSON constellation design file,
% ii) propagates the constellation in time for a full mean solar day 
%iii) compute the number of spacecraft in LOS of each city i at each t
% iv) 3D plots constellation orbits and the Earth w/ cities/coastlines
%
% Requires 3 Completed Functions:
% loadConstellation.m , propogateState.m , testLoS.m 
%
% CJ Kennedy - SI - 109408903
% 12/1/21
% 
% Full run takes about 200+ seconds--mainly from iii)--so commented out for
% publishing pdf
%
%
%
clear all;clc;close all;
%% i) INITIAL/GIVENS 
fname = 'example_constellation.json'; % constellations
% loadConstellation call
[num_launches, num_spacecraft, satellite_list] = loadConstellation(fname);
% Givens
Re = 6371;
J2 = 1.087e-3;
MU = 3.986e5;
t0 = 0;
%% ii) PROPOGATION of ORBITS of SATELLITES
solarday=86400*1;
t = 1:30:solarday; % mean solar day length with intervals in 30 sec
time = length(t);
numSats = 4; % number of sattlies
x = zeros(6,time*numSats); % intitialize orbit vectors (6,size of time)
for j=1:numSats
    for i=1:time % propogate to get full orbits
        x(1:6,i+(j-1)*time) = propagateState(satellite_list(j).oe0,t(i),t0,MU,J2,Re);    
    end
end

%% iii) CITY SATELLITES LoS 
cities = readtable('worldcities.csv');
sizes = 41001; % length of cities
ReV = Re*ones(sizes,1); 
table = zeros(sizes,time); % rows: city index; columns: timestamp; values: num cities in LoS
longCities = zeros(sizes,1);
latCities = zeros(sizes,1);
population = zeros(sizes,1);
index = 1;
for i=1:sizes
    if ismissing(str2double(table2array(cities(i,10))))==0
        longCities(index) = str2double(table2array(cities(index,4))); 
        latCities(index) = str2double(table2array(cities(index,3))); 
        population(index) = str2double(table2array(cities(index,10))); 
        index = index+1;
    end
end

longCities = deg2rad(longCities); % long and lat in radians
latCities = deg2rad(latCities);
[p q r] = sph2cart(longCities,latCities,ReV); % convert spherical to cartesian coords
% for k=1:numSats % progress through satellites
%     hm = (k-1)*time; % variable to progress through dif sattelites
%     for j =1:time % progress through time
%         for i=1:sizes % progress through cities
%             if testLoS([p(i) q(i) r(i)],[x(1,j+hm) x(2,j+hm) x(3,j+hm)],(15*pi/180))==1
%                 table(i,j)=table(i,j)+1; % add to table count if in LoS
%             end
%         end
%     end
% end
%% iv) COASTLINE
coast = load('world_coastline_low.txt');
sizes = length(coast(:,1)); % size of coastline data
ReV = Re*ones(sizes,1); 
long = ones(sizes,1); % preallocate
lat = ones(sizes,1);
for i = 1:sizes % convert coastline to long and lat
    if coast(i,1) < 0
        long(i) = coast(i,1)+360;
    else
        long(i) = coast(i,1);
    end
end
long = deg2rad(long); % long and lat -> convert to rad
lat = deg2rad(coast(:,2));
[l m n] = sph2cart(long,lat,ReV); % convert spherical to cartesian coords
[L M N] = sphere(18); % create Earth sphere
R = Re*.98; % radius slightly smaller than Earth
%% PLOTTTING iv) and iii)
hold on
axis equal
grid on
set(gca,'Color','k') % black background
ax = gca;
ax.GridColor = '#FFFFFF'; % white grid
xlim([-15000 15000]) % lims of pos
ylim([-15000 15000]) 
zlim([-15000 15000])
xlabel('X (km)') % labels
ylabel('Y (km)')
zlabel('Z (km)')
title('3D Rendering') % title
figure(1)
for i=1:numSats
    index1 = ((i-1)*time+1);
    index2 = (time+(i-1)*time);
    plot3(x(1,index1:index2),x(2,index1:index2),x(3,index1:index2)) % plot orbits
end
% plot3(x(1,1:time),x(2,1:time),x(3,1:time),'Color',	'#800000') % plot 3 orbits
% plot3(x(1,time:time*2),x(2,time:time*2),x(3,time:time*2),'Color','#800080')
% plot3(x(1,time*2:time*3),x(2,time*2:time*3),x(3,1:time*3),'Color','#DB0000')
surf(L*R,M*R,N*R,'edgecolor','#88B2FF','facecolor','#002663') % plot Earth
plot3(l,m,n,'Color','#d9efff') % plot coastlines clearly
plot3(p,q,r,'o','MarkerSize',.2,'Color','#FFFF00') % plot cities
%% INCLUDED COMPLETED FUNCTIONS
function inLoS = testLoS(r_site,r_sc,elevation_limit)
    %DESCRIPTION: Determines whether the spacecraft is within line-of-sight
    %(LoS) of the site given an elevation limit
    %
    %INPUT:
    % r_site            The position vector of the site (km, 3x1)
    % r_sc              The position vector of the spacecraft (km, 3x1)
    % elevation_limit   Lower elevation limit (above the horizon) (rad)
    %
    %OUTPUT: 
    % inLoS             A boolean flag (0 or 1); 1 indicates the spacecraft and
    %                   the site have line-of-sight
    
    %1) Compute whether the site and spacecraft have line of sight (hint, I
    %suggest drawing a picture and writing this constraint as an inequality
    %using a dot product)
    % Create right triangle w/ points r_site, r_sc, and Earth surface
    r_t = r_sc-r_site; 
    % Find theta, angle between sc to Earth and sc to site
    theta = acos(dot(r_t, r_sc)/(norm(r_sc)*norm(r_t))); 
    % Find theta2, angle between horizon and s/c
    theta2 = pi/4-theta; % angle restricted by constraint 
    if theta2>=elevation_limit % return 1 or 0
        inLoS = 1;
    else
        inLoS = 0;
    end
end
function x = propagateState(oe0,t,t_0,MU,J2,Re)
    %DESCRIPTION: Computes the propagated position and velocity in km, km/s
    %accounting for approximate J2 perturbations
    %
    %INPUTS:
    % oe0       Orbit elements [a,e,i,Om,om,f] at time t0 (km,s,rad)
    % t         Current time (s)
    % t0        Time at the initial epoch (s)
    % MU        Central body's gravitational constant (km^3/s^2)
    % J2        Central body's J2 parameter (dimensionless)
    % Re        Radius of central body (km)
    %
    %OUTPUTS:
    % x         Position and velocity vectors of the form [r; rdot] (6x1) at
    %             time t
    
    
    %make sure that function has outputs
    x = NaN(6,1);
    
    %1) Compute the mean orbit elements oe(t) at time t due to J2 perturbations
    a = oe0(1);
    e = oe0(2);
    i = oe0(3);
    Om = oe0(4);
    om = oe0(5);
    f1 = oe0(6);
    
    p = a*(1-e^2);
    n = sqrt(MU/a^3);
    Omdot = -(3/2)*n*J2*(Re/p)*cos(i);
    omdot = (3/2)*n*J2*(Re/p)^2*(2-5/2*sin(i)^2);
    Om = Om+Omdot*(t-t_0);
    om = om+omdot*(t-t_0);
    
    M = n*t;
    
    %2) Solve the time-of-flight problem to compute the true anomaly at tiem t
    
    % Uses Newton's method to solve E -e*sin(E) = M
    % Adopted from Alg 3.1
    table = zeros(1,4);
    if M < pi
        E = M + e/2;
    else
        E = M - e/2;
    end
    count = 0;
    table(1,1) = count;
    table(1,2) = 0;
    table(1,3) = E;
    table(1,4) = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    % Tol
    ratio = 1;
    error = 1*10^-9;
    while abs(ratio) > error
        count = count + 1;
        tmp = table(count,3); % E temp
        ratio = (tmp - e*sin(tmp)-M)/(1-e*cos(tmp));
        tmp = tmp - ratio;
        table(count+1,1) = count; % Iteration
        table(count+1,3) = tmp; % E
        table(count+1,2) = abs(table(count+1,3)-table(count,3)); % delta E
        table(count+1,4) = 2*atan(sqrt((1+e)/(1-e))*tan(tmp/2)); % f
    end
    f = table(count+1,4);
    
    %3) Compute r(t), rdot(t) in the perifocal frame
    % Find Perifocal Frame
    r = p/(1+e*cos(f));
    rV = [r*cos(f); r*sin(f); 0];
    vV = sqrt(MU/p)*[-sin(f); (e+cos(f)); 0];
    %4) Compute r(t), rdot(t) in the ECI frame, save into x
    % Convert to ECI Frame
    Q = [cos(om)*cos(Om)-sin(om)*cos(i)*sin(Om) cos(om)*sin(Om)+sin(om)*cos(i)*cos(Om) sin(om)*sin(i);
        -sin(om)*cos(Om)-cos(om)*cos(i)*sin(Om) -sin(om)*sin(Om)+cos(om)*cos(i)*cos(Om) cos(om)*sin(i);
        sin(i)*sin(Om) -sin(i)*cos(Om) cos(i)]';
    rX = Q*rV;
    vX = Q*vV;
    x = [rX' vX']';
end
function [num_launches, num_spacecraft, satellite_list] = loadConstellation(filename)
%DESCRIPTOIN: Ingests constellation description .json file and parses it
%into a list of structs with full initial orbit elements (km, s, rad) and
%satellite name.
%
%INPUTS:
% filename      A string indicating the name of the .json file to be parsed
%
%OUTPUTS:
% nl            Number of total launches
% ns            Total number of spacecraft between all launches
% satlist       Array of structs with 'name' and 'oe0' properties
%Temporary - just so the function runs the first time you use it.
%You'll need to change all of these!
num_launches = 0;
num_spacecraft = 0;
satellite_list.name = '';
satellite_list.oe0 = NaN(6,1);
%1) extract the constellation structure from the json file
fid = fopen(filename); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
con = jsondecode(str);
%2) read all of the launches and payloads to understand how many launches
% and spacecraft are in the constellation; note, this will be useful in
% Part 2!
insert = 1;
num_launches = length(con.launches);
ae = zeros(6,1);
for i = 1:num_launches
    num_spacecraft = num_spacecraft + length(con.launches(i).payload);
end
%3) RECOMMENDED: Pre-allocate the satellite_list struct
satellite_list(num_spacecraft).name = '';
satellite_list(num_spacecraft).oe0 = NaN(6,1);
%4) Populate each entry in the satellite struct list with its name and
%initial orbit elements [a,e,i,Om,om,f] at time t0
for i = 1:num_launches
    index = length(con.launches(i).payload);
    for j = 1:index
        satellite_list(insert).name = con.launches(i).payload(j).name;
        ae = [con.launches(i).orbit.a con.launches(i).orbit.e ...
            con.launches(i).orbit.i con.launches(i).orbit.Om ...
            con.launches(i).orbit.om con.launches(i).payload(j).f];
        satellite_list(insert).oe0 = ae;
        ae = zeros(6,1);
        insert = insert+1;
    end
end
end