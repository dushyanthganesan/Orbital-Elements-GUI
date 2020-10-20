%% Orbital elements UI control

% This code shows demonstrates the parameters of an object in earth orbit.
% These parameters are the orbital elements. Each of these elements can be
% changed in real time using a slider and the orbit is shown in the GUI
% plot.

%% Constants 
mu          = 398600;                   % gravitational parameter

%% Arbitrary start values
global oe
a           = 40000;                    % semimajor axis
e           = 0.5;                      % eccentricity
Omega       = deg2rad(30);              % Longitude of the ascending node
inc         = deg2rad(28.5);            % Inclination
omega       = deg2rad(10);              % Argument of the periapsis
oe          = [a;e;Omega;inc;omega;0];  % Orbital Elements
nu          = 0:0.03:2*pi;              % Range of nu

%% Figure setup

f           = figure(1);
panel       = uipanel(f,'Position',[0.061,0.275,0.132,0.490],'BackgroundColor',[1 1 1]);
subplot(1,5,2:5)
earthSphere(50)
xlim([-6e4 5e4])
ylim([-5e4 5e4])
zlim([-3e4 4e4])
ax           = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

%% Semi-Major axis Slider

aSlider     = uicontrol(panel,'Style','slider','Position',[8.9,333.6,164.8,20],'Value',a,'min',8e3, 'max',500e3,'String','a','BackgroundColor',[1 1 1]);
aText       = uicontrol(panel,'Style','text','String',['a = ' sprintf('%g',a) ' km'],'Position',[50.8,355.2,80,20],'BackgroundColor',[1 1 1]);
addlistener(aSlider,'Value','PostSet',@(s,e) listenerVal(oe,aSlider,panel));

%% Eccentricity Slider

eSlider     = uicontrol(panel,'Style','slider','Position',[8.9,278,164.8,20],'Value',e,'min',0, 'max',0.8,'String','e','BackgroundColor',[1 1 1]);
eText       = uicontrol(panel,'Style','text','String',['e = ' sprintf('%g',e)],'Position',[50.8,300,80,20],'BackgroundColor',[1 1 1]);
addlistener(eSlider,'Value','PostSet',@(s,e) listenerVal(oe,eSlider,panel));

%% Longitude of the ascending node, Omega slider

OmSlider    = uicontrol(panel,'Style','slider','Position',[8.9,223,164.8,20],'Value',Omega,'min',0, 'max',2*pi,'String','O','BackgroundColor',[1 1 1]);
OmText      = uicontrol(panel,'Style','text','String',['Omega = ' sprintf('%3.1f',rad2deg(Omega)) ' deg'],'Position',[43.8,245,100,20],'BackgroundColor',[1 1 1]);
addlistener(OmSlider,'Value','PostSet',@(s,e) listenerVal(oe,OmSlider,panel));

%% Inclination Slider

incSlider 	= uicontrol(panel,'Style','slider','Position',[8.9,168,164.8,20],'Value',inc,'min',0, 'max',2*pi,'String','i','BackgroundColor',[1 1 1]);
incText     = uicontrol(panel,'Style','text','String',['inc = ' sprintf('%3.1f',rad2deg(inc)) ' deg'],'Position',[43.8,190,100,20],'BackgroundColor',[1 1 1]);
addlistener(incSlider,'Value','PostSet',@(s,e) listenerVal(oe,incSlider,panel));

%% Argument of the periapsis slider

omegaSlider = uicontrol(panel,'Style','slider','Position',[8.9,113,164.8,20],'Value',omega,'min',0, 'max',2*pi,'String','o','BackgroundColor',[1 1 1]);
omegaText   = uicontrol(panel,'Style','text','String',['omega = ' sprintf('%3.1f',rad2deg(omega)) ' deg'],'Position',[43.8,135,100,20],'BackgroundColor',[1 1 1]);
addlistener(omegaSlider,'Value','PostSet',@(s,e) listenerVal(oe,omegaSlider,panel));

%% Pop Up Menu

popup    = uicontrol(panel,'Style','popupmenu','Position',[18,58,150,20]);
popup.String = {'Preset Orbits','Molniya','Lunar','Tundra','Low-Earth','Medium-Earth','Geosynchronous','High-Earth','Polar','Near Equatorial'};
%popup.Callback = @(s,e) popval;

%% Functions

function oe = rv2oe_Ganesan_Dushyanth(rv,vv,mu)

% Inputs:                                                                %-
%    rPCI:  Cartesian planet-centered inertial (PCI) position (3 by 1)   %-
%    vPCI:  Cartesian planet-centered inertial (PCI) velocity (3 by 1)   %-
%    mu:    gravitational parameter of centrally attacting body.         %-
% Outputs:  orbital elements                                             %-
%    oe(1): semi-major axis.                                             %-
%    oe(2): eccentricity.                                                %-
%    oe(3): longitude of the ascending node (rad)                        %-
%    oe(4): inclination (rad)                                            %-
%    oe(5): argument of the periapsis (rad)                              %-
%    oe(6): true anomaly (rad)                                           %-
%% Inertial Reference frame
Ix = [1 0 0].';
Iy = [0 1 0].';
Iz = [0 0 1].';
%% Orbit properties
r       = norm(rv);                     % magnitude: position
hv      = cross(rv,vv);                 % vector: angular momentum 
h       = norm(hv);                     % magnitude: angular momentum  
p       = h^2/mu;                       % semi-latus rectum

%% eccentricity, e
ev      = cross(vv,hv)/mu - rv/r;       % vector: eccentricity
e       = norm(ev);                     % magnitude: eccentricity

%% Semimajor axis, a
a       = p/(1-e^2);                    % semi-major axis

%% Longitude of the ascending node, Omega
nv      = cross(Iz,hv);                 % Ascending node
n       = norm(nv);                     % magnitude of n

Omega   = atan2(-dot(nv,Iy),-dot(nv,Ix)) + pi; % Longitude of the ascending node

%% Inclination, inc
bv      = cross(nv,Iz);                 % Vector in the orbit plane orthogonal to n
b       = norm(bv);                     % magnitude of vector b (~= 0)
uv      = bv/b;                         % Unit vector in the orbit plane orthogonal to n

inc     = atan2(dot(hv,uv),dot(hv,Iz)); % Inclination of the orbit

%% Argument of the periapsis, omega
% Constructing a right-handed orthonormal basis U
u_n     = nv/n;
u_h     = hv/h;
u_hn    = cross(u_h,u_n);

omega   = atan2(-dot(ev,u_hn),-dot(ev,u_n)) + pi; % Argument of the periapsis

%% True anomaly, nu
% Perifocal basis
px      = ev/e;                         % unit vector in the direction of eccentricity
pz      = hv/h;                         % unit vector along angular momentum
py      = cross(pz,px);

nu      = atan2(-dot(rv,py),-dot(rv,px)) + pi; % True anomaly

%% Orbital elements
oe = [a; e; Omega; inc; omega; nu];

%-%--- IMPORTANT: THE OUTPUT oe MUST BE A COLUMN VECTOR OF LENGTH SIX ---%-
end

function [rPCI,vPCI]  = oe2rv_Ganesan_Dushyanth(oe,mu)

% Input:  orbital elements             (6 by 1 column vector)          
%   oe(1): Semi-major axis.                                            
%   oe(2): Eccentricity.                                               
%   oe(3): Longitude of the ascending node (rad)                       
%   oe(4): Inclination (rad)                                           
%   oe(5): Argument of the periapsis (rad)                             
%   oe(6): True anomaly (rad)                                          
%   mu:    Planet gravitational parameter     (scalar)                 

% Outputs:                                                             
%   rPCI:  Planet-Centered Inertial (PCI) Cartesian position           
%          (3 by 1 column vector)                                      
%   vPCI:  Planet-Centered Inertial (PCI) Cartesian inertial velocity  
%          (3 by 1 column vector)                                      

%% Orbital elements
a       = oe(1);                    % Semi-major axis
e       = oe(2);                    % Eccentricity
Omega   = oe(3);                    % Longitude of the ascending node
inc     = oe(4);                    % Inclination
omega   = oe(5);                    % Argument of the periapsis
nu      = oe(6);                    % True anomaly

%% Transformation 1 about "3"-axis via Omega
% Basis N expressed in I
nxI = [cos(Omega) sin(Omega) 0].';  % nx expressed in I
nyI = [-sin(Omega) cos(Omega) 0].'; % ny expressed in I
nzI = [0 0 1].';                    % nz expressed in I

T_N2I = cat(2,nxI,nyI,nzI);         % Transformation from N to I

%% Transformation 2 about "1"-axis via inc
% Basis Q expressed in N
qxN = [1 0 0].';                    % qx expressed in N
qyN = [0 cos(inc) sin(inc)].';      % qy expressed in N
qzN = [0 -sin(inc) cos(inc)].';     % qz expressed in N

T_Q2N = cat(2,qxN,qyN,qzN);         %  Transformation from Q to N

%% Transformation 3 about "3"-axis via omega
% Basis P expressed in Q
pxQ = [cos(omega) sin(omega) 0].';  % px expressed in Q
pyQ = [-sin(omega) cos(omega) 0].'; % py expressed in Q
pzQ = [0 0 1].';                    % pz expressed in Q

T_P2Q = cat(2,pxQ,pyQ,pzQ);         % Transformation from P to Q

%% Transformation from the basis P to the Inertial basis I
T_P2I = T_N2I * T_Q2N * T_P2Q ;     % Transformation from P to I

%% Position in the perifocal basis
p = a*(1-e^2);                      % semi latus rectum
r = p/(1+e*cos(nu));                % position magnitude from orbit eq

r_p = [r*cos(nu) r*sin(nu) 0];      % position in the perifocal basis

%% Velocity in the perifocal basis
v_p = sqrt(mu/p)*[-sin(nu) e+cos(nu) 0]; % velocity in the perifocal basis

%% Position and velocity in the inertial frame
rv = T_P2I * r_p.' ;                % position in the inertial frame
vv = T_P2I * v_p.' ;                % velocity with respect to the inertial frame

%% Output vectors
rPCI = rv;
vPCI = vv;

end


%Callback
function listenerVal(globaloe,ui,panel,hObject,event)
% General purpose callback function for ui control sliders
global oe
[az,el]     = view;
val         = ui.Value;

mu          = 398600;           % gravitational parameter
a           = oe(1);            % semimajor axis
e           = oe(2);              % eccentricity
Omega       = oe(3);      % Longitude of the ascending node
inc         = oe(4);    % Inclination
omega       = oe(5);      % Argument of the periapsis
if length(ui.String) == 1
    if ui.String == 'e'
        e       = val;
        eText   = uicontrol(panel,'Style','text','String',['e = ' sprintf('%g',e)],'Position',[50.8,300,80,20],'BackgroundColor',[1 1 1]);
    elseif ui.String == 'a'
        a       = val;
        aText   = uicontrol(panel,'Style','text','String',['a = ' sprintf('%g',a) ' km'],'Position',[50.8,355.2,80,20],'BackgroundColor',[1 1 1]);
    elseif ui.String == 'O'
        Omega   = val;
        OmText   = uicontrol(panel,'Style','text','String',['Omega = ' sprintf('%3.1f',rad2deg(Omega)) ' deg'],'Position',[43.8,245,100,20],'BackgroundColor',[1 1 1]);
    elseif ui.String == 'i'
        inc     = val;
        incText = uicontrol(panel,'Style','text','String',['inc = ' sprintf('%3.1f',rad2deg(inc)) ' deg'],'Position',[43.8,190,100,20],'BackgroundColor',[1 1 1]);
    elseif ui.String == 'o'
        omega   = val;
        omegaText   = uicontrol(panel,'Style','text','String',['omega = ' sprintf('%3.1f',rad2deg(omega)) ' deg'],'Position',[43.8,135,100,20],'BackgroundColor',[1 1 1]);
        
    end
% elseif length(ui.String) == 10
%     if ui.Value == 2  % Molniya
%        % Slider.Value = __
%     
%     end
    
    
end
oe          = [a;e;Omega;inc;omega;0];
nu          = 0:0.02:2*pi;
plotOrbit(oe,nu,mu);
view(az,el)
ax           = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

hold off
end

% popup Callback
function popval(hObject,event)
orbit = get(hObject,'Value');

if orbit == 2  % Molniya
    oe = [];
    
end


end

% Orbit Plotter
function p = plotOrbit(oe,nu,mu)
p           = zeros(length(nu),6);
for ii = 1:length(nu)
    [rPCI,vPCI] = oe2rv_Ganesan_Dushyanth([oe(1:5);nu(ii)],mu);
    p(ii,:)     = [rPCI;vPCI].';
end
earthSphere(50)
hold on
plot3(p(:,1),p(:,2),p(:,3),'k','LineWidth',2)
end

% EarthSphere
function [xx,yy,zz] = earthSphere(varargin)
%EARTH_SPHERE Generate an earth-sized sphere.
%   [X,Y,Z] = EARTH_SPHERE(N) generates three (N+1)-by-(N+1)
%   matrices so that SURFACE(X,Y,Z) produces a sphere equal to 
%   the radius of the earth in kilometers. The continents will be
%   displayed.
%
%   [X,Y,Z] = EARTH_SPHERE uses N = 50.
%
%   EARTH_SPHERE(N) and just EARTH_SPHERE graph the earth as a 
%   SURFACE and do not return anything.
%
%   EARTH_SPHERE(N,'mile') graphs the earth with miles as the unit rather
%   than kilometers. Other valid inputs are 'ft' 'm' 'nm' 'miles' and 'AU'
%   for feet, meters, nautical miles, miles, and astronomical units
%   respectively.
%
%   EARTH_SPHERE(AX,...) plots into AX instead of GCA.
% 
%  Examples: 
%    earth_sphere('nm') produces an earth-sized sphere in nautical miles
%
%    earth_sphere(10,'AU') produces 10 point mesh of the Earth in
%    astronomical units
%
%    h1 = gca;
%    earth_sphere(h1,'mile')
%    hold on
%    plot3(x,y,z)
%      produces the Earth in miles on axis h1 and plots a trajectory from
%      variables x, y, and z
%   Clay M. Thompson 4-24-1991, CBM 8-21-92.
%   Will Campbell, 3-30-2010
%   Copyright 1984-2010 The MathWorks, Inc. 

% Input Handling
[cax,args,nargs] = axescheck(varargin{:}); % Parse possible Axes input
error(nargchk(0,2,nargs)); % Ensure there are a valid number of inputs
                           % Handle remaining inputs.
                           % Should have 0 or 1 string input, 0 or 1 numeric input
j = 0;
k = 0;
n = 50; % default value
units = 'km'; % default value
for i = 1:nargs
    if ischar(args{i})
        units = args{i};
        j = j+1;
    elseif isnumeric(args{i})
        n = args{i};
        k = k+1;
    end
end
if j > 1 || k > 1
    error('Invalid input types')
end
%% Calculations
% Scale factors
Scale = {'km' 'm'  'mile'            'miles'           'nm'              'au'                 'ft';
         1    1000 0.621371192237334 0.621371192237334 0.539956803455724 6.6845871226706e-009 3280.839895};
% Identify which scale to use
try
    myscale = 6378.1363*Scale{2,strcmpi(Scale(1,:),units)};
catch %#ok<*CTCH>
    error('Invalid units requested. Please use m, km, ft, mile, miles, nm, or AU')
end

% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;
x = myscale*cosphi*cos(theta);
y = myscale*cosphi*sintheta;
z = myscale*sin(phi)*ones(1,n+1);
%% Plotting
if nargout == 0
    cax = newplot(cax);
    % Load and define topographic data
    load('topo.mat','topo','topomap1');
    % Rotate data to be consistent with the Earth-Centered-Earth-Fixed
    % coordinate conventions. X axis goes through the prime meridian.
    % http://en.wikipedia.org/wiki/Geodetic_system#Earth_Centred_Earth_Fixed_.28ECEF_or_ECF.29_coordinates
    %
    % Note that if you plot orbit trajectories in the Earth-Centered-
    % Inertial, the orientation of the contintents will be misleading.
    topo2 = [topo(:,181:360) topo(:,1:180)]; % #ok<NODEF>
    
    % Define surface settings
    props.FaceColor= 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = topo2;
    % Create the sphere with Earth topography and adjust colormap
    surface(x,y,z,props,'parent',cax)
    colormap(topomap1)
    % Replace the calls to surface and colormap with these lines if you do 
    % not want the Earth's topography displayed.
    %     surf(x,y,z,'parent',cax)
    %     shading flat
    %     colormap gray
    
    % Refine figure
    axis equal
    xlabel(['X [' units ']'])
    ylabel(['Y [' units ']'])
    zlabel(['Z [' units ']'])
    view(127.5,30)
else
    xx = x; yy = y; zz = z;
end
end
