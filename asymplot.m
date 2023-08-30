
% Copyright: Giovanni Cerchiari, Yannick Weiser, Tommaso Faorlin,
%            Thomas Lafenthaler
%
% e-mail: giovanni.cerchiari@uibk.ac.at
% date : 08/2023
% 
% This file is a free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% at your option) any later version.
% 
% This file is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the other repository files.
% If not, it can be found at <https://www.gnu.org/licenses/>.
%
% If you wish to cite this work, we prepared a citation file
% "selective_suppression.bib" in bibtex format in the repository.
%
% This file uses the function "polarPcolor" which was not written
% by the authors of this file and that it is licenced
% in a different way. Copy of the licence for this function is included 
% in the folder "ECheynet-polarPcolor-242fc16" 
% The full file can be downloaded from 
% https://www.mathworks.com/matlabcentral/fileexchange/49040-pcolor-in-polar-coordinates

% This file was used to generate plots of the signal emerging from an
% asymmetric spatial distribution of a light scatterer.
% We assume that light scattering is suppressed via a hemiphserical mirror. 
% We assume a coordinate system 
% (x,y,z) = (cos(theta) cos(phi), sin(theta) sin(phi), cos(theta))

function asymplot()
    % this is the main function of the file
    
    close all
    addpath('.\ECheynet-polarPcolor-242fc16')
    indfigure = 0;
    % polarization x-axis, object x-axis
    indfigure = indfigure + 1;
    plotasymmetrysignal(indfigure, pi/2, 0, pi/2, 0);
    % polarization x-axis, object y-axis
    indfigure = indfigure + 1;
    plotasymmetrysignal(indfigure, pi/2, 0, pi/2, pi/2);
    % polarization x-axis, object z-axis
    indfigure = indfigure + 1;
    plotasymmetrysignal(indfigure, pi/2, 0, 0, 0);
    
    % polarization z-axis, object x-axis
    indfigure = indfigure + 1;
    plotasymmetrysignal(indfigure, 0, 0, pi/2, 0);
    % polarization z-axis, object y-axis
    indfigure = indfigure + 1;
    plotasymmetrysignal(indfigure, 0, 0, pi/2, pi/2);
    % polarization z-axis, object z-axis
    indfigure = indfigure + 1;
    plotasymmetrysignal(indfigure, 0, 0, 0, 0);
end


function plotasymmetrysignal(indfigure, theta_p, phi_p, theta_1, phi_1)
    % this function presents the plot of the radiated power for an
    % arbitrary orientation of linearly polarized light (theta_p, phi_p) 
    % and arbitrary orientation of the object's asymmetry (theta_1, phi_1)
    dim = 512;
    % angles on the sphere
    thetadeg = linspace(0.0,90,dim); % degrees
    phideg = linspace(0,360,dim); % degrees
    theta = thetadeg*pi/180; % radiants
    phi = phideg*pi/180; % radiants

    [Theta, Phi] = meshgrid(theta, phi);
    % orientation of the object theta=0, phi=0
    Theta1 = 0*Theta + theta_1;
    Phi1 = 0*Phi + phi_1;
    % orientation of the laser polarization theta=pi/2, phi=0
    Thetap = 0*Theta + theta_p;
    Phip = 0*Phi + phi_p;

    % contribution of the shape
    asy = asymmetry(Theta, Phi, Theta1, Phi1);
    % contribution of the polarization
    pol = polarization(Theta, Phi, Thetap, Phip)
    
    figure(indfigure)
    [~,c]=polarPcolor(thetadeg,phideg,pol.*asy);
    ylabel(c,'normalized intensity')
    title( sprintf('thetap=%d, phip=%d, theta1=%d, phi1=%d',...
        theta_p*(180/pi), phi_p*(180/pi), theta_1*(180/pi), phi_1*(180/pi)) );
end

function asym = asymmetry(Theta, Phi, Theta1, Phi1)
    % This function calculates the contribution to scattering
    % due to the asymmetry
    asym = cos(Theta).^2 .* cos(Theta1).^2 + sin(Theta).^2 .* sin(Theta1).^2 .* sin(Phi-Phi1).^2;
    asym = asym';
end

function pol = polarization(Theta, Phi, Thetap, Phip)
    % This function calculates the contribution to scattering
    % due to the polarization of the light field
    pol = (1-(ntimesm(Theta, Phi, Thetap, Phip)).^2)';
end

function scalprod = ntimesm(Thetan, Phin, Thetam, Phim)
    % This function calculates the scalar product of two versors
    % given their angular orientation
    nx = sin(Thetan) .* cos(Phin) .* sin(Thetam) .* cos(Phim);
    ny = sin(Thetan) .* sin(Phin) .* sin(Thetam) .* sin(Phim);
    nz = cos(Thetan) .* cos(Thetam);
    scalprod = nx+ny+nz;
end

