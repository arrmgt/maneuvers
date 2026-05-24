function Venu = frd2enu(att, Vfrd)
%BODYFRD_TO_ENU Convert aircraft body velocity FRD to earth ENU.
%
% Inputs
%   att   [3 x N] : [roll; pitch; trueheading] in radians
%                    roll       positive right wing down
%                    pitch      positive nose up
%                    trueheading clockwise from true north
%
%   Vfrd  [N x 3] : [forward, right, down]
%
% Output
%   Venu  [N x 3] : [east, north, up]
%
% Notes
%   Earth frame is local ENU.
%   Heading is true heading, referenced to geodetic true north.

    % checks
    if size(att,1) ~= 3
        error('att must be 3 x N = [roll; pitch; trueheading].');
    end
    if size(Vfrd,2) ~= 3
        error('Vfrd must be N x 3 = [forward, right, down].');
    end
    if size(att,2) ~= size(Vfrd,1)
        error('att is 3xN, so Vfrd must be Nx3 with the same N.');
    end

    % attitude
    roll    = att(1,:).';   % Nx1
    pitch   = att(2,:).';
    heading = att(3,:).';

    % body velocity components
    F = Vfrd(:,1);
    R = Vfrd(:,2);
    D = Vfrd(:,3);

    % trig
    cphi = cos(roll);
    sphi = sin(roll);
    cth  = cos(pitch);
    sth  = sin(pitch);
    cpsi = cos(heading);
    spsi = sin(heading);

    % FRD -> NED
    N =  cth .* cpsi .* F ...
       + (sphi .* sth .* cpsi - cphi .* spsi) .* R ...
       + (cphi .* sth .* cpsi + sphi .* spsi) .* D;

    E =  cth .* spsi .* F ...
       + (sphi .* sth .* spsi + cphi .* cpsi) .* R ...
       + (cphi .* sth .* spsi - sphi .* cpsi) .* D;

    Dn = -sth .* F ...
       +  sphi .* cth .* R ...
       +  cphi .* cth .* D;

    % NED -> ENU
    Venu = [E, N, -Dn];
end