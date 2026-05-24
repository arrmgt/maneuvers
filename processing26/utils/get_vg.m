function [vg,va,mom]=get_vg(bet0,alp0,tas,att0,omega,arm,bfactor,afactor,rolloff,pitoff,hedoff);
%
% calculate ground- and air-relative velocities
%    corrected for moment arm effect
%
% Outputs
%   vg      ENU Earth referenced  air velocity (EAST, NORTH, UP)
%   va      FRD Aircraft reverence air velocity (FORWARD, RIGHT, DOWN)

att1=zeros(size(att0));
%  apply offsets
att1(1,:)=att0(1,:)+rolloff;
att1(2,:)=att0(2,:)+pitoff;
att1(3,:)=att0(3,:)+hedoff;
%  apply upwash/sidewash factors
bet1=bfactor(1).*bet0 + bfactor(2);
alp1=afactor(1).*alp0 + afactor(2);
%
va=tas2va(tas,alp1,bet1);
[yy,zz]=size(omega);
if size(arm,1) >1
    ARM = arm;
else
    ARM = arm(ones(yy,1),:);
end
mom=cross(omega,ARM);
vg = frd2enu(att1,va-mom );

return
