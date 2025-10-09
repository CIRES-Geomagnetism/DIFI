function [theta_gm, phi_gm, varargout] = gg2gm(theta_gg, phi_gg, varargin)
% transformation between geographic and geomagnetic (dipole) co-ordinates
%          and components
%
% theta_gm                = gg2gm(theta_gg, phi_gg)
% [theta_gm, phi_gm]      = gg2gm(theta_gg, phi_gg)
% [theta_out, phi_out]    = gg2gm(theta_in, phi_in, i_trans)
%                         i_trans = +1: gg -> gm
%                         i_trans = -1: gm -> gg
% [theta, phi, B_theta, B_phi]  = gg2gm(theta_gm, phi_gm, i_trans, B_theta, B_phi)
% [theta, phi, B_theta, B_phi]  = gg2gm(theta_gm, phi_gm, i_trans, B_theta, B_phi, [theta_b phi_b])
% [theta, phi, B_theta, B_phi]  = gg2gm(theta_gm, phi_gm, i_trans, B_theta, B_phi, [g_1^0 g_1^1 h_1^1])
% [theta, phi, R]               = gg2gm(theta_gm, phi_gm)
%
% theta(:) co-latitude [deg]
% phi(:)   longitude [deg]
%
% R(:,2,2) is the matrix that rotates a horizontal vector from GEO to MAG:
%          [B_theta_gm; B_phi_gm] = R*[B_theta_geo; B_phi_geo]
%
% March 2003, Nils Olsen, DTU Space

% September 2004: Output of matrix R added
% September 2011: input argument #6 added, allowing for defining the pole location 
%                 [theta_b phi_b] or [g_1^0 g_1^1 h_1^1]

% June 2015: bug at line 82 corrected (now works when nargin=6)
% September 2025: fixed bug at line 85 (works when varargin{2} is empty)
%
% A. Chulliat, 2025-09-26

rad = pi/180;
if nargin < 3
    i_trans = +1;
else
    i_trans = varargin{1};
end
if abs(i_trans) ~=1; error('i_trans should be +1 or -1'); end

% if nargin > 3
%     if nargin == 5 && nargout == 4
%     else
%         error('5 input and 4 output variables required')
%     end
% end

% Initialization: coordinates of geomagnetic North pole
if nargin > 5
    tmp = varargin{4};
    if length(tmp) == 2; % coordinates theta_b and phi_b of dipole N-pole
        theta_b = tmp(1);
        phi_b = tmp(2);
    elseif length(tmp) == 3; % Gauss coefficients g10, g11, h11
        phi_b = atan2(-tmp(3), -tmp(2))/rad;
        theta_b = atan2(sqrt(tmp(2).^2 + tmp(3).^2), -tmp(1))/rad;
    else
        error('wrong input of input argument #6: neither [theta_b phi_b] nor [g_1^0 g_1^1 h_1^1]');
    end
else
    theta_b = 9.92; phi_b = 287.78;             % epoch 2010.0 (IGRF-11)
end

s_p_b = sin(phi_b*rad);   c_p_b = cos(phi_b*rad);
c_t_b = cos(theta_b*rad); s_t_b = sin(theta_b*rad);

A = [[+c_t_b*c_p_b +c_t_b*s_p_b  -s_t_b]
    [    -s_p_b     +c_p_b     0]
    [+s_t_b*c_p_b +s_t_b*s_p_b  +c_t_b]];
if i_trans == -1; A = A'; end; % gm -> gg

c_t = cos(theta_gg*rad); s_t = sin(theta_gg*rad);
c_p = cos(phi_gg*rad);   s_p = sin(phi_gg*rad);

z = c_t;
x = s_t .* c_p;
y = s_t .* s_p;

x_gm = A(1,1)*x + A(1,2)*y  +A(1,3)*z;
y_gm = A(2,1)*x + A(2,2)*y  +A(2,3)*z;
z_gm = A(3,1)*x + A(3,2)*y  +A(3,3)*z;

theta_gm = 90 - atan2(z_gm, sqrt(x_gm.^2 + y_gm.^2))/rad;
phi_gm = mod(atan2(y_gm, x_gm)/rad, 360);

if nargin >= 5 && ~isempty(varargin{2})
    B_theta = varargin{2};
    B_phi   = varargin{3};
    BE = B_theta.*c_t;
    Bx = BE.*c_p - B_phi.*s_p;
    By = BE.*s_p + B_phi.*c_p;
    Bz = -B_theta.*s_t;
    
    Bx_gm = A(1,1)*Bx + A(1,2)*By  +A(1,3)*Bz;
    By_gm = A(2,1)*Bx + A(2,2)*By  +A(2,3)*Bz;
    Bz_gm = A(3,1)*Bx + A(3,2)*By  +A(3,3)*Bz;
    
    c_t = cos(theta_gm*rad);   s_t = sin(theta_gm*rad);
    c_p = cos(phi_gm*rad);     s_p = sin(phi_gm*rad);
    BE           = Bx_gm.*c_p + By_gm.*s_p;
    varargout{1} = BE.*c_t    - Bz_gm.*s_t; % B_theta
    varargout{2} = By_gm.*c_p - Bx_gm.*s_p; % B_phi    
end

if nargout == 3
    if i_trans ~= 1; error('Calculation of R only for i_trans = 1'); end
    % transformation from GEO/spherical -> GEO/cartesian
    c_t = cos(theta_gg*rad); s_t = sin(theta_gg*rad);
    c_p = cos(phi_gg*rad);   s_p = sin(phi_gg*rad); 
    R_1 = zeros(size(c_p, 1), 3, 3);
    R_1(:,1,:) = [s_t.*c_p  c_t.*c_p  -s_p];
    R_1(:,2,:) = [s_t.*s_p  c_t.*s_p  +c_p];
    R_1(:,3,:) = [c_t       -s_t   zeros(size(c_p))];
    % transformation from MAG/cartesian -> MAG/spherical
    c_t = cos(theta_gm*rad); s_t = sin(theta_gm*rad);
    c_p = cos(phi_gm*rad);   s_p = sin(phi_gm*rad); 
    R_2 = zeros(size(c_p, 1), 3, 3);
    R_2(:,1,:) = [s_t.*c_p  s_t.*s_p  +c_t];
    R_2(:,2,:) = [c_t.*c_p  c_t.*s_p  -s_t];
    R_2(:,3,:) = [-s_p      c_p      zeros(size(c_p))];
    % R_3 is transformation matrix from GEO/spherical -> MAG/spherical
    s_p_b = repmat(s_p_b, size(c_p, 1), 1);
    c_p_b = repmat(c_p_b, size(c_p, 1), 1);
    s_t_b = repmat(s_t_b, size(c_p, 1), 1);
    c_t_b = repmat(c_t_b, size(c_p, 1), 1);
    R_mag_geo = zeros(size(c_p, 1), 3, 3);
    R_mag_geo(:,1,:) = [+c_t_b.*c_p_b +c_t_b.*s_p_b  -s_t_b];
    R_mag_geo(:,2,:) = [    -s_p_b     +c_p_b     zeros(size(c_p, 1), 1)];
    R_mag_geo(:,3,:) = [+s_t_b.*c_p_b +s_t_b.*s_p_b  +c_t_b];
    R_tmp = mat_mul_mat(R_mag_geo, R_1);     
    R_3 = mat_mul_mat(R_2, R_tmp);
    varargout{1} = R_3(:, 2:3, 2:3);
end
