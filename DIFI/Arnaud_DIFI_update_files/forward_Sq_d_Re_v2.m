function [B_1,B_2] = forward_Sq_d_Re_v2(r,theta,phi,t,f107,s)

% [B_1,B_2] = forward_Sq_d_Re_v2(r,theta,phi,t,f107,s)
%
% Calculate the primary (1) and induced (2) Sq magnetic field from a model
% in DIPOLE coordinates (parameterization of Sabaka et al., 2002).
%
% "v2" version: 
% (a) fixes a bug in the calculation of the primary field above currents 
% that was present in earlier versions,
% (b) uses gg2gm_v2.m instead of gg2gm_2010.m,
% (c) uses data chunking to run faster.
%
% Inputs:   
%       r(:)                    radius [km]
%       theta(:), phi(:)        geographic co-latitude, longitude [deg]
%       t(:)                    time [MJD2000]
%       f107(:)                 F10.7 index
%       s is a structure with the following elements                       
%           m_e_d_Re(:,:)           primary Sq model (REAL coefficients)
%           m_i_d_Re(:,:)           secondary Sq model (REAL coefficients)
%           nmax, mmax              max degree and order in dipole coord.
%           p_vec(:)                diurnal wavenumbers
%           s_vec(:)                seasonal wavenumbers
%           theta_NGP               co-latitude of N. Geomagnetic Pole [deg]
%           phi_NGP                 longitude of N. Geomagnetic Pole [deg]
%           h                       altitude of Sq currents [km]
%           N                       Woldf ratio of the F10.7 dependence
%
% Dependendies :
% - "gg2gm_v2", "jd2000", "getmut"
% - "design_SHA_Sq_i_Re_v2" and "design_SHA_Sq_e_Re_v2"
%
% Earlier versions: 2011-04-23, 2016-09-22 ("forward_Sq_d_Re.m")
%
% A. Chulliat, 2025-09-28

% maximum memory
function [B_1,B_2] = forward_Sq_d_Re_v2(r,theta,phi,t,f107,s)
mem = 5e8;      % adjust if needed

% calculate radii in units of reference radius

a = 6371.2;     % [km]
rho = r / a;
rho_Sq = (a + s.h) / a;

% extract data from s structure

nmax = s.nmax;
mmax = s.mmax;
p_vec = s.p_vec;
s_vec = s.s_vec;
theta_NGP = s.theta_NGP;
phi_NGP = s.phi_NGP;

% check inputs are scalars or 1D vectors

if ~isvector(rho) || ~isvector(theta) || ~isvector(phi) || ...
        ~isvector(t) || ~isvector(f107)
    error('Inputs must be scalars or 1D vector.');
end

% convert to column vectors

rho = rho(:);
theta = theta(:);
phi = phi(:);
t = t(:);
f107 = f107(:);

% expand scalars to match the maximum length of non-scalars

max_length = max([length(rho); length(theta); length(phi); ...
    length(t); length(f107)]);

if isscalar(rho),   rho   = rho   * ones(max_length, 1); end
if isscalar(theta), theta = theta * ones(max_length, 1); end
if isscalar(phi),   phi   = phi   * ones(max_length, 1); end
if isscalar(t),     t     = t     * ones(max_length, 1); end
if isscalar(f107),  f107  = f107  * ones(max_length, 1); end

if ~isequal(size(rho), size(theta), size(phi), size(t), size(f107))
    error('Variables must be of equal size (or scalars)');
end

% calculate time in year (season) and MUT

t_1 = datevec(t);
year = t_1(:,1) + 2000;
ndays = zeros(size(year));
for i = 1:length(year)
    ndays(i) = sum(eomday(year(i), 1:12), 2);
end
t_season = (t - jd2000(year, 1, 1, 0)) ./ ndays;

t_mut = getmut(t, theta_NGP, phi_NGP);

% calculate number of coefficients

N_nm       = mmax * (mmax + 2) + (nmax - mmax) * (2 * mmax + 1);
N_sp       = length(p_vec) * length(s_vec);
N_coeff_nm = 2 * N_nm * N_sp;

if (size(s.m_e_d_Re, 1) ~= N_nm || size(s.m_e_d_Re, 2) ~= 2 * N_sp)
    error('wrong number of model coefficients')
end

% convert coefficients to column vectors

m_e_d_Re = reshape(s.m_e_d_Re, N_coeff_nm, 1);
m_i_d_Re = reshape(s.m_i_d_Re, N_coeff_nm, 1);

% calculate number of data and data chunk size

N_data = length(rho);
chunk_size = fix(mem / (8 * 3 * 2 * N_coeff_nm));
% disp(['data chunk size: ', num2str(chunk_size)])

% preallocate intermediate output

B_r_1_tmp     = zeros(N_data, 1);
B_theta_1_tmp = zeros(N_data, 1);
B_phi_1_tmp   = zeros(N_data, 1);

B_r_2_tmp     = zeros(N_data, 1);
B_theta_2_tmp = zeros(N_data, 1);
B_phi_2_tmp   = zeros(N_data, 1);   

% CASE 1: above Sq currents

if (min(rho) > rho_Sq)

    % calculate primary Sq model above Sq currents

    f = zeros(N_nm, 1);     % diagonal of C matrix
    i1 = 1;
    for n = 1:nmax
        nr = min(2 * n + 1, 2 * mmax + 1);
        f(i1:i1 + nr - 1) = -n / (n + 1) * rho_Sq^(2 * n + 1);
        i1 = i1 + nr;
    end

    f = repmat(f, 2 * N_sp, 1);
    m_e_d_Re_C = bsxfun(@times, m_e_d_Re, f);

    % loop over data chunks
    
    for k = 1:chunk_size:N_data

        index = k:min(k+chunk_size-1, N_data);

        % calculate dipolar coordinates + matrix R
        
        [theta_d, phi_d, R] = gg2gm_v2(theta(index), phi(index), ...
            1, [], [], [theta_NGP phi_NGP]);   % [deg]

        % calculate design matrices

        [A_r_i_d, A_theta_i_dd, A_phi_i_dd] = ...
            design_SHA_Sq_i_Re_v2(rho(index), theta_d, phi_d, ...
            t_season(index), t_mut(index), nmax, mmax, p_vec, s_vec);
                
        r11 = R(:,1,1)';
        r21 = R(:,2,1)';
        r12 = R(:,1,2)';
        r22 = R(:,2,2)';
        
        A_theta_i_d = bsxfun(@times, r11, A_theta_i_dd) + ...
                  bsxfun(@times, r21, A_phi_i_dd);
        A_phi_i_d = bsxfun(@times, r12, A_theta_i_dd) + ...
                  bsxfun(@times, r22, A_phi_i_dd);

        % calculate magnetic fields

        B_r_1_tmp(index)     = A_r_i_d' * m_e_d_Re_C;
        B_theta_1_tmp(index) = A_theta_i_d' * m_e_d_Re_C;
        B_phi_1_tmp(index)   = A_phi_i_d' * m_e_d_Re_C;
    
        B_r_2_tmp(index)     = A_r_i_d' * m_i_d_Re;
        B_theta_2_tmp(index) = A_theta_i_d' * m_i_d_Re;
        B_phi_2_tmp(index)   = A_phi_i_d' * m_i_d_Re;   
        
    end

% CASE 2: below Sq currents
      
elseif (max(rho) < rho_Sq)
       
    % loop over data chunks

    for k = 1:chunk_size:N_data

        index = k:min(k+chunk_size-1, N_data);

        % calculate dipolar coordinates + matrix R
        
        [theta_d, phi_d, R] = gg2gm_v2(theta(index), phi(index), ...
            1, [], [], [theta_NGP phi_NGP]);   % [deg]

        % calculate design matrices

        [A_r_i_d, A_theta_i_dd, A_phi_i_dd] = ...
            design_SHA_Sq_i_Re_v2(rho(index), theta_d, phi_d, ...
            t_season(index), t_mut(index), nmax, mmax, p_vec, s_vec);
    
        [A_r_e_d, A_theta_e_dd, A_phi_e_dd] = ...
            design_SHA_Sq_e_Re_v2(rho(index), theta_d, phi_d, ...
            t_season(index), t_mut(index), nmax, mmax, p_vec, s_vec);
    
        r11 = R(:,1,1).';
        r21 = R(:,2,1).';
        r12 = R(:,1,2).';
        r22 = R(:,2,2).';
        
        A_theta_i_d = bsxfun(@times, r11, A_theta_i_dd) + ...
                  bsxfun(@times, r21, A_phi_i_dd);
        A_phi_i_d = bsxfun(@times, r12, A_theta_i_dd) + ...
                  bsxfun(@times, r22, A_phi_i_dd);

        A_theta_e_d = bsxfun(@times, r11, A_theta_e_dd) + ...
                  bsxfun(@times, r21, A_phi_e_dd);
        A_phi_e_d = bsxfun(@times, r12, A_theta_e_dd) + ...
                  bsxfun(@times, r22, A_phi_e_dd);

        % calculate magnetic field
    
        B_r_1_tmp(index)     = A_r_e_d' * m_e_d_Re;
        B_theta_1_tmp(index) = A_theta_e_d' * m_e_d_Re;
        B_phi_1_tmp(index)   = A_phi_e_d' * m_e_d_Re;
    
        B_r_2_tmp(index)     = A_r_i_d' * m_i_d_Re;
        B_theta_2_tmp(index) = A_theta_i_d' * m_i_d_Re;
        B_phi_2_tmp(index)   = A_phi_i_d' * m_i_d_Re;

    end
    
% CASE 3: error
    
else
    
    error('data in both regions (below and above Sq currents)')

end

% correct for F10.7 dependence

w = (1+s.N*f107);

B_1 = bsxfun(@times, [B_r_1_tmp B_theta_1_tmp B_phi_1_tmp], w);
B_2 = bsxfun(@times, [B_r_2_tmp B_theta_2_tmp B_phi_2_tmp], w);











