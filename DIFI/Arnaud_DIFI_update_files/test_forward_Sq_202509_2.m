% test_forward_Sq_202509_2.m
%
% Test the output and speed of the "forward_Sq_d_Re_v2.m" DIFI forward 
% calculation function.
%
% A. Chulliat, 2025-09-28
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

addpath('../')

rad = pi/180;

a = 6371.2;            % reference radius [km]
h = 110;               % height of Sq currents [km]

% filenames

path_model = ...
    '../../Runs_xDIFI_2024/2024-04-03_1044_DIFI_r9_TF_COR_AC_20240402_THRUD';

filename_model = fullfile(path_model, 'PRODUCT', ...
    'SW_TEST_MIO_SHA_2D_20140101T000010_20231231T235700_1001_TEST.DBL');

path_F107 = '../../F107_for_DIFI/data_F107';

filename_f107 = fullfile(path_F107, ...
    'SW_OPER_AUX_F10_2__20060101T000000_20250815T000000_0001.DBL');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test points #1: altitude profile below currents at noon, April 15, 2022

r_1     = a + (0:10:100)';              % assuming a spherical Earth
theta_1 = 45;                           % [deg]
phi_1   = 0;                            % [deg]
t_1     = jd2000(2022, 4, 15, 12);      % [MJD2000]

% test points #2: altitude profile above currents at noon, April 15, 2022

r_2     = a + (120:10:500)';            % assuming a spherical Earth
theta_2 = 45;                           % [deg]
phi_2   = 0;                            % [deg]
t_2     = jd2000(2022, 4, 15, 12);      % [MJD2000]

% test points #3: time variation at BOU observatory on April 10-20, 2022

lat_3 = 90-49.86;               % [deg]
h_3   = 2; %1.682;                  % [km]

[r_3, theta_3] = geod2geoc(lat_3*rad, h_3);                     % [km, rad]
theta_3 = theta_3/rad;                                          % [deg]
phi_3   = 254.76;                                               % [deg]
t_3     = (jd2000(2022,4,10,0):1/24:jd2000(2022,4,20,23))';     % [MJD2000]

% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read model

s = SwarmL2_MIO_SHA_Read_v2(filename_model);

disp(s)

% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read F10.7 data

[t_f107, f107] = SwarmL2_F107_Read(filename_f107);

% figure
% hold on
% grid on
% plot(t_f107, f107, '.-b')
% title('Daily F10.7 from 2006-01-01 to 2025-08-15')

f107_1 = interp1(t_f107, f107, t_1);
f107_2 = interp1(t_f107, f107, t_2);
f107_3 = interp1(t_f107, f107, t_3);

% [~, index_1] = min(abs(t_f107 - t_1));
% [~, index_2] = min(abs(t_f107 - t_2));
% index_3 = find(t_f107 >= t_3(1) & t_f107 <= t_3(end));

figure
hold on
grid on
% plot(t_f107(index_1), f107(index_1), 'o-r')
% plot(t_f107(index_2), f107(index_1), '*-g')
% plot(t_f107(index_3), f107(index_3), '.-b')
plot(t_1, f107_1, 'o-r')
plot(t_2, f107_2, '*-g')
plot(t_3, f107_3, '.-b')
legend('test points #1', 'test points #2', 'test points #3', ...
    'Location', 'NorthEast')
title('Daily F10.7 for selected test points')

% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate Sq field for test points #1 & #2

[B_1_1, B_1_2] = forward_Sq_d_Re(r_1, theta_1, phi_1, t_1, ...
    f107_1, s);
[B_2_1, B_2_2] = forward_Sq_d_Re(r_2, theta_2, phi_2, t_2, ...
    f107_2, s);

[B_1_1_v2, B_1_2_v2] = forward_Sq_d_Re_v2(r_1, theta_1, phi_1, t_1, ...
    f107_1, s);
[B_2_1_v2, B_2_2_v2] = forward_Sq_d_Re_v2(r_2, theta_2, phi_2, t_2, ...
    f107_2, s);

[B_1_1_v2_par, B_1_2_v2_par] = parforward_Sq_d_Re_v2(r_1, theta_1, ...
    phi_1, t_1, f107_1, s);
[B_2_1_v2_par, B_2_2_v2_par] = parforward_Sq_d_Re_v2(r_2, theta_2, ...
    phi_2, t_2, f107_2, s);

% control plots

figure
%
subplot(3,2,1)
hold on
grid on
plot(r_1 - a, B_1_1(:, 1), '.-b');
plot(r_1 - a, B_1_1_v2(:, 1), 'dg');
plot(r_1 - a, B_1_1_v2_par(:, 1), '.k');
plot(r_2 - a, B_2_1(:, 1), '.-b')
plot(r_2 - a, B_2_1_v2(:, 1), 'dg')
plot(r_2 - a, B_2_1_v2_par(:, 1), '.k')
legend('v1', 'v2', 'v2-par', 'Location', 'West')
xlabel('h [km]')
ylabel('B_r [nT]')
title 'primary field'
%
subplot(3,2,2)
hold on
grid on
plot(r_1 - a, B_1_2(:, 1), '.-b');
plot(r_1 - a, B_1_2_v2(:, 1), 'dg');
plot(r_1 - a, B_1_2_v2_par(:, 1), '.k');
plot(r_2 - a, B_2_2(:, 1), '.-b')
plot(r_2 - a, B_2_2_v2(:, 1), 'dg')
plot(r_2 - a, B_2_2_v2_par(:, 1), '.k')
xlabel('h [km]')
ylabel('B_r [nT]')
title 'secondary field'
%
subplot(3,2,3)
hold on
grid on
plot(r_1 - a, B_1_1(:, 2), '.-b')
plot(r_1 - a, B_1_1_v2(:, 2), 'dg')
plot(r_1 - a, B_1_1_v2_par(:, 2), '.k')
plot(r_2 - a, B_2_1(:, 2), '.-b')
plot(r_2 - a, B_2_1_v2(:, 2), 'dg')
plot(r_2 - a, B_2_1_v2_par(:, 2), '.k')
xlabel('h [km]')
ylabel('B_\theta [nT]')
title 'primary field'
%
subplot(3,2,4)
hold on
grid on
plot(r_1 - a, B_1_2(:, 2), '.-b')
plot(r_1 - a, B_1_2_v2(:, 2), 'dg')
plot(r_1 - a, B_1_2_v2_par(:, 2), '.k')
plot(r_2 - a, B_2_2(:, 2), '.-b')
plot(r_2 - a, B_2_2_v2(:, 2), 'dg')
plot(r_2 - a, B_2_2_v2_par(:, 2), '.k')
xlabel('h [km]')
ylabel('B_\theta [nT]')
title 'secondary field'
%
subplot(3,2,5)
hold on
grid on
plot(r_1 - a, B_1_1(:, 3), '.-b')
plot(r_1 - a, B_1_1_v2(:, 3), 'dg')
plot(r_1 - a, B_1_1_v2_par(:, 3), '.k')
plot(r_2 - a, B_2_1(:, 3), '.-b')
plot(r_2 - a, B_2_1_v2(:, 3), 'dg')
plot(r_2 - a, B_2_1_v2_par(:, 3), '.k')
xlabel('h [km]')
ylabel('B_\phi [nT]')
title 'primary field'
%
subplot(3,2,6)
hold on
grid on
plot(r_1 - a, B_1_2(:, 3), '.-b')
plot(r_1 - a, B_1_2_v2(:, 3), 'dg')
plot(r_1 - a, B_1_2_v2_par(:, 3), '.k')
plot(r_2 - a, B_2_2(:, 3), '.-b')
plot(r_2 - a, B_2_2_v2(:, 3), 'dg')
plot(r_2 - a, B_2_2_v2_par(:, 3), '.k')
xlabel('h [km]')
ylabel('B_\phi [nT]')
title 'secondary field'

% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate Sq field for test points #3

[B_3_1, B_3_2] = forward_Sq_d_Re(r_3, theta_3, phi_3, t_3, ...
    f107_3, s);

[B_3_1_v2, B_3_2_v2] = forward_Sq_d_Re_v2(r_3, theta_3, phi_3, t_3, ...
    f107_3, s);

[B_3_1_v2_par, B_3_2_v2_par] = parforward_Sq_d_Re_v2(r_3, theta_3, ...
    phi_3, t_3, f107_3, s);

% control plots

figure
%
subplot(3,2,1)
hold on
grid on
plot(t_3, B_3_1(:,1), '.-b')
plot(t_3, B_3_1_v2(:,1), 'og')
plot(t_3, B_3_1_v2_par(:,1), '.k')
legend('v1', 'v2', 'v2-par', 'Location', 'NorthWest')
xlabel 'time [MJD2000]'
ylabel 'B_r [nT]'
title 'primary field'
%
subplot(3,2,2)
hold on
grid on
plot(t_3, B_3_2(:,1), '.-b')
plot(t_3, B_3_2_v2(:,1), 'og')
plot(t_3, B_3_2_v2_par(:,1), '.k')
xlabel 'time [MJD2000]'
ylabel 'B_r [nT]'
title 'secondary field'
%
subplot(3,2,3)
hold on
grid on
plot(t_3, B_3_1(:,2), '.-b')
plot(t_3, B_3_1_v2(:,2), 'og')
plot(t_3, B_3_1_v2_par(:,2), '.k')
xlabel 'time [MJD2000]'
ylabel 'B_\theta [nT]'
title 'primary field'
%
subplot(3,2,4)
hold on
grid on
plot(t_3, B_3_2(:,2), '.-b')
plot(t_3, B_3_2_v2(:,2), 'og')
plot(t_3, B_3_2_v2_par(:,2), '.k')
xlabel 'time [MJD2000]'
ylabel 'B_r\theta [nT]'
title 'secondary field'
%
subplot(3,2,5)
hold on
grid on
plot(t_3, B_3_1(:,3), '.-b')
plot(t_3, B_3_1_v2(:,3), 'og')
plot(t_3, B_3_1_v2_par(:,3), '.k')
xlabel 'time [MJD2000]'
ylabel 'B_\phi [nT]'
title 'primary field'
%
subplot(3,2,6)
hold on
grid on
plot(t_3, B_3_2(:,3), '.-b')
plot(t_3, B_3_2_v2(:,3), 'og')
plot(t_3, B_3_2_v2_par(:,3), '.k')
xlabel 'time [MJD2000]'
ylabel 'B_\phi [nT]'
title 'secondary field'

% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% performance test

tic
[B_1_1_test, B_1_2_test] = forward_Sq_d_Re(...
    repmat(r_1(1), 5e3, 1), ...
    theta_1, phi_1, t_1, f107_1, s);
toc

tic
[B_2_1_test, B_2_2_test] = forward_Sq_d_Re(...
    repmat(r_2(1), 5e3, 1), ...
    theta_2, phi_2, t_2, f107_2, s);
toc

tic
[B_1_1_v2_test, B_1_2_v2_test] = forward_Sq_d_Re_v2(...
    repmat(r_1(1), 5e3, 1), ...
    theta_1, phi_1, t_1, f107_1, s);
toc

tic
[B_2_1_v2_test, B_2_2_v2_test] = forward_Sq_d_Re_v2(...
    repmat(r_2(1), 5e3, 1), ...
    theta_2, phi_2, t_2, f107_2, s);
toc

tic
[B_1_1_v2_par_test, B_1_2_v2_par_test] = parforward_Sq_d_Re_v2( ...
    repmat(r_1(1), 5e3, 1), ...
    theta_1, phi_1, t_1, f107_1, s);
toc

tic
[B_2_1_v2_par_test, B_2_2_v2_par_test] = parforward_Sq_d_Re_v2( ...
    repmat(r_2(1), 5e3, 1), ...
    theta_2, phi_2, t_2, f107_2, s);
toc
