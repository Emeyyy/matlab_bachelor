% script by Simon Otten (s.j.otten@utwente.nl, simonotten@fastmail.com)
% changed 14 June 2019
% use Matlab r2016b or newer

% strip properties
W = 0.012;          % strip width [m]
Jc = 25e3;          % critical current per unit width [A/m]
Ic = Jc*W;          % critical current [A]
Ec = 1e-4;          % electric field at J=Jc [V/m]
n = 100;             % n-value
N = 100;            % number of elements for numerical calculation

% external magnetic field (perpendicular) and applied current
f = 1;              % frequency [Hz]
B0 = 0.0;           % amplitude [T]
omega = 2*pi*f;
Bext = @(t) B0*sin(omega*t);
Bdot = @(t) omega*B0*cos(omega*t); % time derivative

% applied current
I0_vector = linspace(Ic/30,2*Ic,30); % amplitude vector
gamma = 100;        % feedback constant for voltage source [V/m/A]

% time vector
N_step = 1000;      % time steps per cycle
dt = 1/(N_step*f);  % time step [s]
t = 0:dt:2/f;       % time vector (2 cycles)

%%% END OF SETTINGS %%%

% define left and right boundaries of the elements
%y_vector = linspace(-W/2,W/2,N+1)'; % option 1: linear spacing
y_vector = -W/2*cos((0:N)/N*pi)'; % option 2: smaller elements near edges (better for low amplitudes)
y_left = y_vector(1:end-1);     % left boundaries of each element
y_right = y_vector(2:end);      % right boundaries
y_middle = (y_left+y_right)/2;  % middle point
w = y_right - y_left;           % element width

% calculate the K matrix (used as "mass matrix" for the ode solver)
K = zeros(N);
for i = 1:N
    yi = y_middle(i);
    for j = 1:N
        aj = y_left(j);
        bj = y_right(j);
        K(i,j) = (yi-bj)*(log(abs(yi-bj))-1) -(yi-aj)*(log(abs(yi-aj))-1);
    end
end

% create vector for the AC loss results
Q_vector = zeros(size(I0_vector)); % AC loss per cycle [J/m]

for i = 1:numel(I0_vector)
    % applied current
    I0 = I0_vector(i);    
    Iset = @(t) I0*sin(omega*t);
    
    % electric field function and derivative
    E = @(J) Ec*sign(J).*abs(J./Jc).^n;
    dEdJ = @(J) Ec*(n./Jc).*abs(J./Jc).^(n-1);
    
    % initial current distribution
    J0 = zeros(N,1);
    
    % gradient of the electric potential (voltage source)
    dphi = @(t,J) gamma*(w'*J-Iset(t));
    
    % function to integrate (see equation 4)
    % the factor 5e6 equals 2pi/mu0
    fun = @(t,J) -5e6 * (E(J)-y_middle*Bdot(t)+dphi(t,J));
    
    % jacobian matrix
    jac = @(t,J) -5e6 * (diag(dEdJ(J))+gamma*w');
    
    % numerical integration
    options = odeset();
    options.AbsTol = 1e-3*Jc;
    options.Jacobian = jac;
    options.Mass = K;
    options.MStateDependence = 'none';
    options.RelTol = 1e-6;
    options.Stats = 'on';
    options.Vectorized = 'on';
    tic;
    [T,J] = ode15s(fun,t,J0,options);
    toc
    
    %%% POST-PROCESSING %%%
    
    % calculate the dissipated power at each time step
    P = J.*E(J)*w;
    % calculate the AC loss per cycle by integrating the dissipated power over
    % the last cycle
    Q = trapz(t(end-N_step:end),P(end-N_step:end));
    disp(['Q = ',num2str(Q),' J/m']); disp(' ')
    Q_vector(i) = Q;
    
    % figure with current distribution
    figure(100)
    clf
    plot(y_middle,J(N_step:10:2*N_step,:)')
    title(['I_0 = ',num2str(I0),' A, Q = ',num2str(Q),' J/m'])
    xlabel('y [m]')
    ylabel('J [A/m]')
    drawnow         
end

% plot loss per cycle
figure(101)
clf
semilogy(I0_vector,Q_vector)
hold on
% compare with exact result from Norris
mu0 = 4e-7*pi;
Q_norris = mu0*Ic^2/pi * ( ...
    (1-I0_vector/Ic).*log(1-I0_vector/Ic) + ...
    (1+I0_vector/Ic).*log(1+I0_vector/Ic) - ...
    (I0_vector/Ic).^2 ...
);
plot(I0_vector,Q_norris)
title('Transport AC loss')
xlabel('I_0 [A]')
ylabel('Q [J/m]')
legend('Numeric','Norris')