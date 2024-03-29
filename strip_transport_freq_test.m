% script by Simon Otten (s.j.otten@utwente.nl, simonotten@fastmail.com)
% changed 14 June 2019
% use Matlab r2016b or newer

%Daten genommen von https://htsdb.wimbush.eu/dataset/3759327
clearvars;

%%optionen für das Programm



currentSweep = false;   % Einstellung ob der Strom von 0 bis Ic gesweept werden soll
percentage = 0.25;       % Wenn currentSweep = false: Eingabe für bei welchen I0 = x*Ic untersucht werden soll
heatmap = 1;            % option für Darstellung; 1 für Heatmap, 0 für Standard       
freqSweep = true;     % option für Sweep der Frequenz; true für an, false für aus. Bei aus wird baseFreq als Frequenz benutzt.
baseFreq = 10;           % Startfrequenz für Sweep [Hz], alternativ Frequenz ohne Sweep
endFreq = 1000;          % Endfrequenz für Sweep [Hz]
freqSteps = 50;          % Schritte für den Frequenzsweep

%% strip properties
W=0.012; % Erhältliche Breiten des zu untersuchenden Leiters, Manuell eintragen
Ec = 1e-4;          % electric field at J=Jc [V/m]
N = 100;            % number of elements for numerical calculation

% Konstanten Definieren
mu0 = 4e-7*pi;
% Einlesen und Processing von externen Dateien, Definition von Variablen


Jc = 25e3;     % vektor mit der gleichen länge wie n_var erstellen, damit code funktioniert
n = [1,5,15,50];
Ic = Jc*W;     %Kritischen Strom für momentane Temperatur setzen

gamma = 100;        % feedback constant for voltage source [V/m/A]
N_step = 1000;      % time steps per cycle




% Setzen der Frequenzvariable
if freqSweep == true
    f = linspace(baseFreq,endFreq,freqSteps);
else
    f = baseFreq;
end


q_values = zeros(size(n,2),size(f,2));% Matrix für Abspeicherung der Werte definieren
q_norris_values = zeros(size(q_values));

if currentSweep == true
    I0_vector = linspace(Ic/30,Ic,30);
else
    I0_vector = Ic*percentage;
end

% define left and right boundaries of the elements
%y_vector = linspace(-(W(count3))/2,(W(count2))/2,N+1)'; % option 1: linear spacing
y_vector = -(W)/2*cos((0:N)/N*pi)'; % option 2: smaller elements near edges (better for low amplitudes)
y_left = y_vector(1:end-1);     % left boundaries of each element
y_right = y_vector(2:end);      % right boundaries
y_middle = (y_left+y_right)/2;  % middle point
w = y_right - y_left;           % element width




%%% END OF SETTINGS %%%

for count_n = 1:length(n)

for count_freq = 1:length(f)   % for loop für Frequenzsweep

    % Definition von Zeitvariablen
    B0 = 0.0;
    omega = 2*pi*f(count_freq);
    dt = 1/(N_step*f(count_freq));
    t = 0:dt:2/f(count_freq);
    Bext = @(t) B0*sin(omega*t);
    Bdot = @(t) omega*B0*cos(omega*t);
    

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

% create vector for the AC loss results VIELLEICHT LÖSCHEN
Q_vector = zeros(size(I0_vector)); % AC loss per cycle [J/m]

for i = 1:numel(I0_vector)
    % applied current
    I0 = I0_vector(i);    
    Iset = @(t) I0*sin(omega*t);
    
    % electric field function and derivative
    E = @(J) Ec*sign(J).*abs(J./Jc).^n(count_n);
    dEdJ = @(J) Ec*(n(count_n)./Jc).*abs(J./Jc).^((n(count_n))-1);
    
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

    %bei erster J matrix den größten wert ermitteln, dies dann als limits
    %für y-achse verwenden.
    if count_freq==1 && i==1
        maxy = 1.25*max(J, [], 'all');
    end

    %anfang für heatmap?
    if heatmap == 1
        maxMatrix = max(J, [], 1); %jeweils höchste werte der matrix auslesen
        figure(104)
        clf
        area(y_middle,maxMatrix)
        %ylim([0 maxy]) % NICHT VERGESSEN NOCHMAL ANZUGUCKEN
        title(['Frequenz = ',num2str(f(count_freq)),' Hz, Breite = ',num2str(1000*W),' mm, Ic = ',num2str(Ic),' A'])
        subtitle(['I0 = ',num2str(I0_vector(i)),' A'])
        xlabel('Position im Leiter [m]')
        ylabel('J [A/m]')
        drawnow
        
    else
    %figure with current distribution
        figure(100)
        clf
        plot(y_middle,J(N_step:10:2*N_step,:)') %N_step teil sind anzahl der linien, : ist 100 schritte,, falls änderung der darstellung nicht möglich, vielleicht in der arbeit stattdessen erklären, was der plot bedeutet?
        ylim([-maxy maxy]);                       %y-achse festgesetzt, damit man relativ die Stromstärke erkennen kann
        title(['I_0 = ',num2str(I0),' A, Q = ',num2str(Q),' J/m'])
        xlabel('y [m]')
        ylabel('J [A/m]')
        drawnow
    end

    if currentSweep == true
        q_values(count_freq) = Q_vector(24);     %Fallback auf Verlust bei I0 = 80%Ic 24, da Ic in 30 aufgeteilt ist
    else
        %norris ist nicht von frequenz abhängig, also vielleicht ändern,
        %oder nur für keinen Sweep benutzen
        q_values(count_n,count_freq) = Q_vector(i); 
        q_norris_values(count_n,count_freq) = mu0*Ic^2/pi * ( ...
            (1-I0_vector/Ic).*log(1-I0_vector/Ic) + ...
            (1+I0_vector/Ic).*log(1+I0_vector/Ic) - ...
            (I0_vector/Ic).^2 ...
        );
    end


end
end
end

%% Berechnen und Darstellung der Differenz
norris_sym_diff = q_values - q_norris_values; % Differenz zwischen simulierten und berechneten Werten als Vektor


figure(504)
clf
plot(f,norris_sym_diff)
hold on
title('Differenz Zwischen der simulierten Werte und Norris')
subtitle({'Über der Referenzlinie sind die Verluste der simulierten Werte größer, unterhalb die von Norris', ...
    ['I_0/I_C = ',num2str(percentage)]})
xlabel('Frequenz in Hz')
ylabel('Numeric - Norris in J/m')
yline(0,'-','Referenz')
legend('n=1','n=5','n=15','n=50')

