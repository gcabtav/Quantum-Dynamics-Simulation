function [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)

% Inputs
% 
    % tmax: Maximums integration time
    % level: Discretization level
    % lambda: dt/dx
    % idtype: Selects initial data type
    % idpar: Vector of initial data parameters
    % vtype: Selects of potential type
    % vpar: Vector of potential parameters
% 
% Outputs 
% 
    % x: Vector of x coordinates [nx]
    % t: Vector of t coordinates [nt]
    % psi: Array of computed psi values [nt x nx]
    % psire: psi: Array of computed psi_re values [nt x nx]
    % psiim: Array of computed psi_im values [nt x nx]
    % psimod: Array of computed sqrt(psi psi*) values [nt x nx]
    % prob: Array of computed running integral values [nt x nx]
    % v: Array of potential values [nx] (time-independent)

   % Define mesh and derived parameters
   nx = 2^level + 1;
   dx = 2^(-level);
   dt = lambda * dx;
   nt = round(tmax / dt) + 1;

   x = linspace(0.0, 1.0, nx);
   t = (0 : nt-1) * dt;

    % Initialize solution and set Initial data types:
    psi = zeros(nt, nx);

    % idtype is an integer that chooses initial data type 
    if idtype == 0 % exact family 
        m = idpar(1);
        psi(1, :) = sin(m*pi*x);
    elseif idtype == 1 % boosted gaussian
        x0 = idpar(1);
        delta = idpar(2);
        p = idpar(3);
        psi(1, :) = exp(1i*p*x).*exp(-((x-x0)/delta).^2);
    else
      fprintf('sch_1d_cn: Invalid idtype=%d\n', idtype);
      return
    end

    % vtype is an integer that chooses initial potential type 
    v = zeros(nx);
    if vtype == 0 % no potential
        v = zeros(nx, 1);
    elseif vtype == 1 % rectangular barrier or well
        xmin = vpar(1);
        xmax = vpar(2);
        vc = vpar(3);
        for i = 1 : nx
            if x(i) < xmin || x(i) > xmax 
                v(i) = 0;
            else
                v(i) = vc;
            end
        end
    end  

   % MOL: Initialize storage for sparse matrix and RHS

   % Set up Tridiagonal systems
   dl = ((dx^-2)/2) * ones(nx, 1);
   for j = 1 : nx
    d(j, 1) = (1i*dt^-1 - dx^-2 - (v(j))/2);
   end
   du = dl;
    
   dl1 = -dl;
   for j = 1 : nx
    d1(j, 1) = (1i*dt^-1 + dx^-2 + (v(j))/2);
   end
   du1 = dl1;

   % Fix up boundary cases
   dl(nx-1) = 0;
   d(1) = 1;
   d(nx) = 1;
   du(2) = 0;

   dl1(nx-1) = 0;
   d1(1) = 1;
   d1(nx) = 1;
   du1(2) = 0;

   % Define sparce matrix
   A = spdiags([dl d du], [-1 0 1], nx, nx);
   B = spdiags([dl1 d1 du1], [-1 0 1], nx, nx);

   % Compute solution using Crank-Nicolson Scheme
    for n = 1 : nt - 1
        % Define RHS of linear system ...
        f = B*psi(n, :).';
        f(1) = 0;
        f(nx) = 0;
        % Solve system, thus updating approximation to next time-step
        psi(n+1, :) = A \ f;
    end

    psire = real(psi);
    psiim = imag(psi);
    psimod = abs(psi);
    psisquared = psimod.^2;
    

    % Calculating probability using trapesoidal formula (O(h^2))
    prob = zeros(nt, nx);
    for n = 1 : nt
        for j = 1 : nx - 1
            prob(n, j+1) = prob(n, j) + (1/2)*(psisquared(n, j) + psisquared(n, j + 1))*(x(j+1) - x(j));
        end
    end 

end





