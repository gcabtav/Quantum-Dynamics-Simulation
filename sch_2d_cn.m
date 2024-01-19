function [x, y, t, psi, psire, psiim, psimod, prob, v] = sch_2d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)

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
    % y: Vector of y coordinates [ny]
    % t: Vector of t coordinates [nt]
    % psi: Array of computed psi values [nt x nx x ny]
    % psire: psi: Array of computed psi_re values [nt x nx x ny]
    % psiim: Array of computed psi_im values [nt x nx x ny]
    % psimod: Array of computed sqrt(psi psi*) values [nt x nx x ny]
    % prob: Array of computed running integral values [nt x nx x ny]
    % v: Array of potential values [nx x ny] (time-independent)

   % Define mesh and derived parameters ...
   nx = 2^level + 1;
   ny = nx;
   dx = 2^(-level);
   dy = dx;
   dt = lambda * dx;
   nt = round(tmax / dt) + 1;

   x = linspace(0.0, 1.0, nx);
   y = linspace(0.0, 1.0, ny);
   t = (0 : nt-1) * dt;

    % Initialize solution and set Initial data types:
    psi = zeros(nt, nx, ny);

    % idtype is an integer that chooses initial data type 
    if idtype == 0 % exact family 
        mx = idpar(1);
        my = idpar(2);
        psi(1, :, :) = sin(mx*pi*x).*sin(my*pi*x).';
    elseif idtype == 1 % boosted gaussian
        x0 = idpar(1);
        y0 = idpar(2);
        delta_x = idpar(3);
        delta_y = idpar(4);
        px = idpar(5);
        py = idpar(6);
        psi(1, :, :) = exp(1i*px*x).*exp(-((x-x0)/delta_x).^2).*(exp(1i*py*y).*exp(-((y-y0)/delta_y).^2)).';
    else
      fprintf('sch_2d_cn: Invalid idtype=%d\n', idtype);
      return
    end

    % vtype is an integer that chooses initial potential type 
    v = zeros(nx, ny);
    if vtype == 0 % no potential
        v = zeros(nx, ny);
    elseif vtype == 1 % rectangular barrier or well
        xmin = vpar(1);
        xmax = vpar(2);
        ymin = vpar(3);
        ymax = vpar(4);
        vc = vpar(5);
        for i = 1 : nx
            for j = 1 : ny
                if (x(i) < xmin || x(i) > xmax) || (y(j) < ymin || y(j) > ymax)
                    v(i, j) = 0;
                else
                    v(i, j) = vc;
                end
            end
        end
    elseif vtype == 2 % double slit
        x1 = vpar(1);
        x2 = vpar(2);
        x3 = vpar(3);
        x4 = vpar(4);
        vc = vpar(5);
        j_p = (ny - 1)/4 + 1;
        for i = 1 : nx
            if ((x1 <= x(i)) && (x(i) <= x2)) || ((x3 <= x(i)) && (x(i) <= x4)) 
                v(i, j_p) = 0;
                v(i, j_p + 1) = 0;
            else
                v(i, j_p) = vc;
                v(i, j_p + 1) = vc;
            end
        end
    elseif vtype == 3 % two double slit barriers
        x1 = vpar(1);
        x2 = vpar(2);
        x3 = vpar(3);
        x4 = vpar(4);
        vc = vpar(5);
        j_p = (ny - 1) / 4 + 1;
        j_p1 = 2*j_p; 

        % Create masks for the double slit condition
        mask1 = (x >= x1) & (x <= x2);
        mask2 = (x >= x3) & (x <= x4);

        % Define the potential matrix based on masks
        v(:, j_p:j_p+1) = vc;
        v(:, j_p1:j_p1+1) = vc;
        v(mask1 | mask2, j_p:j_p+1) = 0;
        v(mask1 | mask2, j_p1:j_p1+1) = 0;
    elseif vtype == 4 % Circular barrier or well
        xmin = vpar(1);
        xmax = vpar(2);
        ymin = vpar(3);
        ymax = vpar(4);
        vc = vpar(5);
        x_center = (xmin + xmax) / 2;
        y_center = (ymin + ymax) / 2;
        radius = min((xmax - xmin), (ymax - ymin)) / 2; % Using minimum side length as the radius

        % Create a meshgrid around the circular region
        [X, Y] = meshgrid(x, y);

        % Calculate the distance from the center for each point in the meshgrid
        distances = sqrt((X - x_center).^2 + (Y - y_center).^2);

        % Define the potential based on the circular shape
        v(distances <= radius) = vc;
        v(distances > radius) = 0;   
    elseif vtype == 5 % Gaussian potential well or barrier
        xmin = vpar(1);
        xmax = vpar(2);
        ymin = vpar(3);
        ymax = vpar(4);
        vc = vpar(5);
        x_center = (xmin + xmax) / 2;
        y_center = (ymin + ymax) / 2;
        sigma_x = (xmax - xmin) / 2; % Adjust the sigma values as needed for the desired width
        sigma_y = (ymax - ymin) / 2;
    
        % Create a meshgrid around the defined region
        [X, Y] = meshgrid(x, y);
    
        % Define the Gaussian potential based on the meshgrid
        v = vc * exp(-(X - x_center).^2 / (2 * sigma_x^2) - (Y - y_center).^2 / (2 * sigma_y^2));
    else
      fprintf('sch_2d_cn: Invalid vtype=%d\n', vtype);
      return
    end


   % Initialize storage 

   f  = zeros(nx, ny);
   X  = zeros(nx, ny);
   psi_nhalf = zeros(nx, ny);

   % Tridiagonal systems

   % A
   dl = (-1i*dt*(dx^-2)/2) * ones(nx, 1);
   d = (1 + 1i*dt*dx^(-2)) * ones(nx, 1);
   du = dl;

   % B 
   dlB = -dl;
   dB = (1 - 1i*dt*dx^(-2)) * ones(nx, 1);
   duB = -dl;
   
   % C
   % Initializing an empty cell array to store sparse matrices
   Cs = cell(1, nx);
   
   dlC = (1i*dt*(dy^-2)/2) * ones(nx, 1);
   duC = dlC;
   dlC(nx-1) = 0;
   duC(2) = 0;
   for i = 1 : nx
       dC = (1 - 1i*dt*dy^(-2) - 1i*dt*(v(i,:).')/2); 
       dC(1) = 1;
       dC(nx) = 1;     
       Cs{i} = spdiags([dlC dC duC], [-1 0 1], nx, nx);
   end

    % D   
    Ds = cell(1, nx);
    
    dlD = (-1i*dt*(dy^-2)/2) * ones(nx, 1);
    duD = dlD;
    dlD(nx-1) = 0;
    duD(2) = 0;
    
    for i = 1 : nx  
        dD = (1 + 1i*dt*dy^(-2) + 1i*dt*(v(i,:).')/2); 
        dD(1) = 1;
        dD(nx) = 1; 
        Ds{i} = spdiags([dlD dD duD], [-1 0 1], nx, nx);
    end


   % Fix up boundary cases
   dl(nx-1) = 0;
   d(1) = 1;
   d(nx) = 1;
   du(2) = 0;

   dlB(nx-1) = 0;
   dB(1) = 1;
   dB(nx) = 1;
   duB(2) = 0;

   % Define sparce matrix
   A = spdiags([dl d du], [-1 0 1], nx, nx);
   B = spdiags([dlB dB duB], [-1 0 1], nx, nx);  

   % Computing solution using ADI Scheme 
    for n = 1 : nt - 1
        for i = 1 : nx 
            X(i, :) = reshape(Cs{i} * squeeze(psi(n, i, :)), 1, []);
        end
        
        for j = 1 : ny
            f(:, j) = B*X(:, j); 
        end

        % Boundary Conditions
        f(1,:) = 0;     
        f(:, 1) = 0;
        f(nx, :) = 0;
        f(:, ny) = 0; 
        % f(1,1) = 0;     
        % f(nx, 1) = 0;
        % f(1, nx) = 0;
        % f(nx, ny) = 0; 
        
        for j = 1 : ny
            psi_nhalf(:, j) = A \ f(:, j);
        end
        
        for i = 1 : nx
            psi(n+1, i, :) = Ds{i} \ psi_nhalf(i, :).';
        end
    end    

    psire = real(psi);
    psiim = imag(psi);
    psimod = abs(psi);
    psisquared = psimod.^2;


    % Calculating probability using trapesoidal formula
    prob = zeros(nt, nx, ny);
    for n = 1 : nt
        for i = 1 : nx - 1
            for j = 1 : ny - 1
                prob(n, i+1, j+1) = prob(n, i, j) + (1/2)*(psisquared(n, i, j) + psisquared(n, i+1, j+1))*(x(i+1) - x(i))*(y(j+1) - y(j));
            end
        end
    end    
end
    