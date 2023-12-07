

params = [4, 4, 12, 10, 100, 0];
%fuel_mass_4series(params)

min_guess = [0, 0, 11, 3, 50, 0];
max_guess = [9, 9, 30, 50, 250, 1];

a = fmincon(@fuel_mass_4series, params, [], [], [], [], min_guess, max_guess, @constraint)
%a = fmincon(@fuel_mass_4series, params, [], [], [], [], min_guess, max_guess)

function E_co2 = fuel_mass_4series(v)
    global CL
    global LpD
    global rho_inf
    global w0
    global velocity
    global NACA34
   
    NACA1 = round(round(v(1)));
    NACA2 = round(round(v(2)));
    NACA34 = round(round(v(3)));
    c = v(4);
    b = v(5);
    fuel_frac = v(6);

    FCF = 3.16; % Fuel consumption factor
    LC = 89; % baseline emissions
    AR = b/c; % aspect ratio
    rho_inf = 0.4135; % density of air at 35,000 ft in kg/m^3
    g = 9.81; % gravity

    range_act = 5837;%*1000;
    cruise_alt = 35000;
    cruise_mach = .65;
    n_t = 43.10;
    wb = 170000*4.4482216153;
    span_e = .9; % span efficiency factor

    alpha = 0:10;
    speed_sound = 295; % m/s
    velocity = cruise_mach * speed_sound * 3.28084;
    viscosity = 1.458e-5;
    viscosity = 2.995e-7;
    rho_inf = 7.38e-4;
    Re = rho_inf * velocity * c / viscosity % Reynolds number
    %Re = 1e6;

    airfoil = append('NACA', int2str(NACA1), int2str(NACA2), sprintf('%02d', NACA34));
    [pol foil] = xfoil(airfoil,alpha,Re,cruise_mach,'ppar T 0.2', 'ppar n 80','oper iter 2000');

    n = (span_e * AR) / (span_e * AR + 2);

    CL = pol.CL * n;
    CD = pol.CD + CL.^2 / (pi * span_e * AR);

    [LpD, i] = max(CL./CD);
    CL = CL(i);
    CD = CD(i);
    %[LpD, i] = max(pol.CL./pol.CD)

    Sw = b * c;

    thickness = NACA34 / 100;
    W_wing = (228 * AR * (Sw*10.764)^(0.5) * (thickness)^(-0.04)) * 4.4482216153;

    w1 = wb + W_wing;
    d_hf = (42.8 - 20 * fuel_frac^2);
    w0 = exp(g*range_act / (d_hf*n_t * LpD)) * w1;

    W_fuel = 1.06 * (w0 - w1);
    M_fuel = W_fuel/g;

    L = CL * 0.5 * rho_inf * (velocity)^2 * Sw;
    LS_f = (89 - 30 * fuel_frac); % lifetime emissions

    range = d_hf/9.81 * n_t * LpD * log(w0/w1);
    E_co2 = M_fuel * FCF *(LS_f / LC);

    airfoil
    fuel_frac
    E_co2
    range
    w1
    w0
    L
    b
    c

end

function [c,ceq] = constraint(v)
    global CL
    global LpD
    global rho_inf
    global w0
    global velocity
    cc = v(4);
    b = v(5);
    fuel_frac = v(6);

    L = CL * 0.5 * rho_inf * (velocity)^2 * b * cc;

    c = [];
    ceq = CL * 0.5 * rho_inf * (velocity)^2 * b * cc - w0

end