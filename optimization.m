

params = [2, 4, 12, 10, 120, .3];
% fuel_mass_4series(params)

min_guess = [0, 0, 11, 5, 50, 0];
max_guess = [9, 9, 30, 50, 250, 1];

a = fmincon(@fuel_mass_4series, params, [], [], [], [], min_guess, max_guess, @constraint)

function E_co2 = fuel_mass_4series(v)
    global cl
    global cd
    global NACA34
   
    NACA1 = round(v(1));
    NACA2 = round(v(2));
    NACA34 = round(v(3));
    c = v(4);
    b = v(5);
    fuel_frac = v(6);

    FCF = 3.16; % Fuel consumption factor
    LS_f = (89 - 30 * fuel_frac); % lifetime emissions
    LC = 89; % baseline emissions
    AR = b/c; % aspect ratio
    rho_inf = 0.4135; % density of air at 35,000 ft in kg/m^3
    g = 9.81; % gravity

    range_act = 5837;
    cruise_alt = 35000;
    cruise_mach = .65;
    n_t = 43.10;
    wb = 170000;
    span_e = .9; % span efficiency factor

    n = (span_e * AR) / (span_e * AR + 2);

    alpha = -5:20;
    speed_sound = 295; % m/s
    velocity = cruise_mach * speed_sound;
    viscosity = 1.458e-5;
    Re = rho_inf * velocity * c / viscosity % Reynolds number

    airfoil = append('NACA', int2str(NACA1), int2str(NACA2), sprintf('%02d', NACA34));
    [pol foil] = xfoil(airfoil,alpha,Re,cruise_mach,'ppar n 100','oper iter 200');

    CL = pol.CL * n;
    CD = pol.CD + CL.^2 / (pi * span_e * AR);

    [LpD, i] = max(CL./CD);
    CL = CL(i);
    CD = CD(i);

    d_hf = (42.8 - 20 * fuel_frac^2);
    Sw = b * c;

    speed_sound = 295; % m/s
    L = CL * 0.5 * rho_inf * (cruise_mach * speed_sound)^2 * Sw;

    thickness = NACA34 / 100;
    W_wing = 228 * AR * Sw^(0.5) * (thickness)^(-0.04);

    w1 = wb + W_wing;
    w0 = L;

    W_fuel = 1.06 * (w0 - w1);
    M_fuel = W_fuel/g;

    E_co2 = M_fuel * FCF *(LS_f / LC);
    range = d_hf/g * n_t * LpD * log(w0/w1);
    h = E_co2 + (L - w0)^2;

    range
    Re
    c
    b
end

function [cc,ceq] = constraint(v)
    global cl
    global cd
    global NACA34
    c = v(1);
    b = v(2);
    fuel_frac = v(3);

    FCF = 3.16; % Fuel consumption factor
    LS_f = (89 - 30 * fuel_frac); % lifetime emissions
    LC = 89; % baseline emissions
    AR = b/c; % aspect ratio
    rho_inf = 0.4135; % density of air at 35,000 ft in kg/m^3
    g = 9.81; % gravity

    range_act = 5837;
    cruise_alt = 35000;
    cruise_mach = .65;
    n_t = 43.10;
    wb = 170000 * 4.4482216153;
    span_e = .9; % span efficiency factor

    n = (span_e * AR) / (span_e * AR + 2);

    CL = cl * n;
    CD = cd + CL^2 / (pi * span_e * AR);

    LpD = CL/CD;
    %[LpD, i] = max(CL./CD);

    d_hf = (42.8 - 20 * fuel_frac^2);
    Sw = b * c;

    speed_sound = 295; % m/s
    L = CL * 0.5 * rho_inf * (cruise_mach * speed_sound)^2 * Sw;

    thickness = NACA34 / 100;
    W_wing = 228 * AR * Sw^(0.5) * (thickness)^(-0.04);

    w1 = wb + W_wing;
    w0 = L;

    W_fuel = 1.06 * (w0 - w1);
    M_fuel = W_fuel/g;

    E_co2 = M_fuel * FCF *(LS_f / LC);
    range = d_hf/9.81 * n_t * LpD * log(w0/w1);
    h = E_co2 + (L - w0)^2;

    cc = [];
    ceq = range - range_act;
end