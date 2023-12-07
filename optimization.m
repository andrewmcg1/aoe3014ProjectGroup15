
params = [2, 4, 12, 10, .33];
fuel_mass_4series(params)


% input [NACA-1st digit, NACA-2nd digit, NACA-3rd and 4th digits, chord length, fuel fraction]
% solves for a span that satisfies all constraints (not pressure constraint)
% returns the the value of the emitted CO2, if returned value is Inf, then
% the input parameters are not possible
function E_co2 = fuel_mass_4series(v)
    global cl
    global cd
    global LpD
    global rho_inf
    global w0
    global velocity
    global NACA34
    global c
    global fuel_frac
    global L
    global w0
   
    NACA1 = round(round(v(1)));
    NACA2 = round(round(v(2)));
    NACA34 = round(round(v(3)));
    c = v(4);
    fuel_frac = v(5);

    FCF = 3.16; % Fuel consumption factor
    LC = 89; % baseline emissions
    rho_inf = 0.4135; % density of air at 35,000 ft in kg/m^3
    g = 9.81; % gravity

    cruise_alt = 35000;
    cruise_mach = .65;
    n_t = 43.10;
    wb = 170000*4.4482216153;

    alpha = 0:10;
    speed_sound = 295; % m/s
    velocity = cruise_mach * speed_sound;
    viscosity = 1.458e-5;
    Re = rho_inf * velocity * c / viscosity % Reynolds number
    %Re = 1e6;

    airfoil = append('NACA', int2str(NACA1), int2str(NACA2), sprintf('%02d', NACA34));
    [pol foil] = xfoil(airfoil,alpha,Re,cruise_mach,'ppar T 0.2', 'ppar n 80','oper iter 10000');

    cl = pol.CL;
    cd = pol.CD;

    b = fmincon(@solveSpan, 100, [], [], [], [], 50, 250)
    if L >= w0*(1-1e-5)
        LS_f = (89 - 30 * fuel_frac); % lifetime emissions
        Sw = b * c;
        AR = b/c; % aspect ratio
        thickness = NACA34 / 100;
        W_wing = (228 * AR * (Sw*10.764)^(0.5) * (thickness)^(-0.04)) * 4.4482216153;

        w1 = wb + W_wing;
        W_fuel = 1.06*(w0 - w1);
        M_fuel = W_fuel / g;

        d_hf = (42.8 - 20 * fuel_frac^2);
        range = d_hf/9.81 * n_t * LpD * log(w0/w1);
        E_co2 = M_fuel * FCF *(LS_f / LC);
    else
        E_co2 = Inf;
    end

end

function out = solveSpan(v)
    global cl
    global cd
    global LpD
    global LpD
    global rho_inf
    global w0
    global velocity
    global NACA34
    global c
    global fuel_frac
    global L
    global w0
    b = v;

    Sw = b * c;
    AR = b/c; % aspect ratio

    span_e = .9; % span efficiency factor
    n = (span_e * AR) / (span_e * AR + 2);
    n_t = 43.10;
    wb = 170000*4.4482216153;
    range_act = 5837;%*1000;
    g = 9.81; % gravity

    cruise_mach = .65;
    speed_sound = 295; % m/s
    velocity = cruise_mach * speed_sound;

    CL = cl * n;
    CD = cd + CL.^2 / (pi * span_e * AR);

    [LpD, i] = max(CL./CD);
    CL = CL(i);
    CD = CD(i);
    %[LpD, i] = max(pol.CL./pol.CD)

    thickness = NACA34 / 100;
    W_wing = (228 * AR * (Sw*10.764)^(0.5) * (thickness)^(-0.04)) * 4.4482216153;

    w1 = wb + W_wing;
    d_hf = (42.8 - 20 * fuel_frac^2);
    w0 = exp(g*range_act / (d_hf*n_t * LpD)) * w1;

    L = CL * 0.5 * rho_inf * (velocity)^2 * Sw;

    out = abs(L - w0);

end