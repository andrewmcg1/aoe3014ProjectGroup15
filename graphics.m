airfoil = 'NACA2412';
alpha = -15:0.1:15;
Re = 1e6;
cruise_mach = 0;


[pol foil] = xfoil(airfoil,alpha,Re,cruise_mach,'ppar T 0.2', 'ppar n 160','oper iter 200','oper vpar xtr 0.1 1')

figure
hold on
plot(pol.alpha,pol.CL)
plot(pol.alpha,pol.CD)
hold off

figure
plot(pol.alpha,pol.Cm)
title('Cm vs Angle of Attack');
xlabel('Angle of Attack (deg)');
ylabel('Cm');
saveas(gcf,'cm_vs_alpha.png');


alpha = 5;
i = find(foil.alpha==alpha);
[foil.x(:,i), foil.Dstar(:,i)]
figure
plot(foil.x(:,i),foil.Dstar(:,i))


uiwait(gcf);