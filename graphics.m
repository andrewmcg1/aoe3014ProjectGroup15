airfoil = 'NACA2412';
alpha = -15:0.1:15;
Re = 1e6;
cruise_mach = 0;


[pol foil] = xfoil(airfoil,alpha,Re,cruise_mach,'ppar T 0.2', 'ppar n 160','oper iter 200','oper vpar xtr 0.1 1');

figure
hold on
plot(pol.alpha,pol.CL)
plot(pol.alpha,pol.CD)

[min_val, i] = min(abs(pol.CL));
cl_0 = pol.alpha(i);
[min_val, i] = min(abs(pol.CD));
cd_min = pol.alpha(i);
str = {sprintf('{\\alpha}_{L=0} = %.1f ^{\\circ}',cl_0),sprintf('{\\alpha}_{D, min} = %.1f ^{\\circ}',cd_min)};
text(-12,1.5,str)

title(append('{c}_{l} vs \alpha and {c}_{d} vs \alpha (', airfoil, ')'));
legend('{c}_{l}','{c}_{d}')
xlabel('^{x}/_{c}')
saveas(gcf,'cdcl_vs_alpha.png');
hold off

figure
plot(pol.alpha,pol.Cm)
title('{c}_{m,c/4} vs \alpha');
title(append('{c}_{m,c/4} vs \alpha (', airfoil, ')'));
xlabel('\alpha (deg)');
ylabel('{c}_{m,c/4}');
saveas(gcf,'cm_vs_alpha.png');


alpha = 5;
i = find(foil.alpha==alpha);
figure
plot(foil.x(:,i),foil.Dstar(:,i))
title(append('Displacement Thickness vs ^{x}/_{c} (', airfoil, ', \alpha = ', num2str(alpha), '^{\circ})'));
xlabel('^{x}/_{c}');
ylabel('\delta^{*}');
saveas(gcf,'boundary_vs_x.png');


uiwait(gcf);