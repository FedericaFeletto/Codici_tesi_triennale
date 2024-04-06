function demo_adaptive_LS(example)

%--------------------------------------------------------------------------
% OBJECT
%--------------------------------------------------------------------------
% This demo shows how to estimate the Lebesgue constant and its relative
% error at degree "deg" of LS-AM pointset "pts" in a certain
% domain, with a fixed tolerance.
%--------------------------------------------------------------------------
% INPUT
%--------------------------------------------------------------------------
%
% example: the variable defines the following domains
%
%     example 1: 'sun';
%     example 2: 'polygon';
%     example 3: 'cardioid';
%     example 4: 'curvpolygon';
%     example 5: 'uniondisks';
%     example 6: 'lune';
%     example 7: 'hypocycloid';
%     example 8: 'epicycloid';
%     example 9: 'epitrochoid';
%     example 10: 'limacon';
%     example 11: 'ellipse';
%     example 12: 'lissajous';
%     example 13: 'deltoid'
%     example 14: 'rhodonea'
%     example 15: 'habenicht_clover';
%     example 16: 'egg'
%     example 17: 'bifolium'
%     example 18: 'torpedo'
%
% deg: Admissible Mesh degree
%
%--------------------------------------------------------------------------
% DATES:
%--------------------------------------------------------------------------
% Written on March 19, 2024 (by F. Feletto and A. Sommariva).
% Checked on March 20, 2024 (by G. Elefante).
%
% Authors involved:
% L. Leokadia Białas-Cież, G. Elefante, F. Feletto, D.J. Kenne,
% A. Sommariva, M. Vianello.
%--------------------------------------------------------------------------
% LICENSE
%--------------------------------------------------------------------------
% Copyright (C) 2023-
% L. Leokadia Białas-Cież, G. Elefante, F. Feletto, D.J. Kenne,
% A. Sommariva, M. Vianello.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%
% Authors:
%
% G. Elefante
% F. Feletto
% L. Leokadia Białas-Cież,
% D.J. Kenne
% Alvise Sommariva <alvise@math.unipd.it>
% Marco Vianello   <marcov@math.unipd.it>
%
% Date: March 19, 2024.
%--------------------------------------------------------------------------

% SETTINGS.
if nargin < 1, example=1; end

warning off;
tol= 10^(-2); % Lebesgue constant relative error tolerance.
mindeg=1;
maxdeg=50;
save_domain_plots=1; % 1. save figs with domain


% -------------------------- Main Code Below ------------------------------

% DEFINE DOMAIN.
domain= define_domain(example);

% COMPUTE "m" SO THAT "tol=c_m-1" IMPLIES "err < tol".
den=2*acos(1/(1+tol));
m=ceil(pi/den);

% EXPERIMENTS RUN.
degV=mindeg:maxdeg;
AM_cost=zeros(degV(end),1);
lebV_LB=zeros(degV(end),1);
lebV_UB=zeros(degV(end),1);
lebV_est=zeros(degV(end),1);



% Actual experiments.
for deg=degV

    % Below:
    % leb_LB, leb_UB: Lebesgue constant bounds
    % lebL_est      : (leb_LB+leb_UB)/2
    % leb_LB  : fine approximation of Lebesgue constant.

    pts=complex_AM(domain,deg,m);
    [leb_LB,zAM,C,domain]=evaluate_leb_const(domain,pts,deg,m);
    leb_UB=C*leb_LB;
    lebL_est=(leb_LB+leb_UB)/2; % estimating Lebesgue constant
    lebV_est(deg,1)=lebL_est;
    lebV_LB(deg,1)=leb_LB;
    lebV_UB(deg,1)=leb_UB;
    AM_cost(deg,1)=C;

    fprintf('\n \t -> deg: %3.0f leb: [%1.3e , %1.3e]  m: %3.0f  AM#: %7.0f',...
        deg,leb_LB,leb_UB,m,length(zAM));

end

% Results.

fprintf('\n \n \t .............. Summary ...............');
fprintf('\n \t Pointset AM');
fprintf('\n \t Example                    : %2.0f',example);
fprintf('\n \t Max Degree                 : %3.0f',degV(end));
fprintf('\n \t Required tolerance (tol)   : %1.5f',tol);
fprintf('\n \t ......................................');
fprintf('\n \n');

% STATISTICS.

fprintf('\n \t ..................LS..................');
fprintf('\n \t Max Degree                  : %3.0f',degV(end));
fprintf('\n \t Leb. const. est. for maxdeg : %3.6f',lebV_est(maxdeg,1));
fprintf('\n \t Leb. const. LB   for maxdeg : %3.6f',lebV_LB(maxdeg,1));
fprintf('\n \t Leb. const. UB   for maxdeg : %3.6f',lebV_UB(maxdeg,1));
fprintf('\n \n');

% PLOT.

% A. Plot pointset
clear_figure(1);
fig=figure(1);
plot(real(zAM),imag(zAM),'.','color',[0.9290 0.6940 0.1250],'MarkerSize',4); hold on;
axis equal;
axis off;
axis tight;
set(gca,'xtick',[])
set(gca,'ytick',[])
hold off;

if save_domain_plots
print(fig,'example_LS.eps','-depsc');
savefig('example_LS.fig');
end

% B. Lebesgue constants plots
clear_figure(2);
fig=figure(2);
semilogy(degV(:),lebV_est(:,1),'color',[0.6350 0.0780 0.1840],'Linewidth',2); hold on;
title('Lebesgue constant estimate');
xlabel('Degree'); ylabel('Lebesgue constants approximation');
grid on;
hold on;

if save_domain_plots
print(fig,'leb_const_est_LS.eps','-depsc');
savefig('leb_const_est_LS.fig');
end



function clear_figure(nfig)

h=figure(nfig);
f_nfig=ishandle(h)&&strcmp(get(h,'type'),'figure');
if f_nfig
    clf(nfig);
end
