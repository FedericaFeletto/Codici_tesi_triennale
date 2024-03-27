
function demo_adaptive(example)

%--------------------------------------------------------------------------
% OBJECT:
%--------------------------------------------------------------------------
% This demo shows how to estimate the Lebesgue constant and its relative
% error at degree "deg" of the interpolation at pointset "pts" in a certain
% domain, with a fixed tolerance.
%--------------------------------------------------------------------------
% INPUT:
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
if nargin < 1, example=14; end

warning off;
tol= 10^(-2); % Lebesgue constant relative error tolerance.
m_pts=4; % Constant of the AM set used for intp. points extraction.
maxdeg=20;
save_domain_plots=1; % 1. save figs with domain


% -------------------------- Main Code Below ------------------------------

% DEFINE DOMAIN.
domain= define_domain(example);

% COMPUTE "m" SO THAT "tol=c_m-1" IMPLIES "err < tol".
den=2*acos(1/(1+tol));
m=ceil(pi/den);

% EXPERIMENTS RUN.
degV=1:maxdeg;
lejaV=0:2;
AM_cost=zeros(degV(end),3);
lebV_est=zeros(degV(end),3);
lebV_LB=zeros(degV(end),3);
lebV_UB=zeros(degV(end),3);
pts_leja={}; 

for leja=lejaV

    % Displays.
    fprintf('\n \t * Running tests on: ');
    switch leja
        case 0
            set_leja_L='AFP';
        case 1
            set_leja_L='DLP';
        case 2
            set_leja_L='PLP';
    end
    disp(set_leja_L);
    zAM_forall_deg={};
    pts_forall_deg={};

    % Actual experiments.
    for deg=degV
        pts=interp_pointset(domain,leja,deg,m_pts);
        [leb_constant,zAM,C,domain]=evaluate_leb_const(domain,pts,deg,m);
        leb_LB=leb_constant;
        leb_UB=C*leb_LB;   
        lebL_est=(leb_LB+leb_UB)/2; % estimating Lebesgue constant for AFP/DLP/PLP
        lebV_est(deg,leja+1)=lebL_est;
        lebV_LB(deg,leja+1)=leb_LB;
        lebV_UB(deg,leja+1)=leb_UB;
        AM_cost(deg,leja+1)=C;
        zAM_forall_deg{end+1}=zAM;
        pts_forall_deg{end+1}=pts;
        

        fprintf('\n \t -> deg: %2.0f leb: %1.3e  m_pts: %2.0f m: %2.0f  AM#: %6.0f',...
            deg,leb_constant,m_pts,m,length(zAM));

    end

    zAM_forall_deg{end+1}=zAM_forall_deg; 
    pts_leja{end+1}=pts_forall_deg;

    % Results
    fprintf('\n \n \t ........................ Summary ........................');

    fprintf('\n \t Example                    : %2.0f',example);
    fprintf('\n \t Pointset                   : %s',set_leja_L);
    fprintf('\n \t Max Degree                 : %2.0f',degV(end));
    fprintf('\n \t Required tolerance (tol)   : %3.6f',tol);
    fprintf('\n \t AM factor (m_pts)          : %2.0f',m_pts);
    fprintf('\n \t AM card for maxdeg         : %5.0f',length(zAM_forall_deg{1,maxdeg}));
    fprintf('\n \t .........................................................');
    fprintf('\n \n');

end


% STATISTICS.
fprintf('\n \t ...................AFP..................');
fprintf('\n \t Max Degree                  : %3.0f',degV(end));
fprintf('\n \t Leb. const. est. for maxdeg : %3.6f',lebV_est(maxdeg,1));
fprintf('\n \t Leb. const. LB   for maxdeg : %3.6f',lebV_LB(maxdeg,1));
fprintf('\n \t Leb. const. UB   for maxdeg : %3.6f',lebV_UB(maxdeg,1));
fprintf('\n \n');

fprintf('\n \t ...................DLP..................');
fprintf('\n \t Max Degree                  : %3.0f',degV(end));
fprintf('\n \t Leb. const. est. for maxdeg : %3.6f',lebV_est(maxdeg,2));
fprintf('\n \t Leb. const. LB   for maxdeg : %3.6f',lebV_LB(maxdeg,2));
fprintf('\n \t Leb. const. UB   for maxdeg : %3.6f',lebV_UB(maxdeg,2));
fprintf('\n \n');

fprintf('\n \t ...................PLP..................');
fprintf('\n \t Max Degree                  : %3.0f',degV(end));
fprintf('\n \t Leb. const. est. for maxdeg : %3.6f',lebV_est(maxdeg,3));
fprintf('\n \t Leb. const. LB   for maxdeg : %3.6f',lebV_LB(maxdeg,3));
fprintf('\n \t Leb. const. UB   for maxdeg : %3.6f',lebV_UB(maxdeg,3));
fprintf('\n \n');


% Type of pointset in figure 1:
Leja=2;  % 1: AFP, 2: DLP, 3: PLP
set_Leja=pts_leja{Leja};
pts=set_Leja{end};


% PLOT.

% A. Lebesgue constants plots
clear_figure(1);
fig=figure(1);
plot(real(zAM),imag(zAM),'.','color',[0.9290 0.6940 0.1250],'MarkerSize',4); hold on;
plot(real(pts),imag(pts),'o','MarkerEdgeColor','k','MarkerFaceColor','c','MarkerSize',6);
axis equal;
axis off;
axis tight;
set(gca,'xtick',[])
set(gca,'ytick',[])
hold off;

if save_domain_plots
print(fig,'example.eps','-depsc');
savefig('example.fig');
end


% B. Plot domain and pointsets
clear_figure(2);
fig=figure(2);
plot(degV(:),lebV_est(:,1),'m','Linewidth',2); hold on;
plot(degV(:),lebV_est(:,2),'color',[0.4940 0.1840 0.5560],'Linewidth',2); hold on;
plot(degV(:),lebV_est(:,3),'g','Linewidth',2); hold on;
title('Lebesgue constant estimate');
xlabel('Degree'); ylabel('Lebesgue constants approximation');
legend('AFP','DLP','PLP','Location','northwest');
grid on;
hold on;

if save_domain_plots
print(fig,'leb_const_est.eps','-depsc');
savefig('leb_const_est.fig');
end


function clear_figure(nfig)

h=figure(nfig);
f_nfig=ishandle(h)&&strcmp(get(h,'type'),'figure');
if f_nfig
    clf(nfig);
end

