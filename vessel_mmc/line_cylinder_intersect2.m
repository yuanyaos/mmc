function H = line_cylinder_intersect2(P0,P1,E0,E1,rt)
% simple interpolation method
% P0, P1: start and end points of the ray
% E0, E1: start and end points of the vessel edge
% rt: radius of the vessel

epsilon = 10^-9;
Lmove = norm(P1-P0);
Ldir = (P1-P0)/Lmove;

% compute unitary vector of vessel
u = E1-E0;
u = u/norm(u);

% project P0 to the bottom of cylinderu
OP = P0-E0;
st = dot(OP,u);
Pt = u*st;
PP0 = OP-Pt;
DIS0 = norm(PP0);

% project P1 to the bottom of cylinder
OP = P1-E0;
st = dot(OP,u);
Pt = u*st;
PP1 = OP-Pt;
DIS1 = norm(PP1);

if (DIS0>rt && DIS1<rt) || (DIS0<rt && DIS1>rt)
    if DIS0>rt && DIS1<rt
        K = (DIS0-rt)/(rt-DIS1);
    else
        K = (rt-DIS0)/(DIS1-rt);
    end
    H = P0+(K/(1+K))*Lmove*Ldir;
else
%     sprintf('No intersection')
    H = [Inf Inf Inf];
end

end