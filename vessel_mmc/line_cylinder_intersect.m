function H = line_cylinder_intersect(P0,P1,E0,E1,rt)
% projection method
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
    theta = dot(PP0,PP1);
    theta = theta/(DIS0*DIS1+epsilon);
    
    pt0(1) = DIS0;
    pt0(2) = 0;
    pt1(1) = DIS1*theta;
    pt1(2) = DIS1*sqrt(1-theta*theta);
    dx = pt1(1)-pt0(1);
    dy = pt1(2)-pt0(2);
    dr2 = dx*dx+dy*dy;
    drr2 = 1/dr2;
    Dt = pt0(1)*pt1(2);
    delta = sqrt(rt*rt*dr2-Dt*Dt);
    if dy>=0
        sgn = 1;
    else
        sgn = -1;
    end
    xa = Dt*dy*drr2;
    xb = sgn*dx*delta*drr2;
    ya = -Dt*dx*drr2;
    yb = abs(dy)*delta*drr2;
    ph0(1) = xa+xb;
    ph0(2) = ya+yb;
    ph1(1) = xa-xb;
    ph1(2) = ya-yb;
    Dp = (pt1(1)-pt0(1))*(pt1(1)-pt0(1))+(pt1(2)-pt0(2))*(pt1(2)-pt0(2));
    if DIS0>rt && DIS1<rt
        Lratio = (ph1(1)-pt0(1))*(ph1(1)-pt0(1))+(ph1(2)-pt0(2))*(ph1(2)-pt0(2));
        if(Dp>Lratio)
            Lratio = sqrt(Lratio/Dp);
        else
            Lratio = (ph0(1)-pt0(1))*(ph0(1)-pt0(1))+(ph0(2)-pt0(2))*(ph0(2)-pt0(2));
            Lratio = sqrt(Lratio/Dp);
        end
    else
        Lratio = (ph1(1)-pt1(1))*(ph1(1)-pt1(1))+(ph1(2)-pt1(2))*(ph1(2)-pt1(2));
        if(Dp>Lratio)
            Lratio = 1-sqrt(Lratio/Dp);
        else
            Lratio = (ph0(1)-pt1(1))*(ph0(1)-pt1(1))+(ph0(2)-pt1(2))*(ph0(2)-pt1(2));
            Lratio = 1-sqrt(Lratio/Dp);
        end
    end
    H = P0+(Lratio*Lmove)*Ldir;
else
%     sprintf('No intersection')
    H = [Inf Inf Inf];
end

end