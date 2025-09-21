function f = borehole(rw, r, Tu, Hu, Tl, Hl, L, Kw)
%BOREHOLE  Borehole function (flow rate through a borehole).
%
%   f = BOERHOLE(rw, r, Tu, Hu, Tl, Hl, L, Kw)
%
%   Implements:
%       f = (2*pi*Tu*(Hu - Hl)) / ( ln(r/rw) * ( 1 + (2*L*Tu)/(ln(r/rw)*rw^2*Kw) + Tu/Tl ) )
%
%   Notes
%   -----
%   • Fully vectorized: any input can be a scalar or an array; sizes must be
%     compatible for implicit expansion.
%   • Uses natural logarithm (log).
%
%   Inputs (typical parameter names/units in literature)
%   ----------------------------------------------------
%     rw : radius of borehole (m)
%     r  : radius of influence (m)
%     Tu : transmissivity of upper aquifer (m^2/yr)
%     Hu : potentiometric head of upper aquifer (m)
%     Tl : transmissivity of lower aquifer (m^2/yr)
%     Hl : potentiometric head of lower aquifer (m)
%     L  : length of borehole (m)
%     Kw : hydraulic conductivity of borehole (m/yr)
%
%   Output
%   ------
%     f  : discharge/flow rate (m^3/yr)

    % Basic sanity checks (lightweight)
%    if any(rw <= 0, 'all') || any(r <= 0, 'all')
%        error('rw and r must be positive.');
%    end
%    if any(Kw <= 0, 'all') || any(Tl == 0, 'all')
%        error('Kw must be positive and Tl must be nonzero.');
%    end

    % Common term
    log_term = log(r ./ rw);

    % Denominator inner bracket
    bracket = 1 + (2 .* L .* Tu) ./ (log_term .* rw.^2 .* Kw) + (Tu ./ Tl);

    % Final function
    f = (2 .* pi .* Tu .* (Hu - Hl)) ./ (log_term .* bracket);
end