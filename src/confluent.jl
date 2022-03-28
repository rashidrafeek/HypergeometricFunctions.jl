"""
Compute Kummer's confluent hypergeometric function `M(a, b, z) = ₁F₁(a; b; z)`.
"""
function _₁F₁(a, b, z)
    if real(z) ≥ 0
        return _₁F₁maclaurin(a, b, z)
    else
        # return exp(z)*_₁F₁(b-a, b, -z)
        # See https://github.com/scipy/scipy/blob/3cb6d1901f8a67cd46f96bbd0fe7a3ec6b84851b/scipy/special/specfun/specfun.f#L5352
        a = b - a
        a0 = a
        zp = abs(z)
        hg = 0.0

        nl,la = 0,0
        if a >= 2.0
            nl = 1
            la = Int(a)
            a = a-la-1.0
        end

        y0, y1  = 0.0, 0.0

        for k in 0:nl
            if a0 >= 2.0
                a = a+1.0
            end
            if zp <= (30.0 + abs(b)) || a < 0.0
                hg, rg = 1.0, 1.0
                for j in 1:700
                    rg *= (a+j-1.0) / (j*(b+j-1.0))*zp
                    hg += rg
                    if (hg != 0.0 && abs(rg/hg) <= 1e-15)
                        hg *= exp(z)
                        @goto l1
                    end
                end
            else
                lgb, lga, lgba = loggamma(b), loggamma(a), loggamma(b-a)
                sum1, sum2 = 1.0, 1.0
                r1, r2 = 1.0, 1.0
                for i in 1:16
                    r1 = -r1 * (a+i-1.0) * (a-b+i) / (zp*i)
                    r2 = -r2 * (b-a+i-1.0) * (a-i) / (zp*i)
                    sum1 += r1
                    sum2 += r2
                end

                s1 = exp(lgb-lgba+z) * (zp^(-a)*cos(π*a)) * sum1
                s2 = exp(lgb-lga) * (zp^(a-b)) * sum2

                hg = s1+s2
            end
            @label l1
            if k == 0
                y0 = hg
            else
                y1 = hg
            end
        end
        if a0 >= 2.0
            for i in 1:(la-1)
                hg = ((2.0*a-b+zp)*y1 + (b-a)*y0)/a
                y0 = y1
                y1 = hg
                a = a + 1.0
            end
        end
        return hg
    end
end

"""
Compute Kummer's confluent hypergeometric function `M(a, b, z) = ₁F₁(a; b; z)`.
"""
const M = _₁F₁

"""
Compute Tricomi's confluent hypergeometric function `U(a, b, z) ∼ z⁻ᵃ ₂F₀([a, a-b+1]; []; -z⁻¹)`.
"""
function U(a, b, z)
    return z^-a*pFq([a, a-b+1], typeof(z)[], -inv(z))
end
