"""
    ExtraSpecialFunctions

So far only the [Chapman function](https://en.wikipedia.org/wiki/Chapman_function)

"""
module ExtraSpecialFunctions

export Chapman

using SpecialFunctions

"""
    Chapman(Z, χ)

asymptotic expansion, see  Huestis, David L. (2001), "Accurate evaluation of the Chapman function for atmospheric attenuation",
J. Quantitative Spectroscopy & Radiative Transfer, 69, 709–721, doi:10.1016/S0022-4073(00)00107-2.

Z    dimensionless height = (R + z)/H
R    radius of the planet/moon
H    scale iheight (unit is typically km)
"""
function Chapman(Z, χ, method::Val{Symbol}=Val(:Huestis))
    sinχ = sin(χ)
    zf = Z(1 - sinχ)

    √(πZ/(1 + sinχ))*exp(zf)*erfc(√(zf))
end

end # module
