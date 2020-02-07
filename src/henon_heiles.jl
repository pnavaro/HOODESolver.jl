#=
henonHeiles:
- Julia version: 
- Author: ymocquar
- Date: 2019-10-22
=#
"""
 henonHeiles(u, p) : Henon-Heiles differential equation with parameters such that
it can be solved by the Julia solver ODEProblem and solve

Parameters :
    - u : the vector
    - p : the parameter of this version of Henon-Heiles equation, that is epsilon
    - t : the time that is unused.
Ouput : 
    - the derivative vector
"""
function henon_heiles(u,p)
    return [u[3]/p, u[4], -u[1]/p-2u[1]*u[2], -u[2]-u[1]^2+u[2]^2]
end
