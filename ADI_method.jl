using Printf
using LinearAlgebra
using Plots
#using Distributions
#using TimerOutputs
using Profile
using ProgressMeter
using DelimitedFiles
#using NLsolve

# Функция Формирующая трехдиагональную матрицу A
function create_mat(γ, M)
    # нижняя диагональ
    diag_low = fill(-γ, M - 1)
    diag_low[end] = -2γ    # меняем значение последнего элемента вектора

    # главная диагональ
    diag_main = fill(1 + 2γ, M)

    # верхняя диагональ
    diag_up = fill(-γ, M - 1)
    diag_up[1] = -2γ

    return Tridiagonal(diag_low, diag_main, diag_up)
end

function ADI_method(Domain_Length_x, Domain_Length_y, Time_Period, M_x, M_y, N, u_init::Union{Function, Matrix{Float64}}, D, f)

    Δt = Time_Period / (N - 1)
    Δx = (Domain_Length_x) / (M_x - 1)
    Δy = (Domain_Length_y) / (M_y - 1)

    α = (Δt/2 * D) / Δx^2
    β = (Δt/2 * D) / Δy^2

    x = range(-Domain_Length_x/2, Domain_Length_x/2, length = M_x)
    y = range(-Domain_Length_y/2, Domain_Length_y/2, length = M_y)
    t = range(0, Time_Period, length = N)

    # applying initial conditions
    if isa(u_init, Matrix{Float64})
        u = copy(u_init)
    elseif isa(u_init, Function)
        u = [u_init(i, j) for i=x, j=y]
    end
    u_new = copy(u)

    A_x = create_mat(α, M_x)  # матрица A для первого шага
    A_y = create_mat(β, M_y)  # матрица A для второго шага
    B = ones(M_x, M_y)

    for n = 1:N

        # матрица правых частей B
        B[:, 1] = @. (1 - 2β) * u[:, 1] + 2β * u[:, 2] + Δt/2 * [f(t[n], i, 1) for i = 1:M_x]
        B[:, M_y] = @. 2β * u[:, M_y - 1] + (1 - 2β) * u[:, M_y] + Δt/2 * [f(t[n], i, M_y) for i = 1:M_x]
        for j = 2:M_y - 1
            B[:, j] = @. β * u[:, j - 1] + (1 - 2β) * u[:, j] + β * u[:, j + 1] + Δt/2 * [f(t[n], i, j) for i = 1:M_x]
        end

        for j = 1:M_y
            u_new[:, j] = A_x \ B[:, j]
        end

        B[1, :] = @. (1 - 2α) * u_new[1, :] + 2α * u_new[2, :] + Δt/2 * [f(t[n], 1, j) for j = 1:M_y]
        B[M_x, :] = @. 2α * u_new[M_x - 1, :] + (1 - 2α) * u_new[M_x, :] + Δt/2 * [f(t[n], M_x, j) for j = 1:M_y]
        for i = 2:M_x - 1
            B[i, :] = @. α * u_new[i - 1, :] + (1 - 2α) * u_new[i, :] + α * u_new[i + 1, :] + Δt/2 * [f(t[n], i, j) for j = 1:M_y]
        end

        for i = 1:M_x
            u[i, :] = A_y \ B[i, :]
        end

    end

    return u
end


const D = 1

const Domain_Length_x = 5
const Domain_Length_y = 5

const Time_Period = 3.5
const N = 200
const M_x = 101
const M_y = 101

# 1) Область [-a,a]x[-b,b]. Краевые условия: адиабатические стенки (не пропускают тепло).
# При t = 0 Нормальное распределение:
u_init_func(x, y) = exp(- (x^2 + y^2) / 2)

x = range(-Domain_Length_x/2, Domain_Length_x/2, length = M_x)
y = range(-Domain_Length_y/2, Domain_Length_y/2, length = M_y)

u_init = [u_init_func(i, j) for i=x, j=y]

# Источниковый член в уравнении
f(t, x, y) = 0


#Juno.@profiler
u = ADI_method(Domain_Length_x, Domain_Length_y, Time_Period, M_x, M_y, N, u_init_func, D, f)

#open("u_init.txt", "w") do io
#    writedlm(io, u_init)
#end

#=
open("u.txt", "w") do io
    writedlm(io, u)
end
=#


# Writing data into the file.

# https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/
#open("delim_file.txt", "w") do io
#    writedlm(io, [x y])
#end

#matrix = readdlm("delim_file.txt", '\t', Float64, '\n')
