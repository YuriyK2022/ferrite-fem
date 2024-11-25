# Завантаження необхідних бібліотек та залежностей
using Ferrite, Tensors, SparseArrays, LinearAlgebra, Printf

# Визначення J₂-пластичністі матеріалу, що містить параметри матеріалу та пружну жорсткість Dᵉ
struct J2Plasticity{T, S <: SymmetricTensor{4, 3, T}}
    G::T  # Модуль зсуву (Shear modulus)
    K::T  # Об'ємний модуль (Bulk modulus)
    σ₀::T # Начальне допустиме напруження текучості (Yield Strength)
    H::T  # Модуль зміцнення (Hardening modulus)
    Dᵉ::S # Тензор пружної жорсткості (Elastic stiffness tensor)
end;

# Визначення конструктора для екземпляра матеріалу
function J2Plasticity(E, ν, σ₀, H)
    δ(i,j) = i == j ? 1.0 : 0.0 # helper function
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)

    Isymdev(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    Dᵉ = SymmetricTensor{4, 3}(temp)
    return J2Plasticity(G, K, σ₀, H, Dᵉ)
end;

# Визначення конструктора для збереження стану матеріала у точках Гауса
struct MaterialState{T, S <: SecondOrderTensor{3, T}}
    # Store "converged" values
    ϵᵖ::S               # пластичні деформації (plastic strain)
    σ::S                # напруження (stress)
    k::T                # змінна зміцнення (hardening variable)
end

# Конструктор для ініціалізації начального стану матеріала. Кожна величина встановлена ​​на нуль.
function MaterialState()
    return MaterialState(
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                0.0)
end

# Для подальшого використання на етапі постобробки визначимо функцію для обчислення ефективного напруження Мізеса.
function vonMises(σ)
    s = dev(σ)
    return sqrt(3.0/2.0 * s ⊡ s)
end;

# Метод, який обчислює напругу та дотичну жорсткість матеріалу в даній точці інтегрування.
# Вхідними даними є поточна деформація та стан матеріалу за попередній крок за часом.
function compute_stress_tangent(ϵ::SymmetricTensor{2, 3}, material::J2Plasticity, state::MaterialState)
    # зчитування параметров матеріалу
    G = material.G
    H = material.H

    
    σᵗ = material.Dᵉ ⊡ (ϵ - state.ϵᵖ) 
    sᵗ = dev(σᵗ)                       # девіаторна частина тензору напружень
    J₂ = 0.5 * sᵗ ⊡ sᵗ                # другий інваріант для sᵗ
    σᵗₑ = sqrt(3.0*J₂)                 # ефективні напруження (von Mises stress)
    σʸ = material.σ₀ + H * state.k     # попередній yield limit

    φᵗ  = σᵗₑ - σʸ                     # значення для поверхні текучості

    if φᵗ < 0.0                        # пружне навантаження
        return σᵗ, material.Dᵉ, MaterialState(state.ϵᵖ, σᵗ, state.k)
    else                               # пластичне навантаження
        h = H + 3G
        μ =  φᵗ / h                    # plastic multiplier

        c1 = 1 - 3G * μ / σᵗₑ
        s = c1 * sᵗ                    # оновлення девіаторних напружень
        σ = s + vol(σᵗ)                # оновлення тензору напружень

        # Обчислення тангенційної жорсткості
        κ = H * (state.k + μ)          # drag stress
        σₑ = material.σ₀ + κ           # оновлення поверхні текучості

        δ(i,j) = i == j ? 1.0 : 0.0
        Isymdev(i,j,k,l)  = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
        Q(i,j,k,l) = Isymdev(i,j,k,l) - 3.0 / (2.0*σₑ^2) * s[i,j]*s[k,l]
        b = (3G*μ/σₑ) / (1.0 + 3G*μ/σₑ)

        Dtemp(i,j,k,l) = -2G*b * Q(i,j,k,l) - 9G^2 / (h*σₑ^2) * s[i,j]*s[k,l]
        D = material.Dᵉ + SymmetricTensor{4, 3}(Dtemp)

        # Обчислення нового стану матеріала
        Δϵᵖ = 3/2 * μ / σₑ * s        # plastic strain
        ϵᵖ = state.ϵᵖ + Δϵᵖ           # plastic strain
        k = state.k + μ               # hardening variable
        return σ, D, MaterialState(ϵᵖ, σ, k)
    end
end

# Далі наведено функції для складання та розв’язання FE-задачі
function create_values(interpolation)
    
    qr      = QuadratureRule{RefTetrahedron}(2)           # задання правила квадратур
    facet_qr = FacetQuadratureRule{RefTetrahedron}(3)

    # формування ячійок та фасетів
    cellvalues_u = CellValues(qr, interpolation)
    facetvalues_u = FacetValues(facet_qr, interpolation)

    return cellvalues_u, facetvalues_u
end;

# Формування ступенів свободи (Degrees for Dreedom)
function create_dofhandler(grid, interpolation)
    dh = DofHandler(grid)
    add!(dh, :u, interpolation) # додавання полів переміщень для троьх компонентів
    close!(dh)
    return dh
end

# Формування граничних умов (Boundary Condition)
function create_bc(dh, grid)
    dbcs = ConstraintHandler(dh)
    # Вводимо закріплення балки зліва
    dofs = [1, 2, 3]
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x,t) -> [0.0, 0.0, 0.0], dofs)
    add!(dbcs, dbc)
    close!(dbcs)
    return dbcs
end;

# Збірка матриць системи
function doassemble!(K::SparseMatrixCSC, r::Vector, cellvalues::CellValues, dh::DofHandler,
                     material::J2Plasticity, u, states, states_old)
    assembler = start_assemble(K, r)
    nu = getnbasefunctions(cellvalues)
    re = zeros(nu)     # вектор нев'язки елемента
    ke = zeros(nu, nu) # тангенціна матриця елемента

    for (i, cell) in enumerate(CellIterator(dh))
        fill!(ke, 0)
        fill!(re, 0)
        eldofs = celldofs(cell)
        ue = u[eldofs]
        state = @view states[:, i]
        state_old = @view states_old[:, i]
        assemble_cell!(ke, re, cell, cellvalues, material, ue, state, state_old)
        assemble!(assembler, eldofs, ke, re)
    end
    return K, r
end

 
function assemble_cell!(Ke, re, cell, cellvalues, material,
                        ue, state, state_old)
    n_basefuncs = getnbasefunctions(cellvalues)
    reinit!(cellvalues, cell)

    for q_point in 1:getnquadpoints(cellvalues)
        # For each integration point, compute stress and material stiffness
        ϵ = function_symmetric_gradient(cellvalues, q_point, ue) # Total strain
        σ, D, state[q_point] = compute_stress_tangent(ϵ, material, state_old[q_point])

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δϵ = shape_symmetric_gradient(cellvalues, q_point, i)
            re[i] += (δϵ ⊡ σ) * dΩ # add internal force to residual
            for j in 1:i # loop only over lower half
                Δϵ = shape_symmetric_gradient(cellvalues, q_point, j)
                Ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
            end
        end
    end
    symmetrize_lower!(Ke)
end

# Допоміжна функція для симетризації тангенційної матриці матеріалу
function symmetrize_lower!(K)
    for i in 1:size(K,1)
        for j in i+1:size(K,1)
            K[i,j] = K[j,i]
        end
    end
end;

function doassemble_neumann!(r, dh, facetset, facetvalues, t)
    n_basefuncs = getnbasefunctions(facetvalues)
    re = zeros(n_basefuncs)                      # element residual vector
    for fc in FacetIterator(dh, facetset)
        # Add traction as a negative contribution to the element residual `re`:
        reinit!(facetvalues, fc)
        fill!(re, 0)
        for q_point in 1:getnquadpoints(facetvalues)
            dΓ = getdetJdV(facetvalues, q_point)
            for i in 1:n_basefuncs
                δu = shape_value(facetvalues, q_point, i)
                re[i] -= (δu ⋅ t) * dΓ
            end
        end
        assemble!(r, celldofs(fc), re)
    end
    return r
end

# СТВОРЕННЯ ОСНОВНОЇ ФУНКЦІЇ ДЛЯ ОБЧИСЛЕНЬ
function solve()
    # Задання параметрів матеріалу
    E = 200.0e9         # [Pa]
    H = E/20            # [Pa]
    ν = 0.3             # [-]
    σ₀ = 200e6          # [Pa]
    material = J2Plasticity(E, ν, σ₀, H)

    L = 10.0 # довжина балки length [m]
    w = 1.0  # товщина балки width [m]
    h = 1.0  # висота балки  height[m]
    n_timesteps = 10
    u_max = zeros(n_timesteps)
    traction_magnitude = 1.e7 * range(0.5, 1.0, length=n_timesteps)

    # Створення геометрії балки за розмірами, накладення dofs та граничних умов (boundary conditions)
    n = 7
    nels = (10n, n, 2n) # кількість елементів у кожному просторовому напрямку
    P1 = Vec((0.0, 0.0, 0.0))  # start point for geometry
    P2 = Vec((L, w, h))        # end point for geometry
    grid = generate_grid(Tetrahedron, nels, P1, P2)
    interpolation = Lagrange{RefTetrahedron, 1}()^3

    dh = create_dofhandler(grid, interpolation)
    dbcs = create_bc(dh, grid) # створення граничних умов Діріхле

    cellvalues, facetvalues = create_values(interpolation)

    # Попереднє розподілення векторів рішення
    n_dofs = ndofs(dh)  # сумарне число ступенів свободи dofs
    u  = zeros(n_dofs)  # вектор рішення для переміщень
    Δu = zeros(n_dofs)  # приріст переміщень
    r = zeros(n_dofs)   # вектор нев'язок
    K = allocate_matrix(dh); # тангенційна матриця жорсткості

    # Створення стану матеріалу. Створюється окремий масив для кожної ячійки
    # Кожен елемент є масивом матеріальних станів - по одному для кожної точки інтегрування (integration points)
    nqp = getnquadpoints(cellvalues)
    states = [MaterialState() for _ in 1:nqp, _ in 1:getncells(grid)]
    states_old = [MaterialState() for _ in 1:nqp, _ in 1:getncells(grid)]

    # Цикл Ньютона-Рафсона
    NEWTON_TOL = 1
    print("\n Starting Netwon-Raphson Iterations:\n")

    for timestep in 1:n_timesteps
        t = timestep
        traction = Vec((0.0, 0.0, traction_magnitude[timestep]))
        newton_itr = -1
        print("\n Time Step @time = $timestep:\n")
        update!(dbcs, t)        
        apply!(u, dbcs)         

        while true; newton_itr += 1

            if newton_itr > 8
                error("Досягнуто максимум ітерацій Ньютона-Рафсона, рішення зупинено")
                break
            end
            
            doassemble!(K, r, cellvalues, dh, material, u, states, states_old);
            
            doassemble_neumann!(r, dh, getfacetset(grid, "right"), facetvalues, traction)
            norm_r = norm(r[Ferrite.free_dofs(dbcs)])

            print("Iteration: $newton_itr \tresidual: $(@sprintf("%.8f", norm_r))\n")
            if norm_r < NEWTON_TOL
                break
            end

            apply_zero!(K, r, dbcs)
            Δu = Symmetric(K) \ r
            u -= Δu
        end

        # Оновлення попереднього стану матеріалу при успішній ітерації на наступному кроку за часом
        states_old .= states

        u_max[timestep] = maximum(abs, u) # максимальні переміщення на поточному кроці за часом
    end

    # Отримання результатів розрахунків
    # Oстворюємо vtu-файл з фінальнми результатами розрахунків
    mises_values = zeros(getncells(grid))
    κ_values = zeros(getncells(grid))
    for (el, cell_states) in enumerate(eachcol(states))
        for state in cell_states
            mises_values[el] += vonMises(state.σ)
            κ_values[el] += state.k*material.H
        end
        mises_values[el] /= length(cell_states) # середні значення напружень Мізеса (von Mises Stress)
        κ_values[el] /= length(cell_states)     # середні напруження зсуву (Shear Stress)
    end
    VTKGridFile("plasticity", dh) do vtk
        write_solution(vtk, dh, u)              # поле пеерміщень
        write_cell_data(vtk, mises_values, "von Mises [Pa]")
        write_cell_data(vtk, κ_values, "Drag stress [Pa]")
    end

    return u_max, traction_magnitude
end

u_max, traction_magnitude = solve();

# Візуалізація результатів розрахунків за допомогою біліотекі Plots
using Plots
plot(
    vcat(0.0, u_max),                            # додавання центру координат у вигляді крапки
    vcat(0.0, traction_magnitude),
    linewidth=2,
    title="Traction-Displacement Curve",
    label=nothing,
    markershape=:auto
    )
ylabel!("Traction Force [Pa]")
xlabel!("Maximum deflection [meters]")