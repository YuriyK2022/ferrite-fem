# Моделювання гіперпружного деформування матеріала засобами бібліотеки Ferrite.jl

# Завантаження необхідних біліотек та залежностей
using Ferrite, Tensors, TimerOutputs, ProgressMeter, IterativeSolvers

struct NeoHooke
    μ::Float64
    λ::Float64
end
# Створення функції для опису моделі гіперпружності Нео-Хукина
function Ψ(C, mp::NeoHooke)
    μ = mp.μ
    λ = mp.λ
    Ic = tr(C)
    J = sqrt(det(C))
    return μ / 2 * (Ic - 3 - 2 * log(J)) + λ / 2 * (J - 1)^2
end

# Функція для обчислення всіх похідних
function constitutive_driver(C, mp::NeoHooke)
    # Compute all derivatives in one function call
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, mp), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end;

# Функція збірки матриць елементів
function assemble_element!(ke, ge, cell, cv, fv, mp, ue, ΓN)
    # Повторне ініціалізування значеннь в ячійках сітки
    reinit!(cv, cell)
    fill!(ke, 0.0)
    fill!(ge, 0.0)

    b = Vec{3}((0.0, -0.5, 0.0)) # Об'ємне зусилля
    tn = 0.1 
    ndofs = getnbasefunctions(cv)

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        # Обчислення градієнту деформації F і тензора Коші-Гріна C
        ∇u = function_gradient(cv, qp, ue)
        F = one(∇u) + ∇u
        C = tdot(F) # F' ⋅ F
    
        S, ∂S∂C = constitutive_driver(C, mp)
        P = F ⋅ S
        I = one(S)
        ∂P∂F =  otimesu(I, S) + 2 * otimesu(F, I) ⊡ ∂S∂C ⊡ otimesu(F', I)

        
        for i in 1:ndofs
            
            δui = shape_value(cv, qp, i)
            ∇δui = shape_gradient(cv, qp, i)
            
            ge[i] += ( ∇δui ⊡ P - δui ⋅ b ) * dΩ

            ∇δui∂P∂F = ∇δui ⊡ ∂P∂F 
            for j in 1:ndofs
                ∇δuj = shape_gradient(cv, qp, j)
                
                ke[i, j] += ( ∇δui∂P∂F ⊡ ∇δuj ) * dΩ
            end
        end
    end

    # Обчислення фнтегралу п поверхні для зовнішніх зусилль
    for facet in 1:nfacets(cell)
        if (cellid(cell), facet) in ΓN
            reinit!(fv, cell, facet)
            for q_point in 1:getnquadpoints(fv)
                t = tn * getnormal(fv, q_point)
                dΓ = getdetJdV(fv, q_point)
                for i in 1:ndofs
                    δui = shape_value(fv, q_point, i)
                    ge[i] -= (δui ⋅ t) * dΓ
                end
            end
        end
    end
end;

function assemble_global!(K, g, dh, cv, fv, mp, u, ΓN)
    n = ndofs_per_cell(dh)
    ke = zeros(n, n)
    ge = zeros(n)

    # Збірка матриць K и g
    assembler = start_assemble(K, g)

    # Цикл по всім ячійкам сітки
    @timeit "assemble" for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] # елементні ступені свободи
        @timeit "element assemble" assemble_element!(ke, ge, cell, cv, fv, mp, ue, ΓN)
        assemble!(assembler, global_dofs, ke, ge)
    end
end;

# Функція для рішення задачі
function solve()
    reset_timer!()

    # Генерація розрахункової сітки
    N = 10
    L = 1.0
    left = zero(Vec{3})
    right = L * ones(Vec{3})
    grid = generate_grid(Tetrahedron, (N, N, N), left, right)

    # Влсастивості матеріалу
    E = 10.0
    ν = 0.3
    μ = E / (2(1 + ν))
    λ = (E * ν) / ((1 + ν) * (1 - 2ν))
    mp = NeoHooke(μ, λ)

    # Параметри кінцевих елементів
    ip = Lagrange{RefTetrahedron, 1}()^3
    qr = QuadratureRule{RefTetrahedron}(1)
    qr_facet = FacetQuadratureRule{RefTetrahedron}(1)
    cv = CellValues(qr, ip)
    fv = FacetValues(qr_facet, ip)

    # DofHandler
    dh = DofHandler(grid)
    add!(dh, :u, ip) # Додавання переміщень
    close!(dh)

    # Функція для моделювання прикладення скручуючого моменту к моделі
    function rotation(X, t)
        θ = pi / 3 # 60°
        x, y, z = X
        return t * Vec{3}((
            0.0,
            L/2 - y + (y-L/2)*cos(θ) - (z-L/2)*sin(θ),
            L/2 - z + (y-L/2)*sin(θ) + (z-L/2)*cos(θ)
        ))
    end

    dbcs = ConstraintHandler(dh)
    # Додавання граничних умов для закріплення моделі
    dbc = Dirichlet(:u, getfacetset(grid, "right"), (x,t) -> [0.0, 0.0, 0.0], [1, 2, 3])
    add!(dbcs, dbc)
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x,t) -> rotation(x, t), [1, 2, 3])
    add!(dbcs, dbc)
    close!(dbcs)
    t = 0.5
    Ferrite.update!(dbcs, t)

    # Граничні умови Неймана
    ΓN = union(
        getfacetset(grid, "top"),
        getfacetset(grid, "bottom"),
        getfacetset(grid, "front"),
        getfacetset(grid, "back"),
    )

    # Попереднє виділення векторів для рішення та приростів Ньютона
    _ndofs = ndofs(dh)
    un = zeros(_ndofs) # попередній вектор рішення
    u  = zeros(_ndofs)
    Δu = zeros(_ndofs)
    ΔΔu = zeros(_ndofs)
    apply!(un, dbcs)

    # Створення розрідженої матриці та вектору нев'язки
    K = allocate_matrix(dh)
    g = zeros(_ndofs)

    # Обчислення ітерацій за Ньютоном
    newton_itr = -1
    NEWTON_TOL = 1e-8
    NEWTON_MAXITER = 30
    prog = ProgressMeter.ProgressThresh(NEWTON_TOL; desc = "Solving:")

    # Цикл іткрацій за методом Ньютона-Рафсона
    while true; newton_itr += 1
        # Створення кроку прирісту за деформаціями
        u .= un .+ Δu
        # Обчислення нев'язок та тангенційної матриці жорсткості на поточном кроку деформацій
        assemble_global!(K, g, dh, cv, fv, mp, u, ΓN)
        # Прикладення граничніих умов
        apply_zero!(K, g, dbcs)
        # Обчислення норми нев'язок та порівняння с заданим допуском (tolerance)
        normg = norm(g)
        ProgressMeter.update!(prog, normg; showvalues = [(:iter, newton_itr)])
        if normg < NEWTON_TOL
            break
        elseif newton_itr > NEWTON_MAXITER
            error("Досягнуто максимальне число ітерацій Ньютона-Рафсона! Рішення зупинено.")
        end

        # Обчислення ітерацій за допомогою Conjugate Gradients Solver
        @timeit "linear solve" IterativeSolvers.cg!(ΔΔu, K, g; maxiter=1000)

        apply_zero!(ΔΔu, dbcs)
        Δu .-= ΔΔu
    end

    # Збереженя результатів обчислень у файлі VTK для подальшої візуалізації
    @timeit "export" begin
        VTKGridFile("hyperelasticity", dh) do vtk
            write_solution(vtk, dh, u)
        end
    end

    print_timer(title = "Analysis with $(getncells(grid)) elements", linechars = :ascii)
    return u
end

u = solve();
