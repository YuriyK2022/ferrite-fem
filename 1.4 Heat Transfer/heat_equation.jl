# Розрахунок температурногого поля у 2D-постановці
# з внутрішнім рівномірним істочником тепла,
# яке було реалізовано через однорідні граничні умови Діріхле.

using Ferrite, SparseArrays

# Була створена розрахункова сітка на основі чотирьохкутних елементів
# розміром 20x20 із застосуванням функції generate_grid.
# Генератор за замовчуванням використовує квадратну геометричну область,
# тому не треба вказувати кути розрахункової області.

println("Enter mesh dimension dx: ")
input_dx = readline()
dx = parse(Int, input_dx)

println("Enter mesh dimension dy: ")
input_dy = readline()
dy = parse(Int, input_dy)

grid = generate_grid(Quadrilateral, (dx, dy))

# CellValues ​​спрощує процес оцінки значень та градієнтів тестових і пробних функцій.
# Нам потрібно вказати інтерполяційний простір для функцій форми.
# Ми використовуємо функції Лагранжа, засновані на двовимірному чотирикутнику
# Ми також визначаємо квадратурне правило на основі того ж
# опорного елемента. Ми об'єднуємо інтерполяцію і правило квадратури в об'єкт CellValues.
ip = Lagrange{RefQuadrilateral, 1}()
qr = QuadratureRule{RefQuadrilateral}(2)
cellvalues = CellValues(qr, ip);

# Далі нам потрібно визначити DofHandler, який візьме на себе нумерацію
# і розподіл ступенів свободи для наших наближених полів.
# Потім додаємо одне скалярне поле з ім'ям :u
# на основі нашого ip інтерполяції, визначеного вище.
# Нарешті закриваємо функцією close! наш DofHandler.
# Тепер dof розподілені для всіх елементів.
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh);

# Ця функція повертає розріджену матрицю
K = allocate_matrix(dh)

# Граничні умови оброблюються обробником ConstraintHandler.
ch = ConstraintHandler(dh);

# Далі нам потрібно додати обмеження до ch.
# Для цього завдання ми визначаємо однорідні граничні умови Диріхле
∂Ω = union(
    getfacetset(grid, "left"),
    getfacetset(grid, "right"),
    getfacetset(grid, "top"),
    getfacetset(grid, "bottom"),
);

dbc = Dirichlet(:u, ∂Ω, (x, t) -> 0)
add!(ch, dbc);

# Нарешті нам також потрібно застосувати функцію close!
# Коли ми викликаємо close! dofs, що відповідають нашим
# обмеженням, обчислюються і зберігаються в нашому об'єкті ch.
close!(ch)


# Збірка лінійної системи
# Складання елементів
# Ми визначаємо функцію assemble_element! (Див. нижче),
# яка обчислює вклад кожного елемента.
function assemble_element!(Ke::Matrix, fe::Vector, cellvalues::CellValues)
    n_basefuncs = getnbasefunctions(cellvalues)
    # Reset to 0
    fill!(Ke, 0)
    fill!(fe, 0)
    # Цикл по квадратурних точках
    for q_point in 1:getnquadpoints(cellvalues)
        # Видобування квадратурних мір
        dΩ = getdetJdV(cellvalues, q_point)
        # Цикл по тестовим функціям
        for i in 1:n_basefuncs
            δu  = shape_value(cellvalues, q_point, i)
            ∇δu = shape_gradient(cellvalues, q_point, i)
            # Додавання розподілу по fe
            fe[i] += δu * dΩ
            # Цикл по триал-функціям
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q_point, j)
                # Додавання розподілу по Ke
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
    return Ke, fe
end

# Глобальна збірка матриць
# Ми визначаємо функцію assemble_global для перебору елементів
# і виконання глобальної збірки
function assemble_global(cellvalues::CellValues, K::SparseMatrixCSC, dh::DofHandler)
    # Створення матриці жорсткості елемента та вектора сили елемента
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    # Створення глобального вектору сили
    f = zeros(ndofs(dh))
    # Створення асемблера
    assembler = start_assemble(K, f)
    # Цикл по всім ячійкам
    for cell in CellIterator(dh)
        # Инициалізація cellvalues
        reinit!(cellvalues, cell)
        # Обчислення розподілу елементів
        assemble_element!(Ke, fe, cellvalues)
        # Збірка матриць Ke і fe в K і f
        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return K, f
end

# Отримання рішення системи
# Останній крок – рішення нашої системи рівняннь. Спочатку ми викликаємо збірник аssemble_global
# для отримання глобальної матриці жорсткості K та вектора сили f.
K, f = assemble_global(cellvalues, K, dh);

# Для обліку граничних умов ми використовуємо функцію apply!
# Це змінює елементи K і f відповідно таким чином,
# що ми можемо отримати правильний вектор переміщень u.

apply!(K, f, ch)
u = K \ f;

# Експорт рішення в файл VTK
# Щоб візуалізувати результат, ми експортуємо сітку та наше поле u
# у VTK-файл, який можна переглянути, наприклад, у ParaView.
VTKGridFile("heat_equation_ $dx x $dy", dh) do vtk
    write_solution(vtk, dh, u)
end