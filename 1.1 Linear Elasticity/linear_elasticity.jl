# Завантажуемо необхідні біліотеки та залежності
using Ferrite, FerriteGmsh, SparseArrays

# В якості розрахункової області будемо використовувати квадратну область зі стороною рівною 1
# Будемо використовувати готову сітку, яку завантажимо за допомогою біліотеки FerriteGmsh.jl
using Downloads: download
logo_mesh = "logo.geo"
asset_url = "https://raw.githubusercontent.com/Ferrite-FEM/Ferrite.jl/gh-pages/assets/"
isfile(logo_mesh) || download(string(asset_url, logo_mesh), logo_mesh)
FerriteGmsh.Gmsh.initialize()
FerriteGmsh.Gmsh.gmsh.option.set_number("General.Verbosity", 2)
grid = togrid(logo_mesh);
FerriteGmsh.Gmsh.finalize();

# У згенерованій сітці відсутні фасети для меж розрахункової області,
# тому ми додаємо їх за допомогою функції addfacetset з Ferrite.jl
# Він дозволяє нам додавати фасети до сітки на основі координат.
# Зверніть увагу, що приблизне порівняння з 0.0 не дуже добре працює, тому ми використовуємо допуск.
addfacetset!(grid, "top",    x -> x[2] ≈ 1.0) 
addfacetset!(grid, "left",   x -> abs(x[1]) < 1e-6)
addfacetset!(grid, "bottom", x -> abs(x[2]) < 1e-6);

dim = 2                                             # 2D-задача
order = 1                                           # лінійна інтерполяція
ip = Lagrange{RefTriangle, order}()^dim;            # векторна інтерполяція

qr = QuadratureRule{RefTriangle}(1)                 # перша квадратурна точка
qr_face = FacetQuadratureRule{RefTriangle}(1);      # правило квадратур для фасетів

cellvalues = CellValues(qr, ip)                     # збірка набору ячійок
facetvalues = FacetValues(qr_face, ip);             # збірка набору фасетів

# СТУПЕНІ СВОБОДИ
# Для розподілу ступенів свободи ми визначаємо DofHandler.
# DofHandler знає, що u має два ступені свободи для кожного вузла,
# тому що ми векторизували інтерполяцію.
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh);

# ГРАНИЧНІ УМОВИ ДИРИХЛЕ ТА НЕЙМАНА
# Ми задаємо граничні умови Діріхлe, фіксуючи нормаль руху
# до нижньої та лівої границі розрахункової області. Останній аргумент Dirichlet визначає,
# які компоненти поля повинні бути обмежені. Якщо аргумент не вказано,
# всі компоненти обмежені за замовченням.
ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> 0.0, 2))
add!(ch, Dirichlet(:u, getfacetset(grid, "left"),   (x, t) -> 0.0, 1))
close!(ch);

# Граничні умови Неймана
# На верхнє ребро розрахункової області додамо розподілене
# навантаження tN(x) = (20e3)x1e2 (N/mm2).
traction(x) = Vec(0.0, 20e3 * x[1]);

# ІТЕРАТИВНИЙ ПЕРЕБІР ВСІХ ФАСЕТ У НАБОРІ
# На правій границі розрахунковї області ми нічого не робимо, в результаті чого отримуємо
# граничну умову Неймана з нульовим значенням.
# Для того, щоб зібрати зовнішнє зусилля, що діє на верхню лінію
# нашої розрахункової області, необхідно перебрати всі фасети у відомому
# набірі фасетів. Це робиться за допомогою наступної функції:
function assemble_external_forces!(f_ext, dh, facetset, facetvalues, prescribed_traction)
    # Створення тимчасового масиву для локальних вкладів грані у вектор зовнішньої сили
    fe_ext = zeros(getnbasefunctions(facetvalues))
    for facet in FacetIterator(dh, facetset)
        # Оновлення значення фасетів до правильного номера фасета
        reinit!(facetvalues, facet)
        # Тимчасовий масив для наступного фасета
        fill!(fe_ext, 0.0)
        # Доступ до координат ячійок (cells)
        cell_coordinates = getcoordinates(facet)
        for qp in 1:getnquadpoints(facetvalues)
            # Обчислення глобальних координат квадратурних точок.
            x = spatial_coordinate(facetvalues, qp, cell_coordinates)
            tₚ = prescribed_traction(x)
            # Отримання вісів інтергрування для поточної квадратурної точки.
            dΓ = getdetJdV(facetvalues, qp)
            for i in 1:getnbasefunctions(facetvalues)
                Nᵢ = shape_value(facetvalues, qp, i)
                fe_ext[i] += tₚ ⋅ Nᵢ * dΓ
            end
        end
        # Додавання  локальних внесків до правильних індексів у векторі глобальної зовнішньої сили
        assemble!(f_ext, celldofs(facet), fe_ext)
    end
    return f_ext
end

# Параметри матеріалу
Emod = 190e3 # Модуль Юнга [MPa]
ν = 0.27      # Коєфіціент Пуасона [-]

Gmod = Emod / (2(1 + ν))  # Модуль зсуву (Shear modulus)
Kmod = Emod / (3(1 - 2ν)) # Об'ємний модуль (Bulk modulus)

# Далі, запускаємо автоамтичне диференціювання за допомогою модуля Tensors.jl
# для розрахунку тензора пружної жорсткості C
C = gradient(ϵ -> 2 * Gmod * dev(ϵ) + 3 * Kmod * vol(ϵ), zero(SymmetricTensor{2,2}));

# ОБЧИСЛЕННЯ ГЛОБАЛЬНОЇ МАТРИЦІ ЖОРСТКОСТІ
# Щоб розрахувати глобальну матрицю жорсткості Kij, створимо
# спеціальну функцію, яка обчислює локальні матриці жорсткості ke
# для кожного елемента і потім збирає з них глобальну матрицю жорсткості.
# Tензор пружної жорсткості C є постійним.
# Таким чином, його необхідно обчислити лише один раз і потім
# використовувати для всіх точок інтегрування.
function assemble_cell!(ke, cellvalues, C)
    for q_point in 1:getnquadpoints(cellvalues)
        # Отримання вісів інтергрування для квадратурної точки
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:getnbasefunctions(cellvalues)
            # Градієнт тестової функції
            ∇Nᵢ = shape_gradient(cellvalues, q_point, i)
            for j in 1:getnbasefunctions(cellvalues)
                # Симетричний градіент для trial-функції
                ∇ˢʸᵐNⱼ = shape_symmetric_gradient(cellvalues, q_point, j)
                ke[i, j] += (∇Nᵢ ⊡ C ⊡ ∇ˢʸᵐNⱼ) * dΩ
            end
        end
    end
    return ke
end

# ГЛОБАЛЬНA ЗБІРКА ЗАДАЧИ
# Ми визначаємо функцію assemble_global для перебору елементів та
# виконання глобальної збірки. Функція приймає попередньо
# виділену розріджену матрицю K, наш DofHandler dh, наші cellvalues
# і тензор жорсткості С як вхідні аргументи і обчислює глобальну матрицю жорсткості K:
function assemble_global!(K, dh, cellvalues, C)
    # Виділиння матриці жорсткості елемента
    n_basefuncs = getnbasefunctions(cellvalues)
    ke = zeros(n_basefuncs, n_basefuncs)
    # Створення ассемблера матриць
    assembler = start_assemble(K)
    # Цикл проходження по всім ячійкам розрахункової сітки
    for cell in CellIterator(dh)
        # Оновлення градіентів функцій форми на основі координат ячійок
        reinit!(cellvalues, cell)
        # Матриця жорсткості елемента
        fill!(ke, 0.0)
        # Обчислення розподілу елементів
        assemble_cell!(ke, cellvalues, C)
        # Збірка глобальної матриці жорскості K
        assemble!(assembler, celldofs(cell), ke)
    end
    return K
end

# РІШЕННЯ
# Останній крок – рішення системи. Спочатку виділяємо глобальну матрицю жорсткості K і збираємо її.
K = allocate_matrix(dh)
assemble_global!(K, dh, cellvalues, C);

# Потім виділяємо та збираємо вектор зовнішньої сили.
f_ext = zeros(ndofs(dh))
assemble_external_forces!(f_ext, dh, getfacetset(grid, "top"), facetvalues, traction);

# Для обліку граничних умов Діріхле ми використовуємо функцію apply!
# Це змінює елементи в матрицях K і f таким чином, що ми можемо отримати
# коректний вектор рішення u, використовуючи рішення системи алгебраїчних лінійних рівнянь.
apply!(K, f_ext, ch)
u = K \ f_ext;

# # ПОСТРАБОТКА
# У цьому випадку ми хочемо проаналізувати переміщення u,
# а також поле напруженнь. Розраховуємо напрження в кожній квадратурній
# точці, а потім експортуємо його двома різними способами:
# 1. Беремо постійну величину в кожній ячійці (збігається з наближенням постійних
# деформацій в кожному елементі). При цьому слід зазначити, що поточне
# обмеження полягає в тому, що дані ячійок для тензорів другого
# порядку повинні бути експортовані за компонентами.
# 2. Інтерполяція за допомогою функцій лінійного анзацу Лагранжа через L2Projector
function calculate_stresses(grid, dh, cv, u, C)
    qp_stresses = [
        [zero(SymmetricTensor{2,2}) for _ in 1:getnquadpoints(cv)]
        for _ in 1:getncells(grid)]
    avg_cell_stresses = tuple((zeros(getncells(grid)) for _ in 1:3)...)
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        cell_stresses = qp_stresses[cellid(cell)]
        for q_point in 1:getnquadpoints(cv)
            ε = function_symmetric_gradient(cv, q_point, u, celldofs(cell))
            cell_stresses[q_point] = C ⊡ ε
        end
        σ_avg = sum(cell_stresses) / getnquadpoints(cv)
        avg_cell_stresses[1][cellid(cell)] = σ_avg[1, 1]
        avg_cell_stresses[2][cellid(cell)] = σ_avg[2, 2]
        avg_cell_stresses[3][cellid(cell)] = σ_avg[1, 2]
    end
    return qp_stresses, avg_cell_stresses
end

qp_stresses, avg_cell_stresses = calculate_stresses(grid, dh, cellvalues, u, C);

# Тепер ми використовуємо L2Projector для проектування поля напруженнь
# на шматково-лінійний простір кінцевих елементів, яке ми використовували для розв'язання задачі.
proj = L2Projector(Lagrange{RefTriangle, 1}(), grid)
stress_field = project(proj, qp_stresses, qr);

color_data = zeros(Int, getncells(grid))         
colors = [                                       
    "1" => 1, "5" => 1, # purple                 
    "2" => 2, "3" => 2, # red                    
    "4" => 3,           # blue                   
    "6" => 4            # green                  
    ]                                            
for (key, color) in colors                       
    for i in getcellset(grid, key)               
        color_data[i] = color                    
    end                                          
end                                              

# Для візуалізації результату ми експортуємо наші результати до VTK-файлу.
# Створюється неструктурований файл *.vtu, який можна
# переглянути, наприклад, у ParaView.
VTKGridFile("linear_elasticity", dh) do vtk
    write_solution(vtk, dh, u)
    for (i, key) in enumerate(("11", "22", "12"))
        write_cell_data(vtk, avg_cell_stresses[i], "sigma_" * key)
    end
    write_projection(vtk, proj, stress_field, "stress field")
    Ferrite.write_cellset(vtk, grid)
    write_cell_data(vtk, color_data, "colors") #hide
end