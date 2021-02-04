    program Approach_function

   ! implicit none ! отключаем неявную типизацию

    ! Variables
    real :: rLeft, rRight, h ! левая и правая границы отрезка, шаг разбиения
    integer :: N, n0 ! количество разбиений, степень многочлена
    real, dimension (:), allocatable :: X, Y, Polinomial, Vector, Grad, TmpPolinomial ! динамические массивы вещественных чисел
    integer :: i, m = 2, k ! счетчик, количество многочленов
    logical :: bOK ! логическая переменная, True, если метод спуска только начинается, т.е. считаем v(0), иначе - False
    real :: Del, eps = 0.0001 ! дельта, эпсилон
    real :: Sigma ! сигма - среднеквадратическое отклонение
    
    ! вводим данные
    write(*,*) 'Введите левую границу отрезка: '
    read *, rLeft
    write(*,*) 'Введите правую границу отрезка: '
    read *, rRight
    write(*,*) 'Введите количество разбиений: '
    read *, N
    write(*,*) 'Введите степень n0: '
    read *, n0
    
    
    h = (rRight - rLeft)/N ! шаг разбиения
    
    allocate(X(N+1)) ! выделяем место под массив X
    
    ! задаем X(i)
    do i = 1, N+1
        X(i) = rLeft + (i-1)*h
    end do
    
    write(*,*) 'Выберите функцию: '
    write(*,*) '1) Y = 3*cos(X) + X*sin(X*X)'
    write(*,*) '2) Y = (X - 1)*(X - 0.5)*(X - 0.25)'
    write(*,*) '3) Y = (X - 1)*(X - 0.5)*(X - 0.25)**2'
    write (*,*) 'Выбрана функция: '
    read *, i
    
    allocate(Y(N+1)) ! выделяем место под массив Y
    ! выбираем функцию, которую будем аппроксимировать
    select case (i)
        case (1)
            call Value (X, Y, N+1) ! вызываем функцию для вычисления значений функции 1
        case (2)
            call ValuePol1 (X, Y, N+1) ! вызываем функцию для вычисления значений функции 2 (многочлен 1)
        case (3)
            call ValuePol2 (X, Y, N+1) ! вызываем функцию для вычисления значений функции 3 (многочлена 2)
    end select
    
    
    allocate(Polinomial(n0+m)) ! выделяем место под массив Polinomial - многочлен
    allocate(Vector(n0+m)) ! выделяем место под массив Vector
    allocate(Grad(n0+m)) ! выделяем место под массив Grad
    allocate(TmpPolinomial(n0+m)) ! выделяем место под массив TmpPolinomial 
    
    ! обнуляем Polinomial, Grad, Vector, TmpPolinomial
    do i = 1, size(Polinomial)
        Polinomial(i) = 0
        Grad(i) = 0
        Vector(i) = 0
        TmpPolinomial(i) = 0
    end do
    
    !call FastDown(N, n0, m, bOK, TmpPolinomial, Polinomial, Vector)
    
    ! цикл по многочленам, которые строим
    do k = n0, n0 + m - 1
            i = 0
            bOK = .TRUE. 
            Del = 1000
            
            write(*,*), 'n: ', k     ! печатаем степень многочлена
            
            ! цикл пока дельта к-ое >= эпсилон, то есть пока идет метод спуска
            do while (Del >= eps)
                TmpPolinomial = Polinomial ! сохраняем в TmpPolinomial многочлен
                !write(*,*) 'TmpPolinomial = ', TmpPolinomial
                
                i = i + 1
                call GetVector() ! вызываем функцию построения вектора
                !write(*,*) 'Vector = ', Vector
                Polinomial = Polinomial + GetLambda()*Vector ! строим многочлен
                
                ! если метод начался, уже посчитано v(0), то меняем значение логической переменной, иначе - считаем дельту
                if (bOK.eqv..TRUE.) then
                    bOK = .FALSE.
                else
                    Del = GetDel() ! считаем дельту
                end if
                write(*,*) 'Number of iteration: ', i, ': |grad(F)| = ', sqrt(dot_product(Grad, Grad))
               ! печатаем номер итерации и модуль градиента
            end do
            
           
           write(*,*) 'iteration     ', 'X(i)     ', 'Y(i)     ', 'P(i)     ', '|Y(X(i)) - Pn(X(i))|     '    ! печатаем X, Y, P, |Y-P|
            
           do i = 1, N+1
                write(*,*) i-1, X(i), Y(i), Pol(X(i), Polinomial), abs(Y(i) - Pol(X(i), Polinomial))             ! печатаем итерацию, X, Y, P, |Y-P|
           end do
                
           write(*,*) 'A0...An'
           
           ! печатаем коэффициенты многочлена
           do i = 0, k
               write(*,*) 'A', i, ' = ', Polinomial(i+1)
           end do
           
           ! вычисляем и печатаем сигму
           Sigma = sqrt(F(Polinomial)/(N+1))
           write(*,*) 'Sigma = ', Sigma     
       end do        
    
    ! очищаем память
    deallocate(X)
    deallocate(Y)
    deallocate(Polinomial)
    deallocate(Vector)
    deallocate(Grad)
    deallocate(TmpPolinomial)
    
 
    !pause
contains
    
    ! процедура вычисления значения функции в точке
    subroutine Value(X, Y, N)  
        integer :: N, i ! количество разбиений, индекс
        real, dimension (N) :: X, Y ! вещественный массивы X, Y из N элементов
      
        do i = 1, N+1
            Y(i) = 3*cos(X(i)) + X(i)*sin(X(i)*X(i)) ! вычисляем значение функции в X(i)-й точке
        end do
    end subroutine Value
    
    ! процедура вычисления значения функции (многочлена 1) в точке
    subroutine ValuePol1(X, Y, N)
        integer :: N, i ! количество разбиений, индекс
        real, dimension (N) :: X, Y ! вещественный массивы X, Y из N элементов
      
        do i = 1, N+1
            Y(i) = (X(i) - 1)*(X(i) - 0.5)*(X(i) - 0.25) ! вычисляем значение функции в X(i)-й точке
        end do
    
    end subroutine ValuePol1
    
    ! процедура вычисления значения функции (мноочлена 2) в точке
    subroutine ValuePol2(X, Y, N)
        integer :: N, i ! количество разбиений, индекс
        real, dimension (N) :: X, Y ! вещественный массивы X, Y из N элементов
      
        do i = 1, N+1
            Y(i) = (X(i) - 1)*(X(i) - 0.5)*(X(i) - 0.25)**2 ! вычисляем значение функции в X(i)-й точке
        end do
    
    end subroutine ValuePol2
    
    
    ! функция - построение полинома
    real function Pol(variable, Array)
        real :: variable ! переменная в многочлене
        real, dimension (:) :: Array ! массив вещественных чисел
        integer :: i ! индекс
        
        Pol = 0 
        
        do i = 0, k
            Pol = Pol*variable + Array(k - i + 1)
        end do
    end function Pol
    
    
    ! функция F - построение Ф(а)
    real function F(Array)
        real, dimension (:) :: Array ! динамический вещественный массив
        integer :: i ! индекс
        F = 0
        
        do i = 1, N+1
            F = F + (Y(i) - Pol(X(i), Array))**2
        end do
    end function F
    
    
    ! функция получения градиента
    subroutine GetGrad ()
        integer :: i, j ! индексы
        
        do i = 1, k+1 ! общая сумма (по всем Ф')
            Grad(i) = 0
            do j = 1, N+1 ! по всем Х
                Grad(i) = Grad(i) + 2*(X(j)**(i-1))*(Pol(X(j),Polinomial) - Y(j))
            end do
        end do
    end subroutine GetGrad
    
    
    ! функция получения вектора V
    subroutine GetVector()
        real(8) :: rBeta, rTmp 
        
        ! если метод только начался, то считаем вектор v(0)
        if (bOK.eqv..TRUE.) then
            call GetGrad() ! вычисляем градиент
            Vector = (-1)*Grad ! считаем V(0)
            else
                rTmp = dot_product(Grad, Grad) ! сохраняем модуль старого градиента в квадрате
                call GetGrad() ! вычисляем градиент
                rBeta = dot_product(Grad,Grad)/rTmp ! считаем бета
                Vector = (-1)*Grad + rBeta*Vector ! считаем V(K+1)
        end if
    
    end subroutine GetVector
    
    
    ! функция получения лямбда к-го
    real function GetLambda()
        !implicit none
        real :: M1, M2, M3
    
        M1 = F(Polinomial - Vector)
        M2 = F(Polinomial)
        M3 = F(Polinomial + Vector)
        
        GetLambda = (M1 - M3)/(2*(M1 - 2*M2 + M3))
    
    end function GetLambda
    
    
    ! фунция получения дельта к-го
    real function GetDel()
        !implicit none
        integer :: i
        real(8) :: rDelta
        
        
        TmpPolinomial = TmpPolinomial - Polinomial ! записали числитель формулы для дельты к-го ( Ai(k) - Ai(k-1) ) 
        GetDel = 1000
        
        do i = 1, k+1
            ! если i-тая компонента полинома (Ai) не равно 0, то считаем дельту
            if (Polinomial(i) /= 0) then
                rDelta = abs(TmpPolinomial(i)/Polinomial(i)) ! считаем дельту
                if (GetDel > rDelta) then
                    GetDel = rDelta ! ищем максимальное дельта к-ое
                endif
            endif
        end do
       ! write(*,*) 'rDelta = ', rDelta
    end function GetDel
    
    
    end program Approach_function