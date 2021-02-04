    program Approach_function

   ! implicit none ! ��������� ������� ���������

    ! Variables
    real :: rLeft, rRight, h ! ����� � ������ ������� �������, ��� ���������
    integer :: N, n0 ! ���������� ���������, ������� ����������
    real, dimension (:), allocatable :: X, Y, Polinomial, Vector, Grad, TmpPolinomial ! ������������ ������� ������������ �����
    integer :: i, m = 2, k ! �������, ���������� �����������
    logical :: bOK ! ���������� ����������, True, ���� ����� ������ ������ ����������, �.�. ������� v(0), ����� - False
    real :: Del, eps = 0.0001 ! ������, �������
    real :: Sigma ! ����� - �������������������� ����������
    
    ! ������ ������
    write(*,*) '������� ����� ������� �������: '
    read *, rLeft
    write(*,*) '������� ������ ������� �������: '
    read *, rRight
    write(*,*) '������� ���������� ���������: '
    read *, N
    write(*,*) '������� ������� n0: '
    read *, n0
    
    
    h = (rRight - rLeft)/N ! ��� ���������
    
    allocate(X(N+1)) ! �������� ����� ��� ������ X
    
    ! ������ X(i)
    do i = 1, N+1
        X(i) = rLeft + (i-1)*h
    end do
    
    write(*,*) '�������� �������: '
    write(*,*) '1) Y = 3*cos(X) + X*sin(X*X)'
    write(*,*) '2) Y = (X - 1)*(X - 0.5)*(X - 0.25)'
    write(*,*) '3) Y = (X - 1)*(X - 0.5)*(X - 0.25)**2'
    write (*,*) '������� �������: '
    read *, i
    
    allocate(Y(N+1)) ! �������� ����� ��� ������ Y
    ! �������� �������, ������� ����� ����������������
    select case (i)
        case (1)
            call Value (X, Y, N+1) ! �������� ������� ��� ���������� �������� ������� 1
        case (2)
            call ValuePol1 (X, Y, N+1) ! �������� ������� ��� ���������� �������� ������� 2 (��������� 1)
        case (3)
            call ValuePol2 (X, Y, N+1) ! �������� ������� ��� ���������� �������� ������� 3 (���������� 2)
    end select
    
    
    allocate(Polinomial(n0+m)) ! �������� ����� ��� ������ Polinomial - ���������
    allocate(Vector(n0+m)) ! �������� ����� ��� ������ Vector
    allocate(Grad(n0+m)) ! �������� ����� ��� ������ Grad
    allocate(TmpPolinomial(n0+m)) ! �������� ����� ��� ������ TmpPolinomial 
    
    ! �������� Polinomial, Grad, Vector, TmpPolinomial
    do i = 1, size(Polinomial)
        Polinomial(i) = 0
        Grad(i) = 0
        Vector(i) = 0
        TmpPolinomial(i) = 0
    end do
    
    !call FastDown(N, n0, m, bOK, TmpPolinomial, Polinomial, Vector)
    
    ! ���� �� �����������, ������� ������
    do k = n0, n0 + m - 1
            i = 0
            bOK = .TRUE. 
            Del = 1000
            
            write(*,*), 'n: ', k     ! �������� ������� ����������
            
            ! ���� ���� ������ �-�� >= �������, �� ���� ���� ���� ����� ������
            do while (Del >= eps)
                TmpPolinomial = Polinomial ! ��������� � TmpPolinomial ���������
                !write(*,*) 'TmpPolinomial = ', TmpPolinomial
                
                i = i + 1
                call GetVector() ! �������� ������� ���������� �������
                !write(*,*) 'Vector = ', Vector
                Polinomial = Polinomial + GetLambda()*Vector ! ������ ���������
                
                ! ���� ����� �������, ��� ��������� v(0), �� ������ �������� ���������� ����������, ����� - ������� ������
                if (bOK.eqv..TRUE.) then
                    bOK = .FALSE.
                else
                    Del = GetDel() ! ������� ������
                end if
                write(*,*) 'Number of iteration: ', i, ': |grad(F)| = ', sqrt(dot_product(Grad, Grad))
               ! �������� ����� �������� � ������ ���������
            end do
            
           
           write(*,*) 'iteration     ', 'X(i)     ', 'Y(i)     ', 'P(i)     ', '|Y(X(i)) - Pn(X(i))|     '    ! �������� X, Y, P, |Y-P|
            
           do i = 1, N+1
                write(*,*) i-1, X(i), Y(i), Pol(X(i), Polinomial), abs(Y(i) - Pol(X(i), Polinomial))             ! �������� ��������, X, Y, P, |Y-P|
           end do
                
           write(*,*) 'A0...An'
           
           ! �������� ������������ ����������
           do i = 0, k
               write(*,*) 'A', i, ' = ', Polinomial(i+1)
           end do
           
           ! ��������� � �������� �����
           Sigma = sqrt(F(Polinomial)/(N+1))
           write(*,*) 'Sigma = ', Sigma     
       end do        
    
    ! ������� ������
    deallocate(X)
    deallocate(Y)
    deallocate(Polinomial)
    deallocate(Vector)
    deallocate(Grad)
    deallocate(TmpPolinomial)
    
 
    !pause
contains
    
    ! ��������� ���������� �������� ������� � �����
    subroutine Value(X, Y, N)  
        integer :: N, i ! ���������� ���������, ������
        real, dimension (N) :: X, Y ! ������������ ������� X, Y �� N ���������
      
        do i = 1, N+1
            Y(i) = 3*cos(X(i)) + X(i)*sin(X(i)*X(i)) ! ��������� �������� ������� � X(i)-� �����
        end do
    end subroutine Value
    
    ! ��������� ���������� �������� ������� (���������� 1) � �����
    subroutine ValuePol1(X, Y, N)
        integer :: N, i ! ���������� ���������, ������
        real, dimension (N) :: X, Y ! ������������ ������� X, Y �� N ���������
      
        do i = 1, N+1
            Y(i) = (X(i) - 1)*(X(i) - 0.5)*(X(i) - 0.25) ! ��������� �������� ������� � X(i)-� �����
        end do
    
    end subroutine ValuePol1
    
    ! ��������� ���������� �������� ������� (��������� 2) � �����
    subroutine ValuePol2(X, Y, N)
        integer :: N, i ! ���������� ���������, ������
        real, dimension (N) :: X, Y ! ������������ ������� X, Y �� N ���������
      
        do i = 1, N+1
            Y(i) = (X(i) - 1)*(X(i) - 0.5)*(X(i) - 0.25)**2 ! ��������� �������� ������� � X(i)-� �����
        end do
    
    end subroutine ValuePol2
    
    
    ! ������� - ���������� ��������
    real function Pol(variable, Array)
        real :: variable ! ���������� � ����������
        real, dimension (:) :: Array ! ������ ������������ �����
        integer :: i ! ������
        
        Pol = 0 
        
        do i = 0, k
            Pol = Pol*variable + Array(k - i + 1)
        end do
    end function Pol
    
    
    ! ������� F - ���������� �(�)
    real function F(Array)
        real, dimension (:) :: Array ! ������������ ������������ ������
        integer :: i ! ������
        F = 0
        
        do i = 1, N+1
            F = F + (Y(i) - Pol(X(i), Array))**2
        end do
    end function F
    
    
    ! ������� ��������� ���������
    subroutine GetGrad ()
        integer :: i, j ! �������
        
        do i = 1, k+1 ! ����� ����� (�� ���� �')
            Grad(i) = 0
            do j = 1, N+1 ! �� ���� �
                Grad(i) = Grad(i) + 2*(X(j)**(i-1))*(Pol(X(j),Polinomial) - Y(j))
            end do
        end do
    end subroutine GetGrad
    
    
    ! ������� ��������� ������� V
    subroutine GetVector()
        real(8) :: rBeta, rTmp 
        
        ! ���� ����� ������ �������, �� ������� ������ v(0)
        if (bOK.eqv..TRUE.) then
            call GetGrad() ! ��������� ��������
            Vector = (-1)*Grad ! ������� V(0)
            else
                rTmp = dot_product(Grad, Grad) ! ��������� ������ ������� ��������� � ��������
                call GetGrad() ! ��������� ��������
                rBeta = dot_product(Grad,Grad)/rTmp ! ������� ����
                Vector = (-1)*Grad + rBeta*Vector ! ������� V(K+1)
        end if
    
    end subroutine GetVector
    
    
    ! ������� ��������� ������ �-��
    real function GetLambda()
        !implicit none
        real :: M1, M2, M3
    
        M1 = F(Polinomial - Vector)
        M2 = F(Polinomial)
        M3 = F(Polinomial + Vector)
        
        GetLambda = (M1 - M3)/(2*(M1 - 2*M2 + M3))
    
    end function GetLambda
    
    
    ! ������ ��������� ������ �-��
    real function GetDel()
        !implicit none
        integer :: i
        real(8) :: rDelta
        
        
        TmpPolinomial = TmpPolinomial - Polinomial ! �������� ��������� ������� ��� ������ �-�� ( Ai(k) - Ai(k-1) ) 
        GetDel = 1000
        
        do i = 1, k+1
            ! ���� i-��� ���������� �������� (Ai) �� ����� 0, �� ������� ������
            if (Polinomial(i) /= 0) then
                rDelta = abs(TmpPolinomial(i)/Polinomial(i)) ! ������� ������
                if (GetDel > rDelta) then
                    GetDel = rDelta ! ���� ������������ ������ �-��
                endif
            endif
        end do
       ! write(*,*) 'rDelta = ', rDelta
    end function GetDel
    
    
    end program Approach_function