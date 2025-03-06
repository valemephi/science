!main7 - нерабочий
program task11

implicit none 
    real ( kind = 8 ) bi0
    real ( kind = 8 ) bi1
    real ( kind = 8 ) bk0
    real ( kind = 8 ) bk1
    real ( kind = 8 ) di0
    real ( kind = 8 ) di1
    real ( kind = 8 ) dk0
    real ( kind = 8 ) dk1

    real(kind = 8), parameter :: k1 = 1.38 * 10.0**(-23)
    real(kind = 8), parameter :: e =  2.71828182845904, c = 3e8
    real(kind = 8), parameter :: phi0 =  2.06783460999999E-07            !Квант магнитного потока
    real(kind = 8), parameter :: delta = 0.27       !Толщина слоя
    real(kind = 8) :: xi0, xi          !Длина когерентности сверхпроводника
    real(kind = 8), parameter :: lambda0 = 180.0         !Лондоновская глубина проникновения данного материала при Т=0 
    real(kind = 8) :: lambda             !Лондоновская глубина проникновения данного материала при данной температуре
    real(kind = 8) :: T, Tc
    real(kind = 8) :: length, width
    real(kind = 8) :: pi
    real(kind=8) :: tauNa !параметр для процессов рождения и уничтожения,
                          !меняя его, можно регулировать скорость сходимости алгоритмов Монте-Карло
    integer:: iter, N1, N !N- число вихрей в системе, N1- регулируем крайние случаи, iter-число шагов
    real(kind=8):: K0tab(10000), epss, U0 !табулированный массив энергий, собств.энергия вихря и U0
    real(kind=8):: H0, HI, I_ ! Поля и токи

    real(kind=8):: tabul, step, xtabul !Для табулирования
    integer:: itabul
    
    real(kind=8):: Uj, Uj_end, rij, Um_Ut, Usurf
    real(kind=8):: poss1, poss2, poss3, temp, process
    integer:: i, i1, j, k, i_dist
    integer:: i_curve !для петли намагниченности

    real(kind=8):: r1plus(10000,2), r1minus(10000, 2) !массивы вихрей и антивихрей
    integer:: Nplus=1, Nminus=1
    real(kind=8) :: minus_1, plus_1, flag1
    plus_1 = 1.
    minus_1 = -1.
    flag1=1.
    
    tauNa = 100.0
    T = 5.0
    xi0 = 2.0
    iter = 10000

    r1plus(:,:)=100000.
    r1minus(:,:)=100000.
    ! endd = 100

    poss1 = 0.4
    poss2 = 0.65
    poss3 = 0.9

    tabul = 10000.0
    length = 2.0*10.**3 ! по х
    width = 2.0*10.**3   !по у
    Tc = 84.0
    H0 = 0.0
    I_ = 0.

    open(2,file = "res1.txt")
    ! open(33,file = "resanim.txt")
    ! open(4, file = "task1energy.txt")
    open(45, file = "magncurve2.txt")

    xi = xi0 * (1.0 - (T / Tc) ** 3.3) ** (-1.0 / 2.0)
    write(*,*) xi 
    lambda = lambda0 * (1.0 - (T / Tc) ** (3.3)) ** (-1.0 / 2.0)
    write(*,*) lambda
    pi = acos(-1.0)
    HI= 2*pi*I_/c
    write(*,*) H0,HI
    T= T / 11606.0 ! T в эВ
    
    !считаем собственные энергии вихрей
    epss = (phi0 / (4 * pi * lambda * 10 ** (-7.0))) ** 2 * (log(lambda0 / xi0) + 0.52)  * delta * 10 ** (-7.0) *6.2 * 10.0**(11)
    write(*,*) epss
    !Табулирование
    step = 7.0 * lambda / tabul
    do itabul = 1, floor(tabul)
        xtabul = step * itabul / lambda
        call besssel(xtabul,bi0, di0, bi1, di1, bk0, dk0, bk1, dk1)
        K0tab(itabul) = bk0
    end do
    U0 = (phi0**2.0 / (8.0 * pi**2.0 * (lambda * 10. ** (-7.0))**2.0)) * delta * 10.**(-7.0)
    K0tab = K0tab * 6.2 * 10.0 **(11.0)
    
    Uj = 0.0
    ! counter1 = 0
    do i_curve = 1, 1
    write(*,*) 'Iteration number ', i_curve, Nplus, Nminus
    !СЧИТАЕМ НАЧАЛЬНУЮ ЭНЕРГИЮ ТОЧКИ НА КРИВОЙ НАМАГНИЧЕННОСТИ
    Uj=epss*float(Nplus+Nminus-2)

    do j = 1, Nminus-1     !Считаем энергию антивихрей
        call Usurface(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
                Nplus, Nminus, width, length, minus_1, j)
        Uj = Uj + 1./2.*Usurf
        call Usurface1(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
            Nplus, Nminus, width, length, minus_1, j)
        Uj = Uj + 1./2.*Usurf
        call Ucurrent(r1minus(j,1), r1minus(j,1), H0, lambda, length, HI, Um_Ut, minus_1)
        Uj = Uj + Um_Ut*delta*10.**(-7.)
        do k = 1, Nminus-1    
            if( j .NE. k) then
                call distance1(r1minus(j,1), r1minus(k,1), r1minus(j,2), r1minus(k,2), width, rij)
                if(rij < 7.0 * lambda) then 
                    i_dist = ceiling(rij*tabul/(7.0*lambda)) 
            
                    Uj = Uj + 1./2.*U0*K0tab(i_dist)
                end if
            end if
        end do
        do k = 1, Nplus-1    !энергию начальной конфигурации
                call distance1(r1minus(j,1), r1plus(k,1), r1minus(j,2), r1plus(k,2), width, rij)
                if(rij < 7.0 * lambda) then 
                    i_dist = ceiling(rij*tabul/(7.0*lambda)) 
            
                    Uj = Uj - 1./2. *U0*K0tab(i_dist)             
                end if
        end do
    end do
    
    do j = 1, Nplus-1     !Считаем энергию вихрей
        call Usurface(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
                Nplus, Nminus, width, length, plus_1, j)
        Uj = Uj + 1./2.*Usurf
        call Usurface1(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
            Nplus, Nminus, width, length, plus_1, j)
        Uj = Uj + 1./2.*Usurf
        call Ucurrent(r1plus(j,1), r1plus(j,1), H0, lambda, length, HI, Um_Ut, plus_1)
        Uj = Uj + Um_Ut*delta*10.**(-7.)

        do k = 1, Nplus-1    
            if( j .NE. k) then
                call distance1(r1plus(j,1), r1plus(k,1), r1plus(j,2), r1plus(k,2), width, rij)
                if(rij < 7.0 * lambda) then 
                    i_dist = ceiling(rij*tabul/(7.0*lambda)) 
                    Uj = Uj + 1./2.*U0*K0tab(i_dist) 
                end if
            end if
        end do
        do k = 1, Nminus-1    !энергию начальной конфигурации
                call distance1(r1plus(j,1), r1minus(k,1), r1plus(j,2), r1minus(k,2), width, rij)
                if(rij < 7.0 * lambda) then 
                    i_dist = ceiling(rij*tabul/(7.0*lambda)) 
            
                    Uj = Uj - 1./2.*U0*K0tab(i_dist)
                end if
        end do
    end do




    if(H0<1000. .and.  flag1>0.) then
            H0 = H0+400.
        else
            flag1= -1.
            H0=H0-1000.
        end if
    !ВЫПОЛНЯЕМ АЛГОРИТМ ДЛЯ НОВОЙ ТОЧКИ
    write(*,*) Uj;
    ! H0 = H0 + 50.0
    do i = 1, iter
        N = Nplus + Nminus - 2
        if(N>10000) then
        go to 44
        end if
        if(N> 30 .and. N< 9999) then
            N1 = N
        else
            N1 = 30
        end if

        do i1 = 1, N1
            call random_number(process)
            if(process < -0.5) then
                call movement(Uj, lambda, T, length, width, epss, tabul, H0,&
                            HI, U0, K0tab, delta, r1minus, r1plus, Nminus, Nplus, xi)
            
            else if(process < 1.) then
                call creation(lambda, length, width, H0, HI, epss, Nminus, r1minus, &
                        Nplus, r1plus, tabul, delta, tauNa, T, Uj, U0, K0tab, xi)
                
            else if (Nplus+Nminus .NE. 2 .and. process< 0.9) then
                call annihilation(Uj, lambda, tauNa, T, K0tab, width, length, epss,&
                            H0, HI, U0, delta, tabul, r1plus, r1minus, Nplus, Nminus)

            else if (Nplus .NE. 1 .and. Nminus .NE. 1) then
                 call pair_annihilation(lambda, xi, T, tauNa, delta, width, length, tabul, &
                        epss,r1plus,r1minus,Nplus,Nminus, H0,HI, U0, K0tab)
            end if
        end do 
        
    !вывод энергии каждые сколько-то итераций
        if (mod(i,1000)==0) then
            write(*,*) i, 'Uj=', Uj
        end if
    !запись энергии в файл
        ! write(4,*) i, Uj
        ! if (mod(i,100)==0 .and. i<endd) then
        !     counter1 = counter1 + 1
        !     do ianim = 1, N-1
        !         write(33, *) r1(ianim, :)
        !     end do
        !     write(33,*) 'flag'
        ! end if
        
    end do
    if(H0>0.) then
        write(45,*) H0, H0 - float(Nplus+Nminus-2) * phi0 / (length * width * 10.0 ** (-14.0))
    else
        write(45,*) H0, H0 + float(Nplus+Nminus-2) * phi0 / (length * width * 10.0 ** (-14.0))
    end if
    write(*,*) H0, H0 - float(Nplus+Nminus-2) * phi0 / (length * width * 10.0 ** (-14.0));
    
    if(abs(H0+1000.)<50.) then
        flag1=1.
    end if
    end do
    !sectora = sectora / 1000.0
    !write(*,*) sectora
    close(45)
    !H0=H0-400.
44  write(*,*)'in the end', Uj
    Uj_end = 0.0
    ! write(*,*) N
    do j = 1, Nminus-1     !Считаем энергию антивихрей
        call Usurface(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
                Nplus, Nminus, width, length, minus_1, j)
        Uj_end = Uj_end + 1./2.*Usurf
        call Usurface1(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
            Nplus, Nminus, width, length, minus_1, j)
        Uj_end = Uj_end + 1./2.*Usurf
        call Ucurrent(r1minus(j,1), r1minus(j,1), H0, lambda, length, HI, Um_Ut, minus_1)
        Uj_end = Uj_end + Um_Ut*delta*10.**(-7.)
        do k = 1, Nminus-1    
            if( j .NE. k) then
                call distance1(r1minus(j,1), r1minus(k,1), r1minus(j,2), r1minus(k,2), width, rij)
                if(rij < 7.0 * lambda) then 
                    i_dist = ceiling(rij*tabul/(7.0*lambda)) 
            
                    temp = U0*K0tab(i_dist)
                    Uj_end = Uj_end + temp / 2.0
                
                end if
            end if
        end do
        do k = 1, Nplus-1    !энергию начальной конфигурации
                call distance1(r1minus(j,1), r1plus(k,1), r1minus(j,2), r1plus(k,2), width, rij)
                if(rij < 7.0 * lambda) then 
                    i_dist = ceiling(rij*tabul/(7.0*lambda)) 
            
                    temp = U0*K0tab(i_dist)
                    Uj_end = Uj_end - temp / 2.0
                
                end if
        end do
    end do
    
    do j = 1, Nplus-1     !Считаем энергию вихрей
        call Usurface(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
                Nplus, Nminus, width, length, plus_1, j)
        Uj_end = Uj_end + 1./2.*Usurf
        
        call Usurface1(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
            Nplus, Nminus, width, length, plus_1, j)
        Uj_end = Uj_end + 1./2.*Usurf
        call Ucurrent(r1plus(j,1), r1plus(j,1), H0, lambda, length, HI, Um_Ut, plus_1)
        Uj_end = Uj_end + Um_Ut*delta*10.**(-7.)

        do k = 1, Nplus-1    
            if( j .NE. k) then
                call distance1(r1plus(j,1), r1plus(k,1), r1plus(j,2), r1plus(k,2), width, rij)
                if(rij < 7.0 * lambda) then 
                    i_dist = ceiling(rij*tabul/(7.0*lambda)) 
                    Uj_end = Uj_end + 1./2.*U0*K0tab(i_dist) 
                end if
            end if
        end do
        do k = 1, Nminus-1    !энергию начальной конфигурации
                call distance1(r1plus(j,1), r1minus(k,1), r1plus(j,2), r1minus(k,2), width, rij)
                if(rij < 7.0 * lambda) then 
                    i_dist = ceiling(rij*tabul/(7.0*lambda)) 
            
                    temp = U0*K0tab(i_dist)
                    Uj_end = Uj_end - temp / 2.0
                
                end if
        end do
        ! write(*,*) 'Uij+-', Uj_end
    end do
    write(*,*) Uj_end +float(Nplus+Nminus-2) * (epss)
    !write(*,*) r1(1:11,:)
    do i = 1, Nplus-1
        write(2, *) r1plus(i, :)
    end do
    write(2, *) 'flag'
    do i = 1, Nminus-1
        write(2, *) r1minus(i, :)
    end do
    ! open(5, file="heatmap_task1_data.txt")
    ! do i = -200, 200
    !     nomera(i) = i
    ! end do
    ! write(5, *) nomera
    ! do i = -200, 200
    !     !write(5, *) i, sectora(i, :)
    ! end do
    ! close(5)
    
    close(2)
    ! close(33)
    ! close(4)
    !plot 'C:/Users/Valera P. Lenkov/Desktop/Fortran/task1/res2.txt'











    contains
    !Расстояние для одной периодичной и одной непериодичнойт границ
    subroutine distance1(x1, x2, y1, y2, width, rij)
        real(kind=8):: x1,x2,y1,y2,ydiff,xdiff,width,rij
        ydiff = min(abs(y1-y2), width - abs(y1 - y2))
        xdiff = abs(x1-x2)  
        rij = sqrt(ydiff**2 + xdiff**2)
    end subroutine

    !****************
    !Взаимодействие с токами
    !****************
    subroutine Ucurrent(side, x, H0, lambda, length, HI, Um_Ut, signum_)
        real(kind=8):: side, x, H0, lambda, length, HI, Um_Ut, pi, side1, signum_
        real(kind = 8), parameter :: phi0 =  2.06783460999999E-07
        pi = acos(-1.0)
        Um_Ut = 0.
        if (side>0.) then
            side1 = 1.
        else
            side1 = -1.
        end if
        Um_Ut = -(phi0/(4*pi))* &
        (H0*signum_*(1-cosh(x/lambda)/cosh(length/(2.*lambda)))&
        + HI*(sinh(x/lambda)/sinh(length/(2.*lambda))-side1))
        Um_Ut = Um_Ut * 6.2 * 10.**(11.)
        
    end subroutine
    !****************************
    !РОЖДЕНИЕ
    !****************************
    subroutine creation(lambda, length, width, H0, HI, epss, Nminus, r1minus, &
        Nplus, r1plus, tabul, delta, tauNa, T, Uj, U0, K0tab, xi)
        real(kind=8)::r1minus(10000,2), r1plus(10000,2)
        real(kind=8):: lambda, length, width, H0, HI,epss, T, delta, xi
        integer:: Nminus, Nplus

        real(kind = 8):: side, xcr, ycr
        real(kind=8) :: U0, K0tab(10000)
        real(kind=8):: coordcreate(3), dE, Um_Ut, Usurf, Uj
        integer(kind=8)::i_create
        real(kind=8) :: rij, tabul, tauNa, Wc
        real(kind=8):: randomm
        

        call random_number(side)
        if (side>0.5) then
            side = 1.0
        else
            side = -1.0
        end if
        call random_number(xcr) 
        xcr = xcr * lambda * 0.6
        !write(*,*) xch
        if (side>0.) then
            xcr = length / 2. - xcr-0.4*lambda
        else
            xcr = -length / 2. + xcr+0.4*lambda
        end if

        coordcreate(1) = xcr
        call random_number(ycr)
        ycr = ycr * width - width / 2.0
        coordcreate(2) = ycr
        if (H0+HI<0.) then 
            coordcreate(3) = -1.
        else 
            coordcreate(3) = 1.
        end if
        
        dE = epss
        !Новая энергия системы
        do j = 1, Nminus-1       
            call distance1(r1minus(j, 1), coordcreate(1), r1minus(j, 2), coordcreate(2), width, rij)
            !Считаем энергию с помощью табулированного массива
            if(rij < 7.0 * lambda) then 
                i_create = ceiling(rij*tabul/(7.0*lambda)) 
                dE = dE - U0* K0tab(i_create)*coordcreate(3)
                if(rij<4.*xi) then
                    dE = dE+100000.
                end if
            end if
        end do
        do j = 1, Nplus-1       
            call distance1(r1plus(j, 1), coordcreate(1), r1plus(j, 2), coordcreate(2), width, rij)
            !Считаем энергию с помощью табулированного массива
            if(rij < 7.0 * lambda) then 
                i_create = ceiling(rij*tabul/(7.0*lambda)) 
                dE = dE + U0* K0tab(i_create)*coordcreate(3)
                if(rij<4.*xi) then
                    dE = dE+100000.
                end if
            end if
        end do

        call Ucurrent(coordcreate(1), coordcreate(1), H0, lambda, length, HI, Um_Ut, coordcreate(3))
        dE = dE + Um_Ut*delta*10.**(-7.)
        if(coordcreate(3)>0.) then
            r1plus(Nplus,:) = coordcreate(1:2)
            call Usurface2(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
        Nplus, Nminus, width, length, coordcreate(3), Nplus)
        else
            r1minus(Nminus,:) = coordcreate(1:2)
            call Usurface2(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
        Nplus, Nminus, width, length, coordcreate(3), Nminus)
        end if
        
        dE = dE + Usurf
        Wc = tauNa / float(Nplus+Nminus-2 + 1) * exp((-1.0)*(1.0 / T) * dE)
        
        call random_number(randomm)
        if (Nplus+Nminus<9999) then
            if (Wc > randomm) then
                Uj = Uj + dE
                !Определяем, вихрь или антивихрь рождается
                if(coordcreate(3)>0.) then
                    Nplus = Nplus + 1
                else
                    Nminus = Nminus + 1 
                end if
            ! write(*,*)'Nplus=', Nplus-1, 'Nminus=', Nminus-1 
            end if
        end if
        r1plus(Nplus,:)=100000.
        r1minus(Nminus, :) =100000.
    end subroutine


    !***********************
    !УНИЧТОЖЕНИЕ
    !***********************

    subroutine annihilation(Uj, lambda, tauNa, T, K0tab, width, length, epss,&
     H0, HI, U0, delta, tabul, r1plus, r1minus, Nplus, Nminus)
        real(kind = 8):: signum_anihil, random_anihil, coordanihil(3)
        integer :: number, Nplus, Nminus, i_an, j

        real(kind=8):: length, width, lambda, delta, T, K0tab(10000)
        real(kind=8) :: H0, HI, tabul
        real(kind=8)::epss,dE, Usurf, Um_Ut,U0, Uj, tauNa
        real(kind=8):: r1plus(10000, 2), r1minus(10000, 2)

        real(kind=8):: Wa, flag

        flag = -20.
        call random_number(signum_anihil)
        if(signum_anihil>0.5) then
            call random_number(random_anihil)
            number = nint(random_anihil * (Nplus-1))+1
            if (r1plus(number, 1)<10000.) then
            coordanihil(1:2) = r1plus(number, :)
            coordanihil(3) = 1.
            flag = 1.
            end if
        else
            call random_number(random_anihil)
            number = nint(random_anihil * (Nminus-1))+1
            if (r1minus(number, 1)<10000.) then
            coordanihil(1:2) = r1minus(number, :)
            coordanihil(3) = -1.
            flag = 1.
            end if
        end if
        if (flag>0.) then
        if (length/ 2. - abs(coordanihil(1)) < lambda .and. length/ 2. - abs(coordanihil(1))>0.4*lambda) then 
            dE = 0.0+epss
            if (coordanihil(3) >0.) then
                do j = 1, Nplus-1       
                if( j.NE. number) then
                    call distance1(r1plus(j, 1), coordanihil(1), r1plus(j, 2), coordanihil(2), width, rij)
                    if(rij < 7.0 * lambda) then 
                        i_an = ceiling(rij*tabul/(7.0*lambda)) 
                        dE = dE + U0 * K0tab(i_an)
                    end if
                end if
                end do
                do j = 1, Nminus-1       
                        call distance1(r1minus(j, 1), coordanihil(1), r1minus(j, 2), coordanihil(2), width, rij)
                    if(rij < 7.0 * lambda) then 
                        i_an = ceiling(rij*tabul/(7.0*lambda)) 
                        dE = dE - U0 * K0tab(i_an)
                    end if
                end do
            else
                do j = 1, Nplus-1       
                        call distance1(r1plus(j, 1), coordanihil(1), r1plus(j, 2), coordanihil(2), width, rij)
                    if(rij < 7.0 * lambda) then 
                        i_an = ceiling(rij*tabul/(7.0*lambda)) 
                        dE = dE - U0 * K0tab(i_an)
                    end if
                end do
                do j = 1, Nminus-1
                if( j.NE. number ) then       
                        call distance1(r1minus(j, 1), coordanihil(1), r1minus(j, 2), coordanihil(2), width, rij)
                    if(rij < 7.0 * lambda) then 
                        i_an = ceiling(rij*tabul/(7.0*lambda)) 
                        dE = dE + U0 * K0tab(i_an)
                    end if
                end if
                end do
            end if
            
            call Ucurrent(coordanihil(1), coordanihil(1), H0, lambda, length, HI, Um_Ut, coordanihil(3))
            dE = dE + Um_Ut*delta*10.**(-7.)

            call Usurface(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
                    Nplus, Nminus, width, length, coordanihil(3), number)
            dE = dE + Usurf

            Wa = float(Nplus+Nminus-2) / tauNa * exp(-(1.0 / T) * (-dE))
            ! write(*,*) dE, Um_Ut*delta*10.**(-7.), Usurf, Wa, exp(-(1.0 / T) * (-dE))
            call random_number(random_anihil)
            if (Wa > random_anihil) then
                if(coordanihil(3)>0) then
                    r1plus(number, :) = r1plus(Nplus-1, :)
                    r1plus(Nplus-1, :) = 100000. 
                    Nplus = Nplus - 1
                else
                    r1minus(number, :) = r1minus(Nminus-1, :)
                    r1minus(Nminus-1, :) = 100000. 
                    Nminus = Nminus - 1
                end if
                Uj = Uj - dE
            end if
        end if
        end if
    end subroutine

    !****************
    !Движение
    !***************
    subroutine movement(Uj, lambda, T, length, width, epss, tabul, H0,&
                 HI, U0, K0tab, delta, r1minus, r1plus, Nminus, Nplus, xi)
        real(kind=8) :: dE, U0, K0tab(10000), Um_Ut, Usurf, Uj, epss
        real(kind=8) :: length, width, delta, lambda, tabul, xi
        real(kind=8) :: pi, H0, HI, T, possib1

        real(kind=8) :: r1minus(10000, 2), r1plus(10000, 2)
        integer :: Nplus, Nminus

        real(kind=8) :: coordnew(3), signum_move, temp, coordold(2)
        integer :: number_move, j, i_move 

        real(kind=8) :: randd, poschange, angle, xchange, ychange, num_temp, flag
        dE = 0.0
        pi = acos(-1.0)
        call random_number(signum_move)
        flag=-20.
        if(signum_move>0.5) then
            call random_number(num_temp)  !Выбираем случайный вихрь
            number_move = ceiling(num_temp * (Nplus-1))
            if(number_move<1) then
                number_move = 1
            end if
            if(r1plus(number_move,1)<10000.) then           
                coordnew(1) = r1plus(number_move, 1)
                coordnew(2) = r1plus(number_move, 2)
                coordnew(3) = 1.
                flag=1.
            end if
        else
            call random_number(num_temp)  !Выбираем случайный вихрь
            number_move = ceiling(num_temp * (Nminus-1)) 
            if(number_move<1) then
                number_move = 1
            end if
            if (r1minus(number_move, 1)< 10000.) then           
                coordnew(1) = r1minus(number_move, 1)
                coordnew(2) = r1minus(number_move, 2)
                coordnew(3) = -1.
                flag=1.
            end if
        end if
        if (flag>0.) then
        if(number_move<1) then
        write(*,*) 'bibka'
        end if
        !Меняем координаты вихря, используя полярные координаты
        call random_number(randd)
        poschange = randd * lambda / 10.0
        call random_number(randd)
        angle = 2 * pi * randd
        xchange = poschange * cos(angle)
        ychange = poschange * sin(angle)

        !Меняем координату х выбранного случайного вихря
        if (abs(coordnew(1) + xchange) < length/2.) then
            coordnew(1) = coordnew(1) + xchange
            !Меняем координату y выбранного случайного вихря
            if(coordnew(2) + ychange < -width / 2.0) then
                temp = coordnew(2) + ychange + width
                coordnew(2) = temp
            else if (coordnew(2) + ychange > width / 2.0) then
                temp = coordnew(2) + ychange - width
                coordnew(2) = temp
            else 
                temp = coordnew(2) + ychange
                coordnew(2) = temp
            end if
            !Считаем энергию переставленного вихря
            if(coordnew(3)>0.) then
                do j = 1, Nplus-1     
                    if( j.NE. number_move) then
                        call distance1(r1plus(j, 1), coordnew(1), r1plus(j, 2),&
                                 coordnew(2), width, rij)
                        if(rij < 7.0 * lambda) then 
                            i_move = ceiling(rij*tabul/(7.0*lambda)) 
                            dE = dE + U0*K0tab(i_move)
                        end if
                    end if
                end do
                do j = 1, Nminus-1       
                        call distance1(r1minus(j, 1), coordnew(1), r1minus(j, 2), &
                         coordnew(2), width, rij)
                        if(rij < 7.0 * lambda) then 
                            i_move = ceiling(rij*tabul/(7.0*lambda)) 
                            dE = dE - U0*K0tab(i_move)
                            if(rij<4.*xi) then
                                dE = dE+100000.
                            end if
                        end if
                end do
            else
                do j = 1, Nplus-1     
                        call distance1(r1plus(j, 1), coordnew(1), r1plus(j, 2),&
                         coordnew(2), width, rij)
                        if(rij < 7.0 * lambda) then 
                            i_move = ceiling(rij*tabul/(7.0*lambda)) 
                            dE = dE - U0*K0tab(i_move)
                            if(rij<4.*xi) then
                                dE = dE+100000.
                            end if
                        end if
                end do
                do j = 1, Nminus-1       
                    if( j.NE. number_move) then
                        call distance1(r1minus(j, 1), coordnew(1), r1minus(j, 2),&
                         coordnew(2), width, rij)
                        !вычисляем энергию при помощи табулированного массива
                        if(rij < 7.0 * lambda) then 
                            i_move = ceiling(rij*tabul/(7.0*lambda)) 
                            dE = dE + U0*K0tab(i_move)
                        end if
                    end if
                end do
            end if
            



            ! do j = 1, Nplus-1     
            !     if( j.NE. number_move .or. coordnew(3) < 0.) then
            !         call distance1(r1plus(j, 1), coordnew(1), r1plus(j, 2), coordnew(2), width, rij)
            !         !Считаем энергию с помощью табулированного массива
            !         if(rij < 7.0 * lambda) then 
            !             i_move = ceiling(rij*tabul/(7.0*lambda)) 
            !             dE = dE + U0*K0tab(i_move)*coordnew(3)
            !         end if
            !     end if
            ! end do
            ! do j = 1, Nminus-1       
            !     if( j.NE. number_move .or. coordnew(3) > 0.) then
            !         call distance1(r1minus(j, 1), coordnew(1), r1minus(j, 2), coordnew(2), width, rij)
            !         !вычисляем энергию при помощи табулированного массива
            !         if(rij < 7.0 * lambda) then 
            !             i_move = ceiling(rij*tabul/(7.0*lambda)) 
            !             dE = dE - U0*K0tab(i_move)*coordnew(3)
            !         end if
            !     end if
            ! end do



            !Энергия взаимодействия старого вихря
            if(coordnew(3)>0.) then
                do j = 1, Nplus-1     
                    if( j.NE. number_move) then
                        call distance1(r1plus(j, 1), r1plus(number_move,1), r1plus(j, 2),&
                                 r1plus(number_move, 2), width, rij)
                        if(rij < 7.0 * lambda) then 
                            i_move = ceiling(rij*tabul/(7.0*lambda)) 
                            dE = dE - U0*K0tab(i_move)
                        end if
                    end if
                end do
                do j = 1, Nminus-1       
                        call distance1(r1minus(j, 1), r1plus(number_move,1), r1minus(j, 2), &
                         r1plus(number_move,2), width, rij)
                        if(rij < 7.0 * lambda) then 
                            i_move = ceiling(rij*tabul/(7.0*lambda)) 
                            dE = dE + U0*K0tab(i_move)
                        end if
                end do
            else
                do j = 1, Nplus-1     
                        call distance1(r1plus(j, 1), r1minus(number_move,1), r1plus(j, 2),&
                         r1minus(number_move,2), width, rij)
                        if(rij < 7.0 * lambda) then 
                            i_move = ceiling(rij*tabul/(7.0*lambda)) 
                            dE = dE + U0*K0tab(i_move)
                        end if
                end do
                do j = 1, Nminus-1       
                    if( j.NE. number_move) then
                        call distance1(r1minus(j, 1), r1minus(number_move,1), r1minus(j, 2),&
                         r1minus(number_move,2), width, rij)
                        !вычисляем энергию при помощи табулированного массива
                        if(rij < 7.0 * lambda) then 
                            i_move = ceiling(rij*tabul/(7.0*lambda)) 
                            dE = dE - U0*K0tab(i_move)
                        end if
                    end if
                end do
            end if
            

            call Ucurrent(coordnew(1), coordnew(1), H0, lambda, length, HI, Um_Ut, coordnew(3))
            dE = dE + Um_Ut*delta*10.**(-7.)
            
            if(coordnew(3)>0.) then
                call Ucurrent(r1plus(number_move, 1), r1plus(number_move, 1), H0, lambda, length, HI, Um_Ut, coordnew(3))
                dE = dE - Um_Ut*delta*10.**(-7.)
            else
                call Ucurrent(r1minus(number_move, 1), r1minus(number_move, 1), H0, lambda, length, HI, Um_Ut, coordnew(3))
                dE = dE - Um_Ut*delta*10.**(-7.)
            end if

            
            call Usurface(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
                Nplus, Nminus, width, length, coordnew(3), number_move)
            dE = dE - Usurf

            
            !Запишем старые координаты в переменную, новую в массив, чтобы посчитать энергию с границей
            if(coordnew(3)>0.) then
                coordold = r1plus(number_move,:)
                r1plus(number_move,1) =coordnew(1)
                r1plus(number_move,2) =coordnew(2)
            else
                coordold = r1minus(number_move,:)
                r1minus(number_move,:) =coordnew(:2)
            end if
            call Usurface(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
                Nplus, Nminus, width, length, coordnew(3), number_move)
            dE = dE + Usurf

            if (dE<0) then
                    Uj = Uj + dE    
            else
                possib1 = 1.0 / e ** ( dE / T)
                call random_number(randd)
                if(possib1>randd) then
                    Uj = Uj + dE
                    !write(*,*) 'bibi' 
                    
                else
                if(coordnew(3)>0.) then
                    r1plus(number_move,:) = coordold(:2)
                else
                    r1minus(number_move,:) = coordold(:2)
                end if 
                end if
            end if
        !if((i > iter - 11)) then
        !    do i2 = 1, N
        !        a = ceiling(r1(i2,1) / 2.5 * 200.0)
        !        b = ceiling(r1(i2,2) / 2.5 * 200.0)
        
        !        sectora(a,b) = sectora(a, b) + 1
        !    end do
        !end if   
        end if
        else
            if(signum_move>0. .and. number_move>1) then
            write(*, *) number_move, Nplus
            end if
        end if

    end subroutine

    
    !******************
    !Взаимодействие с границей
    !******************
    subroutine Usurface(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
     Nplus, Nminus, width, length, signum_, number)

        real(kind=8):: U0, K0tab(10000), r1plus(10000,2), r1minus(10000,2), ri(3),lambda,tabul
        real(kind=8):: xvort, xdiff, ydiff, ximage, width, length, rij, Usurf
        integer :: Nplus, Nminus, i5, i6, i7, j, number
        real(kind=8):: signum_
        Usurf =0.
        if(signum_>0.) then
            ri(1:2) = r1plus(number, :)
            ri(3) = 1.
        else if (signum_<0.) then
            ri(1:2) = r1minus(number, :)
            ri(3) = -1.
        end if

        xvort = length / 2.0 - abs(ri(1))
        if (2.* xvort < 7.*lambda) then
            i5 = ceiling(2.0*xvort/(7.0*lambda)*tabul)
            Usurf = Usurf - 1.0/2.0*U0*(K0tab(i5)) 
        end if
        if (2.*(length-xvort) < 7.*lambda) then 
            i6 = ceiling(2.0*(length-xvort)/(7.0*lambda)*tabul)
            Usurf = Usurf - 1./2.*U0*K0tab(i6)
        end if
        
        if(signum_>0.) then
            do j =1, Nminus-1
                !r1(j,3)- знак вихря, далее к координате прибавляем 2 расстояния если справа, отнимаеам два расст если слева
                if (r1minus(j,1)>0) then
                    ximage = r1minus(j, 1) + 2.0*(length/2.0-r1minus(j, 1))
                else
                    ximage = r1minus(j, 1) - 2.0*(length/2.0-abs(r1minus(j, 1)))
                end if
                xdiff = abs(ri(1) - ximage)
                ydiff = min(abs(r1minus(j,2)-ri(2)), width - abs(r1minus(j,2) - ri(2)))
                rij = sqrt(xdiff**2.+ydiff**2.)
                if (rij<7.0*lambda) then
                    i7 = ceiling(rij*tabul/(7.0*lambda))
                    Usurf = Usurf + U0 * K0tab(i7) 
                end if
            end do
            do j =1, Nplus-1
            !r1(j,3)- знак вихря, далее к координате прибавляем 2 расстояния если справа, отнимаеам два расст если слева
                if (r1plus(j,1)>0) then
                    ximage = r1plus(j, 1) + 2.0*(length/2.0-r1plus(j, 1))
                else
                    ximage = r1plus(j, 1) - 2.0*(length/2.0-abs(r1plus(j, 1)))
                end if
                xdiff = abs(ri(1) - ximage)
                ydiff = min(abs(r1plus(j,2)-ri(2)), width - abs(r1plus(j,2) - ri(2)))
                rij = sqrt(xdiff**2.+ydiff**2.)
                if (rij<7.0*lambda) then
                    i7 = ceiling(rij*tabul/(7.0*lambda))
                    Usurf = Usurf -  U0 * K0tab(i7)
                end if
            end do
        else
            do j =1, Nminus-1
                !r1(j,3)- знак вихря, далее к координате прибавляем 2 расстояния если справа, отнимаеам два расст если слева
                if (r1minus(j,1)>0) then
                    ximage = r1minus(j, 1) + 2.0*(length/2.0-r1minus(j, 1))
                else
                    ximage = r1minus(j, 1) - 2.0*(length/2.0-abs(r1minus(j, 1)))
                end if
                xdiff = abs(ri(1) - ximage)
                ydiff = min(abs(r1minus(j,2)-ri(2)), width - abs(r1minus(j,2) - ri(2)))
                rij = sqrt(xdiff**2.+ydiff**2.)
                if (rij<7.0*lambda) then
                    i7 = ceiling(rij*tabul/(7.0*lambda))
                    Usurf = Usurf - U0 * K0tab(i7)
                end if
            end do
            do j =1, Nplus-1
            !r1(j,3)- знак вихря, далее к координате прибавляем 2 расстояния если справа, отнимаеам два расст если слева
                if (r1plus(j,1)>0.) then
                    ximage = r1plus(j, 1) + 2.0*(length/2.0-r1plus(j, 1))
                else
                    ximage = r1plus(j, 1) - 2.0*(length/2.0-abs(r1plus(j, 1)))
                end if
                xdiff = abs(ri(1) - ximage)
                ydiff = min(abs(r1plus(j,2)-ri(2)), width - abs(r1plus(j,2) - ri(2)))
                rij = sqrt(xdiff**2.+ydiff**2.)
                if (rij<7.0*lambda) then
                    i7 = ceiling(rij*tabul/(7.0*lambda))
                    Usurf = Usurf + U0 * K0tab(i7)
                end if
            end do
        end if
    
    end subroutine

    subroutine Usurface2(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
     Nplus, Nminus, width, length, signum_, number)

        real(kind=8):: U0, K0tab(10000), r1plus(10000,2), r1minus(10000,2), ri(3),lambda,tabul
        real(kind=8):: xvort, xdiff, ydiff, ximage, width, length, rij, Usurf
        integer :: Nplus, Nminus, i5, i6, i7, j, number
        real(kind=8):: signum_, coord(3)
        Usurf =0.
        if(signum_>0.) then
            ri(1:2) = r1plus(number, :)
            ri(3) = 1.
        else if (signum_<0.) then
            ri(1:2) = r1minus(number, :)
            ri(3) = -1.
        end if

        xvort = length / 2.0 - abs(ri(1))
        if (2.* xvort < 7.*lambda) then
            i5 = ceiling(2.0*xvort/(7.0*lambda)*tabul)
            Usurf = Usurf - 1.0/2.0*U0*(K0tab(i5)) 
        end if
        if (2.*(length-xvort) < 7.*lambda) then 
            i6 = ceiling(2.0*(length-xvort)/(7.0*lambda)*tabul)
            Usurf = Usurf - 1./2.*U0*K0tab(i6)
        end if

        if(signum_>0.) then
            do j =1, Nminus-1
                !r1(j,3)- знак вихря, далее к координате прибавляем 2 расстояния если справа, отнимаеам два расст если слева
                if (r1minus(j,1)>0) then
                    ximage = r1minus(j, 1) + 2.0*(length/2.0-r1minus(j, 1))
                else
                    ximage = r1minus(j, 1) - 2.0*(length/2.0-abs(r1minus(j, 1)))
                end if
                xdiff = abs(ri(1) - ximage)
                ydiff = min(abs(r1minus(j,2)-ri(2)), width - abs(r1minus(j,2) - ri(2)))
                rij = sqrt(xdiff**2.+ydiff**2.)
                if (rij<7.0*lambda) then
                    i7 = ceiling(rij*tabul/(7.0*lambda))
                    Usurf = Usurf + U0 * K0tab(i7) 
                end if
            end do
            do j =1, Nplus
            !r1(j,3)- знак вихря, далее к координате прибавляем 2 расстояния если справа, отнимаеам два расст если слева
                if (r1plus(j,1)>0) then
                    ximage = r1plus(j, 1) + 2.0*(length/2.0-r1plus(j, 1))
                else
                    ximage = r1plus(j, 1) - 2.0*(length/2.0-abs(r1plus(j, 1)))
                end if
                xdiff = abs(ri(1) - ximage)
                ydiff = min(abs(r1plus(j,2)-ri(2)), width - abs(r1plus(j,2) - ri(2)))
                rij = sqrt(xdiff**2.+ydiff**2.)
                if (rij<7.0*lambda) then
                    i7 = ceiling(rij*tabul/(7.0*lambda))
                    Usurf = Usurf - U0 * K0tab(i7)
                end if
            end do
        else
            do j =1, Nminus
                !r1(j,3)- знак вихря, далее к координате прибавляем 2 расстояния если справа, отнимаеам два расст если слева
                if (r1minus(j,1)>0) then
                    ximage = r1minus(j, 1) + 2.0*(length/2.0-r1minus(j, 1))
                else
                    ximage = r1minus(j, 1) - 2.0*(length/2.0-abs(r1minus(j, 1)))
                end if
                xdiff = abs(ri(1) - ximage)
                ydiff = min(abs(r1minus(j,2)-ri(2)), width - abs(r1minus(j,2) - ri(2)))
                rij = sqrt(xdiff**2.+ydiff**2.)
                if (rij<7.0*lambda) then
                    i7 = ceiling(rij*tabul/(7.0*lambda))
                    Usurf = Usurf - U0 * K0tab(i7)
                end if
            end do
            do j =1, Nplus-1
            !r1(j,3)- знак вихря, далее к координате прибавляем 2 расстояния если справа, отнимаеам два расст если слева
                if (r1plus(j,1)>0.) then
                    ximage = r1plus(j, 1) + 2.0*(length/2.0-r1plus(j, 1))
                else
                    ximage = r1plus(j, 1) - 2.0*(length/2.0-abs(r1plus(j, 1)))
                end if
                xdiff = abs(ri(1) - ximage)
                ydiff = min(abs(r1plus(j,2)-ri(2)), width - abs(r1plus(j,2) - ri(2)))
                rij = sqrt(xdiff**2.+ydiff**2.)
                if (rij<7.0*lambda) then
                    i7 = ceiling(rij*tabul/(7.0*lambda))
                    Usurf = Usurf + U0 * K0tab(i7)
                end if
            end do
        end if
    
    end subroutine
    !******************
    !Взаимодействие с границей
    !******************
    subroutine Usurface1(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
     Nplus, Nminus, width, length, signum_, number)

        real(kind=8):: U0, K0tab(10000), r1plus(10000,2), r1minus(10000,2), ri(3),lambda,tabul
        real(kind=8):: xvort, xdiff, ydiff, ximage, width, length, rij, Usurf
        integer :: Nplus, Nminus, i5, i6, i7, j, number
        real(kind=8):: signum_, flag
        flag = -1.
        Usurf=0.
        if(signum_>0.) then
            ri(1:2) = r1plus(number, :)
            ri(3) = 1.
            flag = 1.
        else if (signum_<0.) then
            ri(1:2) = r1minus(number, :)
            ri(3) = -1.
            flag = 1.
        end if
        if(flag>0.) then

        xvort = length / 2.0 - abs(ri(1))
        
        if (2.* xvort < 7.*lambda) then
            i5 = ceiling(2.0*xvort/(7.0*lambda)*tabul)
            Usurf = Usurf - 1.0/2.0*U0*(K0tab(i5)) 
        end if
        if (2.*(length-xvort) < 7.*lambda) then 
            i6 = ceiling(2.0*(length-xvort)/(7.0*lambda)*tabul)
            Usurf = Usurf - 1./2.*U0*K0tab(i6)
        end if
        end if
    end subroutine
    !******************
    !Аннигиляция пары
    !******************
    subroutine pair_annihilation(lambda, xi, T, tauNa, delta, width, length, tabul, &
    epss,r1plus,r1minus,Nplus,Nminus, H0,HI, U0, K0tab)
        real(kind=8):: lambda, xi, epss,r1plus(10000,2),r1minus(10000,2), plusminus
        real(kind=8):: width, length, rplus(2), rminus(2), xdiff, ydiff, rij
        real(kind=8):: T, dE, x1plus, x1minus, x2, y1plus, y1minus, y2, Wa
        real(kind=8)::ydiff1,ydiff2,xdiff1,xdiff2, rij1,rij2, ri(3), tauNa
        real(kind=8)::H0,HI, Um_Ut, Usurf, U0, K0tab(10000), delta, pairanihil
        real(kind=8) :: tabul, sign_, rij_plus_minus, ximage
        integer :: Nplus, Nminus, Nplus_, Nminus_, j, i4

        call random_number(plusminus)
        Nplus_ = nint(plusminus * (Nplus-1))+1
        call random_number(plusminus)
        Nminus_ = nint(plusminus * (Nminus-1))+1  
        if (Nplus_<1 .or. Nplus_==Nplus) then
            Nplus_=1
        end if
        if (Nminus_<1 .or. Nminus_==Nminus) then
            Nminus_=1
        end if
        !Задаем координаты нужных вихрей
        rplus = r1plus(Nplus_, :)
        rminus = r1minus(Nminus_, :)
        !Проверка по расстоянию, достаточно ли близко, чтобы аннигилировать
        call distance1(rplus(1), rminus(1), rplus(2), rminus(2), width, rij_plus_minus)
        if (rij_plus_minus<10.*xi) then
            dE= 2. * epss
            !Считаем взаимодействия обоих вихрей сначала с положительными, затем с отрицательными вихрями
            do j = 1, Nplus-1       !
                if( j.NE. Nplus_ ) then

                    call distance1(rplus(1), r1plus(j, 1), rplus(2), r1plus(j, 2), width, rij)
                    if(rij < 7.0 * lambda) then 
                        i4 = ceiling(rij*tabul/(7.0*lambda)) 
                        dE = dE + U0*K0tab(i4)
                    end if

                    call distance1(rminus(1), r1plus(j, 1), rminus(2), r1plus(j, 2), width, rij)
                    if(rij < 7.0 * lambda) then 
                        i4 = ceiling(rij*tabul/(7.0*lambda)) 
                        dE = dE - U0*K0tab(i4)
                    end if
                end if
            end do

            do j = 1, Nminus-1       !
                if( j.NE. Nminus_ ) then
                    call distance1(rplus(1), r1minus(j, 1), rplus(2), r1minus(j, 2), width, rij)
                    if(rij < 7.0 * lambda) then 
                        i4 = ceiling(rij*tabul/(7.0*lambda)) 
                        dE = dE - U0*K0tab(i4)
                    end if

                    call distance1(rminus(1), r1minus(j, 1), rminus(2), r1minus(j, 2), width, rij)
                    if(rij < 7.0 * lambda) then 
                        i4 = ceiling(rij*tabul/(7.0*lambda)) 
                        dE = dE + U0*K0tab(i4)
                    end if
                end if
            end do
            !Их взаимодействие друг с другом
            i4 = ceiling(rij_plus_minus*tabul/(7.0*lambda))
            dE = dE - U0*K0tab(i4)
            !write(*,*) '++inter', dE
            sign_ = -1.
            call Ucurrent(rminus(1), rminus(1), H0, lambda, length, HI, Um_Ut, sign_)
            dE = dE + Um_Ut*delta*10.**(-7.)
            !write(*,*) '+Uc1', dE
            sign_ = 1.
            call Ucurrent(rplus(1), rplus(1), H0, lambda, length, HI, Um_Ut, sign_)
            dE = dE + Um_Ut*delta*10.**(-7.)
            !write(*,*) '+Uc2', dE
            ri(:2) = rplus
            ri(3) = 1.
            call Usurface(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
        Nplus, Nminus, width, length, ri(3), Nplus_)
            dE = dE + Usurf
            !write(*,*) '+Usurf', dE
            ri(:2) = rminus
            ri(3) = -1.
            call Usurface(U0, K0tab, Usurf, lambda, tabul, r1minus, r1plus,&
        Nplus, Nminus, width, length, ri(3), Nminus_)
            dE = dE + Usurf
            !Два раза учли их вазимодействие с изображениями друг друга
            if (rplus(1)>0.) then
                ximage = rplus(1) + 2.0*(length/2.0-rplus(1))
            else
                ximage = rplus(1) - 2.0*(length/2.0-abs(rplus(1)))
            end if
            xdiff = abs(rminus(1) - ximage)
            ydiff = min(abs(rplus(2)-rminus(2)), width - abs(rplus(2) - rminus(2)))
            rij = sqrt(xdiff**2.+ydiff**2.)
            if (rij<7.0*lambda) then
                i4 = ceiling(rij*tabul/(7.0*lambda))
                dE = dE - 1./2.*U0 * K0tab(i4)
            end if
            !write(*,*) '+Usurf2', dE
            ! write(*,*)  T, dE
            Wa = float(Nplus-1)*float(Nminus-1) /(tauNa * 100.*xi**(2.)*10**(-14.)*&
            length * 10**(-7.)*width*10**(-7.)) * exp(-(1.0 / T) * (-dE))
            ! write(*,*) Wa
            ! read(*,*)
                call random_number(pairanihil)
                if (Wa > pairanihil) then
                    r1plus(Nplus_, :) = r1plus(Nplus-1, :)
                    r1plus(Nplus-1, :) = 100000.
                    Nplus = Nplus-1

                    r1minus(Nminus_, :) = r1minus(Nminus-1, :)
                    r1minus(Nminus-1, :) = 100000.
                    Nminus = Nminus-1

                    Uj = Uj - dE
                end if
        end if
    end subroutine

    subroutine besssel ( x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1 )
            
        !*****************************************************************************80
        !
        !! IK01A compute Bessel function I0(x), I1(x), K0(x), and K1(x).
        !
        !  Discussion:
        !
        !    This procedure computes modified Bessel functions I0(x), I1(x),
        !    K0(x) and K1(x), and their derivatives.
        !
        !  Licensing:
        !
        !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
        !    they give permission to incorporate this routine into a user program 
        !    provided that the copyright is acknowledged.
        !
        !  Modified:
        !
        !    16 July 2012
        !
        !  Author:
        !
        !    Shanjie Zhang, Jianming Jin
        !
        !  Reference:
        !
        !    Shanjie Zhang, Jianming Jin,
        !    Computation of Special Functions,
        !    Wiley, 1996,
        !    ISBN: 0-471-11963-6,
        !    LC: QA351.C45.
        !
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) X, the argument.
        !
        !    Output, real ( kind = 8 ) BI0, DI0, BI1, DI1, BK0, DK0, BK1, DK1, the
        !    values of I0(x), I0'(x), I1(x), I1'(x), K0(x), K0'(x), K1(x), K1'(x).
        !
        
        real ( kind = 8 ), save, dimension ( 12 ) :: a = (/ &
            0.125D+00, 7.03125D-02, &
            7.32421875D-02, 1.1215209960938D-01, &
            2.2710800170898D-01, 5.7250142097473D-01, &
            1.7277275025845D+00, 6.0740420012735D+00, &
            2.4380529699556D+01, 1.1001714026925D+02, &
            5.5133589612202D+02, 3.0380905109224D+03 /)
        real ( kind = 8 ), save, dimension ( 8 ) :: a1 = (/ &
            0.125D+00, 0.2109375D+00, &
            1.0986328125D+00, 1.1775970458984D+01, &
            2.1461706161499D+02, 5.9511522710323D+03, &
            2.3347645606175D+05, 1.2312234987631D+07 /)
        real ( kind = 8 ), save, dimension ( 12 ) :: b = (/ &
            -0.375D+00, -1.171875D-01, &
            -1.025390625D-01, -1.4419555664063D-01, &
            -2.7757644653320D-01, -6.7659258842468D-01, &
            -1.9935317337513D+00, -6.8839142681099D+00, &
            -2.7248827311269D+01, -1.2159789187654D+02, &
            -6.0384407670507D+02, -3.3022722944809D+03 /)
        real ( kind = 8 ) bi0
        real ( kind = 8 ) bi1
        real ( kind = 8 ) bk0
        real ( kind = 8 ) bk1
        real ( kind = 8 ) ca
        real ( kind = 8 ) cb
        real ( kind = 8 ) ct
        real ( kind = 8 ) di0
        real ( kind = 8 ) di1
        real ( kind = 8 ) dk0
        real ( kind = 8 ) dk1
        real ( kind = 8 ) el
        integer ( kind = 4 ) k
        integer ( kind = 4 ) k0
        real ( kind = 8 ) pi
        real ( kind = 8 ) r
        real ( kind = 8 ) w0
        real ( kind = 8 ) ww
        real ( kind = 8 ) x
        real ( kind = 8 ) x2
        real ( kind = 8 ) xr
        real ( kind = 8 ) xr2
        pi = 3.141592653589793D+00
            el = 0.5772156649015329D+00
            x2 = x * x
        
            if ( x == 0.0D+00 ) then
        
                bi0 = 1.0D+00
                bi1 = 0.0D+00
                bk0 = 1.0D+300
                bk1 = 1.0D+300
                di0 = 0.0D+00
                di1 = 0.5D+00
                dk0 = -1.0D+300
                dk1 = -1.0D+300
                return
        
            else if ( x <= 18.0D+00 ) then
        
                bi0 = 1.0D+00
                r = 1.0D+00
                do k = 1, 50
                r = 0.25D+00 * r * x2 / ( k * k )
                bi0 = bi0 + r
                if ( abs ( r / bi0 ) < 1.0D-15 ) then
                    exit
                end if
                end do
        
                bi1 = 1.0D+00
                r = 1.0D+00
                do k = 1, 50
                    r = 0.25D+00 * r * x2 / ( k * ( k + 1 ) )
                    bi1 = bi1 + r
                    if ( abs ( r / bi1 ) < 1.0D-15 ) then
                        exit
                    end if
                end do
        
            bi1 = 0.5D+00 * x * bi1
        
        else
        
            if ( x < 35.0D+00 ) then
            k0 = 12
            else if ( x < 50.0D+00 ) then
            k0 = 9
            else
            k0 = 7
            end if
        
            ca = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
            bi0 = 1.0D+00
            xr = 1.0D+00 / x
            do k = 1, k0
            bi0 = bi0 + a(k) * xr ** k
            end do
            bi0 = ca * bi0
            bi1 = 1.0D+00
            do k = 1, k0
            bi1 = bi1 + b(k) * xr ** k
            end do
            bi1 = ca * bi1
        
        end if
        
        if ( x <= 9.0D+00 ) then
        
            ct = - ( log ( x / 2.0D+00 ) + el )
            bk0 = 0.0D+00
            w0 = 0.0D+00
            r = 1.0D+00
            do k = 1, 50
            w0 = w0 + 1.0D+00 / k
            r = 0.25D+00 * r / ( k * k ) * x2
            bk0 = bk0 + r * ( w0 + ct )
            if ( abs ( ( bk0 - ww ) / bk0 ) < 1.0D-15 ) then
                exit
            end if
            ww = bk0
            end do
        
            bk0 = bk0 + ct
        
        else
        
            cb = 0.5D+00 / x
            xr2 = 1.0D+00 / x2
            bk0 = 1.0D+00
            do k = 1, 8
            bk0 = bk0 + a1(k) * xr2 ** k
            end do
            bk0 = cb * bk0 / bi0
        
        end if
        
        bk1 = ( 1.0D+00 / x - bi1 * bk0 ) / bi0
        di0 = bi1
        di1 = bi0 - bi1 / x
        dk0 = - bk1
        dk1 = - bk0 - bk1 / x
        
        return
        end subroutine
end program