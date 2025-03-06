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





















program task11

implicit none 
    
    !Все вихри одного знака взял
    real ( kind = 8 ) bi0
    real ( kind = 8 ) bi1
    real ( kind = 8 ) bk0
    real ( kind = 8 ) bk1
    
    real ( kind = 8 ) di0
    real ( kind = 8 ) di1
    real ( kind = 8 ) dk0
    real ( kind = 8 ) dk1
    real(kind = 8), parameter :: k1 = 1.38 * 10.0**(-23)
    real(kind = 8), parameter :: e =  2.71828182845904
    real(kind = 8), parameter :: phi0 =  2.06783460999999E-07            !Квант магнитного потока
    real(kind = 8), parameter :: delta = 0.27       !Толщина слоя
    real(kind = 8) :: xi0, xi          !Длина когерентности сверхпроводника
    real(kind = 8), parameter :: lambda0 = 180.0         !Лондоновская глубина проникновения данного материала при Т=0 
    real(kind = 8) :: lambda             !Лондоновская глубина проникновения данного материала при данной температуре
    real(kind = 8) :: T, Tc
    real(kind = 8) :: length, width
    integer :: N = 1, i, j, k, i4, i1, iter, a, b, nomera(-200:200), i2
    real(kind = 8) :: temp, temp1, randd, step
    real(kind = 8) :: Ui, Uj, U0, Ujone, Ujj, Unewone, dE, dE1, Uh !Энергии взаимодействия вихрей
    real(kind = 8) :: poschange, angle, xchange, ychange !на сколько изменится положение вихря
    real(kind = 8) :: possib1 !Вероятность перехода
    real(kind = 8) :: rij !Расстояние между вихрями
    real(kind = 8) :: x, x1, x2, y1, y2, ydiff, xdiff !Координаты вихрей
    real(kind = 8):: r1(10000, 2)!, rmed(1000000,2), countmed
    real(kind = 8) :: pi
    real(kind = 8) :: epss  ! Собственная энергия всех вихрей
    real(kind = 8):: Utab(10000), sectora(-200:200,-200:200), coordnew(2), r3(2,2)
    real(kind = 8) :: number1, tabul, tauNa, H  !будем каждый раз выбирать вихрь по его номеру в массиве
    integer :: number, N1
    real(kind = 8) :: process
    tauNa = 100.0
    T = 1.0
    xi0 = 2.0
    iter = 2000

    tabul = 10000.0
    length = 5.0*10**3 ! по х
    width = 5.0*10**3   !по у
    Tc = 84.0
    H = 400.0
    
    open(2,file = "res1.txt")
    open(3,file = "res2.txt")
    open(4, file = "task1energy.txt")
    xi = xi0 * (1.0 - (T / Tc) ** 3.3) ** (-1.0 / 2.0)
    write(*,*) xi 
    lambda = lambda0 * (1.0 - (T / Tc) ** (3.3)) ** (-1.0 / 2.0)
    write(*,*) lambda
    pi = acos(-1.0)
    !call random_number(r1) 
    T= T / 11606.0 ! T в эВ
    !r1 = (r1 * width - width / 2.0)
    !считаем собственные энергии вихрей
    epss = (phi0 / (4 * pi * lambda * 10 ** (-7.0))) ** 2 * (log(lambda / xi) + 0.52)  * delta * 10 ** (-7.0) *6.2 * 10.0**(11)
    write(*,*) epss
    Uh = (-1.0) * phi0 * H * delta * 10**(-7.0)/ (4.0 * pi) * 6.2 * 10 ** (11.0)
    write(*,*) Uh 


    !Табулирование

    step = 7.0 * lambda / tabul
    do i = 1, floor(tabul)
        x = step * i / lambda
        call besssel(x,bi0, di0, bi1, di1, bk0, dk0, bk1, dk1)
        temp = (phi0**2.0 / (8.0 * pi**2.0 * (lambda * 10 ** (-7.0))**2.0)) * delta * 10**(-7.0)
        Utab(i) = temp * bk0
    end do
    
    Utab = Utab * 6.2 * 10.0 **(11.0)
    Uj = 0.0
    !countmed = 1
    do i = 1, iter
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
                if(process < 0.4) then

                     dE = 0.0
                     call random_number(number1)  !Выбираем
                     !случайный вихрь
                    number = nint(number1 * (N-2))+1            
                    !Меняем координаты вихря, используя полярные координаты
                    coordnew(1) = r1(number, 1)
                    coordnew(2) = r1(number, 2)
                    call random_number(randd)
                    poschange = randd * lambda / 10.0
                    call random_number(randd)
                    angle = 2 * pi * randd
                    xchange = poschange * cos(angle)
                    ychange = poschange * sin(angle)
                     !записали в массив старые значения положения вихрей
                    !Меняем координату х выбранного случайного вихря
                    if(coordnew(1) + xchange < -length / 2.0) then
                        temp = coordnew(1) + xchange + length
                       coordnew(1) = temp
                    else if (coordnew(1) + xchange > length / 2.0) then
                        temp = coordnew(1) + xchange - length
                        coordnew(1) = temp
                    else 
                        temp = coordnew(1) + xchange
                        coordnew(1) = temp
                    end if
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
                    do j = 1, N-1       !
                        if( j.NE. number ) then
                            x1 = r1(j, 1)
                            x2 = coordnew(1)
                            y1 = r1(j, 2)
                            y2 = coordnew(2)
                            ydiff = min(abs(y1-y2), width - abs(y1 - y2))
                            xdiff = min(abs(x1-x2), length - abs(x1 - x2))  
                            rij = sqrt(ydiff**2 + xdiff**2)
                            !Приводим к нужной размерности
                            !Считаем энергию с помощью табулированного массива
                            if(rij < 7.0 * lambda) then 
                                i4 = ceiling(rij*tabul/(7.0*lambda)) 
                                dE = dE + Utab(i4)
                            end if
                        end if
                    end do
                    do j = 1, N-1       
                        if( j.NE. number ) then
                            x1 = r1(number, 1)
                            x2 = r1(j, 1)
                            y1 = r1(number, 2)
                            y2 = r1(j, 2)
                            ydiff = min(abs(y1-y2), width - abs(y1 - y2))
                            xdiff = min(abs(x1-x2), length - abs(x1 - x2))  
                            rij = sqrt(ydiff**2 + xdiff**2)
                            !приводим к нужной размерности 
                            !вычисляем энергию при помощи табулированного массива
                            if(rij < 7.0 * lambda) then 
                                i4 = ceiling(rij*tabul/(7.0*lambda)) 
                                dE = dE - Utab(i4)
                            end if
                        end if
                    end do

                    if (dE<0) then
                        r1(number, 1) = coordnew(1)
                        r1(number, 2) = coordnew(2)
                        Uj = Uj + dE
                
                    else
                        possib1 = 1.0 / e ** ( dE / T)
                        call random_number(randd)
                        if(possib1>randd) then
                            r1(number, 1) = coordnew(1)
                            r1(number, 2) = coordnew(2)
                            Uj = Uj + dE
                        end if
                    end if
                !if((i > iter - 11)) then
                !    do i2 = 1, N
                !        a = ceiling(r1(i2,1) / 2.5 * 200.0)
                !        b = ceiling(r1(i2,2) / 2.5 * 200.0)
                
                !        sectora(a,b) = sectora(a, b) + 1
                !    end do
                !end if   
                else if(process < 0.7) then
                    call creation(Uj, r1, N, lambda, tauNa, T, width, length, Utab, epss,Uh)
                    
                else
                    call annihilation(Uj, r1, N, lambda, tauNa, T, Utab, width, length, epss, Uh)
                end if
            end do 
        
        !вывод жнергии каждые 10 итераций
            if (mod(i,10)==0) then
                write(*,*) i, 'Uj=', Uj
            end if

!            !if (mod(i,10)==0 .and. i>15000) then
!            !    write(77,*) rmed
!            !end if
        !запись энергии в файл
            write(4,*) i, Uj
        
    end do
    sectora = sectora / 1000.0
    !write(*,*) sectora

44    Ujj = 0.0
    write(*,*) N
    do j = 1, N-1     !Считаем
        do k = 1, N-1    !энергию начальной конфигурации
            if( j .NE. k) then
                !Координаты двух вихрей, в массиве r1 100 строк и два стлбца- 100 вихрей и у каждого 2 координаты
                x1 = r1(j, 1)
                x2 = r1(k, 1)
                y1 = r1(j, 2)
                y2 = r1(k, 2)
                !Вычисляем расстояние по х и по у
                ydiff = min(abs(y1-y2), abs(width-max(y1,y2) + min(y1,y2)))
                xdiff = min(abs(x1-x2), abs(length-max(x1,x2) + min(x1,x2)))
                rij = sqrt(ydiff**2.0 + xdiff**2.0)
                !write(*,*) rij
                !Находим энергию из табулированного массива
                if(rij < 7.0 * lambda) then 
                    i4 = ceiling(rij*tabul/(7.0*lambda)) 
            
                    temp = Utab(i4)
                    Ujj = Ujj + temp / 2.0
                
                end if
            end if
        end do
    end do
    
    write(*,*) Ujj +(N-1) * (Uh +epss)
    !write(*,*) r1(1:11,:)
    do i = 1, N-1
        write(2, *) r1(i, :)
    end do
    open(5, file="heatmap_task1_data.txt")
    do i = -200, 200
        nomera(i) = i
    end do
    write(5, *) nomera
    do i = -200, 200
        !write(5, *) i, sectora(i, :)
    end do
    close(5)
    close(2)
    close(3)
    close(4)
    !plot 'C:/Users/Valera P. Lenkov/Desktop/Fortran/task1/res2.txt'











    contains

    !****************************
    !РОЖДЕНИЕ
    !****************************
    subroutine creation(Uj, r1, N, lambda, tauNa, T, width, length, Utab, epss, Uh)
    real(kind=8) :: pc, Wc, dE, Ly, r1(10000,2), tauNa, width, length, Utab(10000), Uj, lambda, T, epss, Uh!Что такое Ly? Ширина по у
    integer :: j, N !lfn - последнее свободное место в массиве вихрей
    real(kind=8) :: ycr, xch, randomm, coordcreate(2), temp1, temp2
    call random_number(xch) 
    !Выбираем право или лево
    !xch = xch * 2 *  lambda - lambda
    !if (xch < 0) then
    !    coordcreate(1) = length / 2.0  + xch
    !else 
    !    coordcreate(1) = -length / 2.0 + xch
    !end if
    !call random_number(xch)
    xch = xch * length - length / 2.0
    coordcreate(1) = xch
    call random_number(ycr)
    ycr = ycr * width - width / 2.0
    coordcreate(2) = ycr
    !write(*,*) coordcreate
    !pc = 1.0 / (lambda * Ly)
    !write(*,*) Uh, epss
    dE = Uh+epss
    !write(*,*) dE, N
    do j = 1, N-1       !Новая энергия системы
                    x1 = r1(j, 1)
                    x2 = coordcreate(1)
                    y1 = r1(j, 2)
                    y2 = coordcreate(2)
                    ydiff = min(abs(y1-y2), width - abs(y1 - y2))
                    xdiff = min(abs(x1-x2), length - abs(x1 - x2))  
                    rij = sqrt(ydiff**2 + xdiff**2)
                    !Приводим к нужной размерности
                    !Считаем энергию с помощью табулированного массива

                    if(rij < 7.0 * lambda) then 
                        i4 = ceiling(rij*tabul/(7.0*lambda)) 
                        dE = dE + Utab(i4)
                    end if
    end do
    !write(*,*) dE
    Wc = tauNa / (N-1.0 + 1.0) * exp((-1.0)*(1.0 / T) * dE)
    call random_number(randomm)
    if (N<9999) then
        if (Wc > randomm) then
            r1(N, 2) = coordcreate(2)
            r1(N, 1) = coordcreate(1)
            N = N + 1
            write(*,*) N
            Uj = Uj + dE
        end if
    end if
end subroutine


!***********************
!УНИЧТОЖЕНИЕ
!***********************

subroutine annihilation(Uj, r1, N, lambda, tauNa, T, Utab, width, length, epss, Uh)
    !А зачем зона уничтожения, если выбираем произвольный вихрь из массива?
    real(kind=8) :: pa, Wa, dE, Uj, r1(10000, 2), lambda, tauNa, Utab(10000), width, length, T, epss, Uh
    real(kind=8) :: x1, x2, rij, y1, ydiff, xdiff
    real(kind=8) :: random_anihil, coordanihil(2)
    integer :: number, j, N
    call random_number(random_anihil)
    number = nint(random_anihil * (N-2))+1
    coordanihil = r1(number, :)
    dE = 0.0+epss+ Uh
    do j = 1, N-1       
        if( j.NE. number ) then
                x1 = r1(j, 1)
                x2 = coordanihil(1)
                y1 = r1(j, 2)
                y2 = coordanihil(2)
                ydiff = min(abs(y1-y2), width - abs(y1 - y2))
                xdiff = min(abs(x1-x2), length - abs(x1 - x2))  
                rij = sqrt(ydiff**2 + xdiff**2)
                !Приводим к нужной размерности
                !Считаем энергию с помощью табулированного массива
            if(rij < 7.0 * lambda) then 
                i4 = ceiling(rij*tabul/(7.0*lambda)) 
                dE = dE + Utab(i4)
            end if
        end if
    end do

    !pa = 1.0 / (2.0 * N)
    
    Wa = (N-1) / tauNa * exp(-(1.0 / T) * (-dE))
    call random_number(random_anihil)
    if (Wa > random_anihil) then
        r1(number, :) = r1(N-1, :)
        r1(N-1, :) = 0
        Uj = Uj - dE
        N = N-1
    end if
end subroutine
end program