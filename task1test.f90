
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














program task1
    
    
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
    real(kind = 8), parameter :: delta = 0.27e-7       !Толщина слоя
    real(kind = 8), parameter :: xi0 = 2e-7           !Длина когерентности сверхпроводника
    real(kind = 8), parameter :: lambda0 = 180e-7          !Лондоновская глубина проникновения данного материала при Т=0 
    real(kind = 8) :: lambda             !Лондоновская глубина проникновения данного материала при данной температуре
    real(kind = 8) :: T, Tc
    integer :: N = 100, i, j, j1, k, i1, i2, i3, i4, sec, sec1, sec10, sec20
    real(kind = 8) :: temp, temp1, randd, step
    real(kind = 8) :: Ui, Uj, U0, Ujone, Uoldone, Unewone !Энергии взаимодействия вихрей
    real(kind = 8) :: poschange, angle, xchange, ychange !на сколько изменится положение вихря
    real(kind = 8) :: possib1 !Вероятность перехода
    real(kind = 8) :: rij !Расстояние между вихрями
    real(kind = 8) :: x, x1, x2, y1, y2, ydiff, xdiff !Координаты вихрей
    real(kind = 8) :: r(100,2)
    real(kind = 8) :: pi
    real(kind = 8) :: epss  ! Собственная энергия всех вихрей
    real(kind = 8):: Ue(100), r1(100,2), r2(100,2), sect(40, 40), Utab(10000)
    real(kind = 8) :: number1 !будем каждый раз выбирать вихрь по его номеру в массиве
    integer :: number
    T = 1.0
    Tc = 84.0
    sect = 0
    Ue = 0
    open(2,file = "res1.txt")
    open(3,file = "res2.txt")
    open(4, file = "task1energy.txt")
    lambda = lambda0 * (1.0 - (T / Tc) ** 3.3) ** (-1.0 / 2.0)
    !write(2,*) lambda
    pi = acos(-1.0)
    call random_number(r) 
    !write(2,*) r * 5, delta, NEW_LINE('a')
    r1 = r * 5
    !считаем собственные энергии вихрей
    epss = (phi0 / (4 * pi * lambda)) ** 2 * (log(lambda0 / xi0) + 0.52)  * delta *6.2 * 10.0**(11)
    write(*,*) epss
     

    !Табулирование

    step = 7.0 * lambda / 10000.0
    do i = 1, 10000
        x = step * i / lambda
        call besssel(x,bi0, di0, bi1, di1, bk0, dk0, bk1, dk1)
        temp = (phi0**2 / (8.0 * pi**2 * lambda**2)) * delta * bk0
        Utab(i) = temp
    end do
    Utab = Utab * 6.2 * 10.0 **(11)









    
    do j = 1, 100         !Считаем
        do k = 1, 100      !энергию начальной конфигурации
            if( j .NE. k) then
                !Координаты двух вихрей, в массиве r1 100 строк и два стлбца- 100 вихрей и у каждого 2 координаты
                x1 = r1(j, 1)
                x2 = r1(k, 1)
                y1 = r1(j, 2)
                y2 = r1(k, 2)
                !Вычисляем расстояние по х и по у
                if(y1>y2) then
                    if(5-y1+y2> y1-y2) then
                        ydiff = y1-y2
                    else
                        ydiff = 5-y1+y2
                    end if
                else
                    if(5-y2+y1> y2-y1) then
                        ydiff = y2-y1
                    else
                        ydiff = 5-y2+y1
                    end if
                end if
                if(x1>x2) then
                    if(5-x1+x2> x1-x2) then
                        xdiff = x1-x2
                    else
                        xdiff = 5-x1+x2
                    end if
                else
                    if(5-x2+x1> x2-x1) then
                        xdiff = x2-x1
                    else
                        xdiff = 5-x2+x1
                    end if
                end if 

                rij = sqrt(ydiff**2 + xdiff**2)
                !Приводим расстояние в СГС
                x = rij*10.0**(-4)
                !Находим энергию из табулированного массива
                if(x < 7 * lambda) then 
                    i4 = ceiling(x*10000.0/(7.0*lambda)) 
            
                    temp = Utab(i4)
                    Ue(j) = Ue(j) + temp
                
                end if
            end if
        end do
    end do
    Uj = sum(Ue) / 2.0
    write(*,*) Uj
    !начальная энергия приблизительно 1.5-2эВ

    !Основной алгоритм
    do i = 1, 5000
        !В каждой итерации N шагов
        do i1 = 1, N
            ! записываем координаты вихрей в еще один массив r2, в нем мы будет переставлять
            ! один из вихрей, после чего, если перестановка удалась, запишем его снова в r1
            r2 = r1
            call random_number(number1)  !Выбираем
            number = floor(number1 * 101) !случайный вихрь
            !Меняем координаты вихря, используя полярные координаты
            call random_number(randd)
            poschange = randd * lambda / 10.0 *10.0**(4)
            call random_number(randd)
            angle = 2 * pi * randd
            xchange = poschange * cos(angle)
            ychange = poschange * sin(angle)
             !записали в массив старые значения положения вихрей
            !Меняем координату х выбранного случайного вихря
            if(r2(number, 1) + xchange < 0) then
                temp = r2(number, 1) + xchange + 5
                r2(number, 1) = temp
            else if (r2(number, 1) + xchange > 5) then
                temp = r2(number, 1) + xchange - 5
                r2(number, 1) = temp
            else 
                temp = r2(number, 1) + xchange
                r2(number, 1) = temp
            end if
            !Меняем координату y выбранного случайного вихря
            if(r2(number, 2) + ychange < 0) then
                temp = r2(number, 2) + ychange + 5
                r2(number, 2) = temp
            else if (r2(number, 2) + ychange > 5) then
                temp = r2(number, 2) + ychange - 5
                r2(number, 2) = temp
            else 
                temp = r2(number, 2) + ychange
                r2(number, 2) = temp
            end if
            Unewone = 0.0
            !Считаем энергию переставленного вихря(используем массив координат r2, где мы только что поменли его координату)
            do j = 1, 100         !
                if( j.NE. number ) then
                    x1 = r2(j, 1)
                    x2 = r2(number, 1)
                    y1 = r2(j, 2)
                    y2 = r2(number, 2)
                    if(y1>y2) then
                        if(5-y1+y2> y1-y2) then
                            ydiff = y1-y2
                        else
                            ydiff = 5-y1+y2
                        end if
                    else
                        if(5-y2+y1> y2-y1) then
                            ydiff = y2-y1
                        else
                            ydiff = 5-y2+y1
                        end if
                    end if
                    if(x1>x2) then
                        if(5-x1+x2> x1-x2) then
                            xdiff = x1-x2
                        else
                            xdiff = 5-x1+x2
                        end if
                    else
                        if(5-x2+x1> x2-x1) then
                            xdiff = x2-x1
                        else
                            xdiff = 5-x2+x1
                        end if
                    end if   
                    rij = sqrt(ydiff**2 + xdiff**2)
                    !Приводим к нужной размерности
                    x = rij*10.0**(-4) 
                    !Считаем энергию с помощью табулированного массива
                    if(x < 7.0 * lambda) then 
                        i4 = ceiling(x*10000.0/(7.0*lambda)) 
                        Unewone = Unewone + Utab(i4)
                    end if
                end if
            end do
            Uoldone = 0.0
            !Считаем энергию вихря до перестановки, используем старый массив координат r1
            !(Не важно, что сначала считали после перестановки, потом до, тк ассивы разные. 
            !Проверял, меняя местами эти циклы, результат один и тот же)
            do j1 = 1, 100         
                if( j1.NE. number ) then
                    x1 = r1(number, 1)
                    x2 = r1(j1, 1)
                    y1 = r1(number, 2)
                    y2 = r1(j1, 2)
                    if(y1>y2) then
                        if(5-y1+y2> y1-y2) then
                            ydiff = y1-y2
                        else
                            ydiff = 5-y1+y2
                        end if
                    else
                        if(5-y2+y1> y2-y1) then
                            ydiff = y2-y1
                        else
                            ydiff = 5-y2+y1
                        end if
                    end if
                    if(x1>x2) then
                        if(5-x1+x2> x1-x2) then
                            xdiff = x1-x2
                        else
                            xdiff = 5-x1+x2
                        end if
                    else
                        if(5-x2+x1> x2-x1) then
                            xdiff = x2-x1
                        else
                            xdiff = 5-x2+x1
                        end if
                    end if   
                    rij = sqrt(ydiff**2 + xdiff**2)
                    !приводим к нужной размерности
                    x = rij*10.0**(-4) 
                    !вычисляем энергию при помощи табулированного массива
                    if(x < 7.0 * lambda) then 
                        i4 = ceiling(x*10000.0/(7.0*lambda)) 
                        Uoldone = Uoldone + Utab(i4)
                    end if
                end if
            end do
            !Главный шаг - вычисляем энергию после перестановки
            Uj = Uj -Uoldone + Unewone
                !Если энергия меньше чем была до перестановки, коныигурация принимается
            if (Uj < Ui) then
                r1 = r2
                U0 = U0+Uj - Ui
                Ui = Uj
            else
                !В противном случае считаем вероятность(Ошибка не в ней, тк пробовал просто убрать ее и считать, что если энергия больше
                !то просто отвергаем конфигурацию, но ошибка никуда не ушла)
            possib1 = 1.0 / e ** ( (Uj-Ui) * 6.2 * 10.0**(11) * 11606.0 / T)
            call random_number(randd)
                if(possib1>randd) then
                    r1 = r2
                    U0 = U0+Uj - Ui
                    Ui = Uj
        
                else
                    !Конфигурация не принята - возврат к исходной энергии
                    Uj = Ui
                    !Пробовал и такое выражение Uj = Uj + Uoldone - Unewone
                    end if
                end if
            end if
        end do

        !вывод жнергии каждые 10 итераций
        if (mod(i1,10)==0) then
            write(*,*) i, Ui , U0 , Uj
        end if
        !запись энергии в файл
        write(4,*) i, Ui * 6.2 * 10.0**(11)
        !наметки тепловой карты, на алгоритм не влияют, пробовал удалять их
        if(i>179 .and. mod(i, 5)==0) then
            do sec1 = 1, 100
                write(2,*) r1(i2,:)
            end do
            write(2, *) NEW_LINE('a'), NEW_LINE('a')    
        end if
        if(i>179) then
            do sec = 1, 100
                sec10 = floor(r1(sec, 1) * 8.0)
                sec20 = floor(r1(sec, 2)*8.0)
                sect(sec10, sec20) = sect(sec10, sec20) + 1
            end do
        end if
    end do
    




    







    !write(*,*) sect 
    do i2 = 1, 100
        write(3,*) r1(i2,:)
    end do
    !100 stop "Ошибка при открытии файла"
    !101 stop "Ошибка при чтении данных"
    !102 stop "Ошибка при создании файла"
    !103 stop "Ошибка при записи данных"
    close(2)
    close(3)
    close(4)
    sect = sect / 21.0
    open(5, file="heatmap_task1_data.txt")
    write(5, *) 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
    do i3 = 1, 20
        write(5, *) i3, sect(i3, :)
    end do
    close(5)

    Ue = 0
    do j = 1, 100         !Считаем
        do k = 1, 100      !начальную энергию этой итерации
            if( j .NE. k) then
                x1 = r1(j, 1)
                x2 = r1(k, 1)
                y1 = r1(j, 2)
                y2 = r1(k, 2)
                if(y1>y2) then
                    if(5-y1+y2> y1-y2) then
                        ydiff = y1-y2
                    else
                        ydiff = 5-y1+y2
                    end if
                else
                    if(5-y2+y1> y2-y1) then
                        ydiff = y2-y1
                    else
                        ydiff = 5-y2+y1
                    end if
                end if
                if(x1>x2) then
                    if(5-x1+x2> x1-x2) then
                        xdiff = x1-x2
                    else
                        xdiff = 5-x1+x2
                    end if
                else
                    if(5-x2+x1> x2-x1) then
                        xdiff = x2-x1
                    else
                        xdiff = 5-x2+x1
                    end if
                end if   
                rij = sqrt(ydiff**2 + xdiff**2)
                x = rij*10.0**(-4)
                if(x < 7 * lambda) then 
                    i4 = ceiling(x*10000.0/(7.0*lambda)) 
            
                    temp = Utab(i4)
                    Ue(j) = Ue(j) + temp
                end if
            end if
        end do
    end do
    Uj = sum(Ue) / 2.0
    write(*,*) Uj
    stop
    
    !plot 'C:/Users/Valera P. Lenkov/Desktop/Fortran/task1/res2.txt'
end program task1