!!!ГЕНЕРАЦИЯ ДЕФЕКТОВ 
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
    real(kind = 8), parameter :: e =  2.71828182845904, c = 3e8
    real(kind = 8), parameter :: phi0 =  2.06783460999999E-07            !Квант магнитного потока
    real(kind = 8), parameter :: delta = 0.27       !Толщина слоя
    real(kind = 8) :: xi0, xi          !Длина когерентности сверхпроводника
    real(kind = 8), parameter :: lambda0 = 180.0         !Лондоновская глубина проникновения данного материала при Т=0 
    real(kind = 8) :: lambda             !Лондоновская глубина проникновения данного материала при данной температуре
    real(kind = 8) :: T, Tc
    real(kind = 8) :: length, width, tabul, defsnet(4000, 4000),stepw,stepl,stepdef
    real(kind=8):: defs(1000,2), Uptab(10000), pi, tauNa, coordx, coordy
    real(kind=8):: x1,x2,y1,y2,xdiff,ydiff,rad, Up, radmin
    integer:: M, idefs, i, j
    
    tauNa = 100.0
    T = 1.0
    xi0 = 2.0


    tabul = 10000.0
    length = 5.0*10.**3 ! по х
    width = 5.0*10.**3   !по у
    Tc = 84.0
    xi = xi0 * (1.0 - (T / Tc) ** 3.3) ** (-1.0 / 2.0)
    lambda = lambda0 * (1.0 - (T / Tc) ** (3.3)) ** (-1.0 / 2.0)
    pi = acos(-1.0)
    T= T / 11606.0 ! T в эВ
    stepw = width / 2000
    stepl = length / 2000
    open(77, file='defs_squared_1.txt')
    read(77, *) M
    do idefs = 1, M
        read(77,*) defs(idefs,:)
    end do
    close(77)

    open(79, file='Uptabpy_10K.txt')
    read(79, *) Uptab
    close(79)
    stepdef = 5.*xi/tabul

    do i = 1, 4000
        do j =1, 4000
            Up = 1.0
            coordx = length / 4000 * i - length / 2.
            coordy = width / 4000 * j - width / 2.
            do idefs = 1, M
                y1 = defs(idefs, 2)
                y2 = coordy
                x1 = defs(idefs, 1)
                x2 = coordx
                ydiff = min(abs(y1-y2), width - abs(y1 - y2))
                xdiff = abs(x1-x2)  
                rad = sqrt(ydiff**2 + xdiff**2)
                if(rad<5*xi) then
                    if(Up>0.) then
                        Up = Uptab(ceiling(rad/stepdef))
                        radmin=rad
                    else if (rad<radmin) then
                    radmin = rad
                    Up = Uptab(ceiling(rad/stepdef))
                    end if
                end if
            end do
            if(Up>0.) then
                defsnet(i, j) = 0.
            else
                defsnet(i,j) = Up
                if(Up<0.) then
                write(*,*) Up
                end if
            end if
        end do
    end do
    open(2, file='defsnet.txt')
    do i = 1, 4000
        write(2, *) defsnet(:, i)
    end do
    close(2)
end program