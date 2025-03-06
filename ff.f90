
program task1
    use bessel
    implicit none
    
    !real, parameter :: phi0 =
    !real, parameter :: delta = 
    !real, parameter :: xi0 = 
    !real, parameter :: lambda0 = 
    real ::b(15), c(100000000, 15), v(100000000)
    integer ::  j , k, co
    integer :: i, m
    real :: d
    
    do i = 1, 100000000
        k = 0
        do while(j<16)
            call random_number(d)
            if(d*10>1) then
                k = k + 1
            else if(k > 0) then
                b(j) = k
                k = 0
                j = j + 1
            end if
        end do
        c(i, :) = b
        b = 0
        j = 1
    end do
    co = 0
    v = maxval(c, 2, c>25)
    do i = 1, 100000000
        
        if(v(i) > 0.1) then
            co = co + 1
        end if
    end do
    
    i = bk0(5)
    write(*,*) co / 100000000.0 
end program task1