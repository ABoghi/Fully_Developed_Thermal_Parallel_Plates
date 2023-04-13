!!!************************************************************
!!!*						                	                *
!!!*                     K Epsilon 1D      	         	    *
!!!*							        	                    *
!!!*                Author: Dr. Andrea Boghi  		        *
!!!*							        	                    *
!!!*								                            *
!!!************************************************************
    
Program main_K_epsilon
    implicit none
    real*8,allocatable :: y(:),U(:),kt(:),eps(:),detady(:),d2etady2(:)
    real*8,allocatable :: U0(:),kt0(:),eps0(:)
    real*8,allocatable :: nut(:),f2(:),dnutdy(:),dUdy(:),dnutdeta(:),dUdeta(:)
    real*8,allocatable :: tau_mu(:),tau_R(:),Pk(:),Tk(:),Dk(:)
    real*8,allocatable :: Peps(:),Teps(:),Deps(:),epseps(:)
    real*8,allocatable :: T(:),Th2(:),T0(:),Th20(:),lambda(:),dlambdadT(:),d2lambdadT2(:)
    real*8,allocatable :: dTdeta(:),dTdy(:),d2Tdeta2(:),d2Tdy2(:),dTh2deta(:),dTh2dy(:),d2Th2deta2(:),d2Th2dy2(:)
    real*8,allocatable :: q_lam(:),q_R(:),q_new(:),P_Th2(:),eps_Th2(:),T_Th2(:),D_Th2(:),H_Th2(:)
    integer j,ny,nhy,iter,niter
    real*8 Re_tau,Pr,Bp,Cp,dy_min,sigmak,sigmae,Cmu,Ce1,Ce2,f1,alphaU,alphaKt,alphaeps
    real*8 resU,resK,resE,resT,resTh2,deta,aU_w,aU_e,sU,aK_w,aK_e,sK,aE_w,aE_e,sE, conv_fac
    real*8 resU_old,resK_old,resE_old,resT_old,resTh2_old
    real*8 aT_e,aT_w,sT,aTh2_e,aTh2_w,sTh2
    real*8 sigmaT,sigmaTh2,alphaT,alphaTh2,U_max
    CHARACTER(len=80)::fname_ke
    CHARACTER(len=80)::fname_th
    CHARACTER(len=80)::fname_res
    logical flag

    open(1,file='imp_ke_th_var.dat')
    read(1,*) flag
    read(1,*) nhy
    read(1,*) niter
    read(1,*) Re_tau
    read(1,*) dy_min
    read(1,*) Pr
    read(1,*) Bp
    read(1,*) Cp
    read(1,*) alphaU
    read(1,*) alphaKt
    read(1,*) alphaeps 
    read(1,*) alphaT
    read(1,*) alphaTh2
    close(1)

    ny = nhy*2

    print*, ' niter =', niter, ' ny =', ny 
    print*, ' Re_tau =', Re_tau, ' flag =', flag 

    allocate(y(1:ny),U(1:ny),kt(1:ny),eps(1:ny),nut(1:ny),dnutdy(1:ny),dUdy(1:ny),f2(1:ny))
    allocate(U0(1:ny),kt0(1:ny),eps0(1:ny),detady(1:ny),d2etady2(1:ny),dnutdeta(1:ny),dUdeta(1:ny))
    allocate(T(1:ny),Th2(1:ny),T0(1:ny),Th20(1:ny),lambda(1:ny),dlambdadT(1:ny),d2lambdadT2(1:ny))
    allocate(dTdeta(1:ny),dTdy(1:ny),d2Tdeta2(1:ny),d2Tdy2(1:ny),dTh2deta(1:ny),dTh2dy(1:ny),d2Th2deta2(1:ny),d2Th2dy2(1:ny))
    allocate(q_lam(1:ny),q_R(1:ny),q_new(1:ny),P_Th2(1:ny),eps_Th2(1:ny),T_Th2(1:ny),D_Th2(1:ny),H_Th2(1:ny))

    call initialization(flag,ny,Re_tau,Pr,Bp,Cp,dy_min,y,detady,d2etady2,U,Kt,eps,T,Th2,deta)

    call nagano_takawa_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu,f1)

    sigmaT = 1.d0
    sigmaTh2 = 1.d0

    conv_fac = 1.d0

    write(fname_res,110)'residualsNT_var_Re_tau=',int(Re_tau),'_log10(Pr)=',int(log10(Pr)),'.csv'
    110 format(a23,i3,a11,i2,a4)

    open(11,file=fname_res)
    write(11,*) '"iter","resU","resK","resE","resT","resTh2"'

    do iter=1,niter

        resU_old = resU
        resK_old = resK
        resE_old = resE
        resT_old = resT
        resTh2_old = resTh2

        U0 = U
        Kt0 = kt
        eps0 = eps
        T0 = T
        Th20 = Th2

        do j=1,ny
            call thernal_diffusivity(lambda(j),dlambdadT(j),d2lambdadT2(j),T(j),Bp,Cp)
        enddo

        call nagano_takawa_k_epsilon_functions(nut,f2,ny,y,kt,eps,Cmu)
        call ddeta(ny,nut,dnutdeta,deta)
        dnutdy = dnutdeta*detady
        call ddeta(ny,U,dUdeta,deta)
        dUdy = dUdeta*detady
        call ddeta(ny,T,dTdeta,deta)
        dTdy = dTdeta*detady
        call ddeta(ny,Th2,dTh2deta,deta)
        dTh2dy = dTh2deta*detady
        call d2deta2(ny,T,d2Tdeta2,deta)
        d2Tdy2 = d2Tdeta2*detady**2.d0 + dTdeta*d2etady2
        call d2deta2(ny,Th2,d2Th2deta2,deta)
        d2Th2dy2 = d2Th2deta2*detady**2.d0 + dTh2deta*d2etady2

        U_max = maxval(U, dim=1, mask=(U>0))

        U(1) = 0.d0
        Kt(1) = 0.d0
        eps(1) = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(1)))+4.d0*dsqrt(dabs(Kt(2)))-dsqrt(dabs(Kt(3))))/(2.d0*deta) )*detady(1) )**2.d0
        call T_coefficients(aT_w,aT_e,sT,nut(1),dnutdy(1),lambda(1),dlambdadT(1),d2lambdadT2(1),dTh2dy(1),d2Th2dy2(1),dTdy(1), &
            Pr,sigmaT,deta,d2etady2(1),detady(1))
        T(1) = sT + aT_e*T(2) + aT_w*( T(2) - 2.d0 * deta * Pr / ( lambda(1) * detady(1) ) )
        !!!T(1) = Pr * Re_tau ! U_max
        Th2(1) = 0.d0

        do j =2,ny-1
            call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdy(j),deta,Re_tau,d2etady2(j),detady(j))
            U(j) =  sU + aU_e*U(j+1) + aU_w*U(j-1)
            call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdy(j),dUdy(j),deta,sigmak,d2etady2(j),detady(j))
            Kt(j) = sK + aK_e*Kt(j+1) + aK_w*Kt(j-1)
            call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdy(j),dUdy(j),deta,sigmae,Ce1,f1,Ce2,f2(j),d2etady2(j), &
            detady(j))
            eps(j) = sE + aE_e*eps(j+1) + aE_w*eps(j-1)
            call T_coefficients(aT_w,aT_e,sT,nut(j),dnutdy(j),lambda(j),dlambdadT(j),d2lambdadT2(j),dTh2dy(j),d2Th2dy2(j),dTdy(j), &
            Pr,sigmaT,deta,d2etady2(j),detady(j))
            T(j) = sT + aT_e*T(j+1) + aT_w*T(j-1)
            call Th2_coefficients(aTh2_w,aTh2_e,sTh2,nut(j),dnutdy(j),lambda(j),dlambdadT(j),d2lambdadT2(j),dTdy(j),d2Tdy2(j), &
            Pr,sigmaTh2,deta,d2etady2(j),detady(j),eps(j),Kt(j))
            Th2(j) = sTh2 + aTh2_e*Th2(j+1) + aTh2_w*Th2(j-1)
        enddo
        U(ny) = 0.d0
        Kt(ny) = 0.d0
        eps(ny) = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(ny)))+4.d0*dsqrt(dabs(Kt(ny-1)))-dsqrt(dabs(Kt(ny-1))))/(2.d0*deta) & 
        )*detady(ny) )**2.d0
        call T_coefficients(aT_w,aT_e,sT,nut(ny),dnutdy(ny),lambda(ny),dlambdadT(ny),d2lambdadT2(ny),dTh2dy(ny),d2Th2dy2(ny), & 
        dTdy(ny), Pr,sigmaT,deta,d2etady2(ny),detady(ny))
        T(ny) = sT + aT_w*T(ny-1) + aT_e*( T(ny-1) + 2.d0 * deta * Pr / ( lambda(ny) * detady(ny) ) ) 
        !!!T(ny) = - T(1)
        Th2(ny) = 0.d0

        call residuals(ny,U,U0,resU)
        call residuals(ny,Kt,Kt0,resK)
        call residuals(ny,eps,eps0,resE)
        call residuals(ny,T,T0,resT)
        call residuals(ny,Th2,Th20,resTh2)
        write(11,102) conv_fac*iter,',',resU,',',resK,',',resE,',',resT,',',resTh2

        U = dabs(alphaU*U +(1.d0-alphaU)*U0)
        Kt = dabs(alphaKt*Kt +(1.d0-alphaKt)*Kt0)
        eps = dabs(alphaeps*eps +(1.d0-alphaeps)*eps0)
        T = dabs(alphaT*T +(1.d0-alphaT)*T0)
        Th2 = dabs(alphaTh2*Th2 +(1.d0-alphaTh2)*Th20)
        print*, ' completed =', 100*real(iter)/real(niter), ' resU = ', resU, ' resK = ', resK, ' resE = ', resE, &
        ' resT = ', resT, ' resTh2 = ', resTh2

        !if(iter > 1) then
        !    call correct_residuals(alphaU,resU,resU_old)
        !    call correct_residuals(alphaKt,resK,resK_old)
        !    call correct_residuals(alphaeps,resE,resE_old)
        !    call correct_residuals(alphaT,resT,resT_old)
        !    call correct_residuals(alphaTh2,resTh2,resTh2_old)
        !endif
        
    enddo
    close(11)

    open(1,file='point_ke_var.dat',form='unformatted')
    write(1) y,detady,d2etady2,U,Kt,eps,nut,f2
    close(1)

    open(1,file='point_T_var.dat',form='unformatted')
    write(1) y,T,Th2
    close(1)

    allocate(tau_mu(1:ny),tau_R(1:ny),Pk(1:ny),Tk(1:ny),Dk(1:ny))
    allocate(Peps(1:ny),Teps(1:ny),Deps(1:ny),epseps(1:ny))

    call output_fields_k_eps(ny,U,Kt,eps,nut,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps,epseps,detady, &
    d2etady2)

    write(fname_ke,111)'momentumNT_var_Re_tau=',int(Re_tau),'_log10(Pr)=',int(log10(Pr)),'.csv'
    111 format(a22,i3,a11,i2,a4)

    open(14,file=fname_ke,form='formatted')
    write(14,*) '"y","U","k","epsilon","nut","tau_mu","tau_R","Pk","Tk","Dk","Peps","Teps","Deps","epsEps"'
    do j=1,ny
       write(14,101) y(j),',',U(j),',',Kt(j),',',eps(j),',',nut(j),',',tau_mu(j),',',tau_R(j),',',Pk(j),',',Tk(j),',',Dk(j),',', &
       Peps(j),',',Teps(j),',',Deps(j),',',epsEps(j)
    enddo
    close(14)

    call output_fields_thermal(ny,T,Th2,lambda,dlambdadT,d2lambdadT2,Kt,eps,nut,detady,d2etady2,deta,sigmaT,sigmaTh2,Pr, &
    q_lam,q_R,q_new,P_Th2,eps_Th2,T_Th2,D_Th2,H_Th2)

    write(fname_th,112)'thermalNT_var_Re_tau=',int(Re_tau),'_log10(Pr)=',int(log10(Pr)),'.csv'
    112 format(a21,i3,a11,i2,a4)

    open(15,file=fname_th,form='formatted')
    write(15,*) '"y","T","Th2","q_lam","q_R","q_new","P_Th2","eps_Th2","T_Th2","D_Th2","H_Th2"'
    do j=1,ny
       write(15,103) y(j),',',T(j),',',Th2(j),',',q_lam(j),',',q_R(j),',',q_new(j),',',P_Th2(j),',',eps_Th2(j),',',T_Th2(j), & 
       ',',D_Th2(j),',',H_Th2(j)
    enddo
    close(15)

    101 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10, &
    A,e18.10)

    102 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    103 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    end

!!!*************************************************
!!!*						         	             *
!!!*            initialization                        *
!!!*								                 *
!!!*************************************************

subroutine initialization(flag,ny,Re_tau,Pr,Bp,Cp,dy_min,y,detady,d2etady2,U,Kt,eps,T,Th2,deta)
    implicit none
    logical, intent(in) :: flag
    integer, intent(in) :: ny
    real*8, intent(in) :: Re_tau,dy_min,Pr,Bp,Cp
    real*8, intent(out) :: y(1:ny),detady(1:ny),d2etady2(1:ny),U(1:ny),kt(1:ny),T(1:ny),Th2(1:ny),eps(1:ny),deta
    integer j
    real*8 Kappa, Cmu,nut(1:ny),f2(1:ny),y_mid

    Kappa = 4.d-1
    Cmu = 9.d-2
    y_mid = 11.635

    !initial conditions
    if (flag) then
        
        print *,'continuation'
  
        open(1,file='point_ke_var.dat',form='unformatted')
        read(1) y(1:ny),detady(1:ny),d2etady2(1:ny),U(1:ny),Kt(1:ny),eps(1:ny),nut(1:ny),f2(1:ny)
        close(1)        
        
        open(1,file='point_T_var.dat',form='unformatted')
        read(1) y(1:ny),T(1:ny),Th2(1:ny)
        close(1)
  
    else
  
        print *,'new run'

        call grid(ny,dy_min,Re_tau,y,detady,d2etady2,deta)

        do j=1,ny/2
            if(y(j)<=y_mid) then
                U(j) = y(j)
                eps(j) = (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)/(Kappa*y_mid)**2.d0
                Kt(j) = dsqrt( (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)*eps(j)/Cmu )
            else
                U(j) = (1.d0/Kappa)*dlog(y(j)) +5.5d0 
                eps(j) = (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)/(Kappa*y(j))**2.d0
                Kt(j) = dsqrt( (y(j)*Kappa -y(j)*y(j)*Kappa/Re_tau)*eps(j)/Cmu )
            endif
        enddo

        do j=ny/2+1,ny
            U(j) = U(ny+1-j)
            eps(j) = eps(ny+1-j)
            Kt(j) = Kt(ny+1-j)
        enddo

        do j=1,ny
            Th2(j) = Kt(j)
        enddo

        !!!********************************************************
        !!!
        !!! dy^+ = -K^+(T)dT^+
        !!!
        !!!********************************************************
  
        do j=1,ny/2
            T(j) = - Pr*( U(j) - U(ny/2) )
        enddo

        do j=ny/2+1,ny
            T(j) = -T(ny+1-j)
        enddo

    endif

    end

!!!*************************************************
!!!*						         	             *
!!!*                    ddeta                        *
!!!*								                 *
!!!*************************************************
    
subroutine  ddeta(ny,A,DA,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),deta
    real*8, intent(out) :: DA(1:ny)
    integer j

    DA(1) =  -(3.d0*A(1) -4.d0*A(2) +A(3))/(2.d0*deta)

    do j=2,ny-1
        DA(j) = (A(j+1) -A(j-1))/(2.d0*deta)
    enddo

    DA(ny) = (3.d0*A(ny) -4.d0*A(ny-1) +A(ny-2))/(2.d0*deta)

    end

!!!*************************************************
!!!*						         	             *
!!!*                    d2deta2                        *
!!!*								                 *
!!!*************************************************

subroutine  d2deta2(ny,A,D2A,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),deta
    real*8, intent(out) :: D2A(1:ny)
    real*8 deta2
    integer j

    deta2 = deta*deta
    
    D2A(1) =  (12.d0*a(1) -30.d0*a(2) +24.d0*a(3) -6.d0*a(4))/(6.d0*deta2)
    
    do j=2,ny-1
        D2A(j) = (a(j+1) - 2.d0*a(j) + a(j-1))/deta2
    enddo
    
    D2A(ny) = (12.d0*a(ny) -30.d0*a(ny-1) +24.d0*a(ny-2) -6.d0*a(ny-3))/(6.d0*deta2)
    
    end

!!!***************************************************
!!!*						         	               *
!!!*       Nagano Takawa K - Epsilon Constants 	       	   *
!!!*								                   *
!!!***************************************************
    
subroutine  nagano_takawa_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu,f1)
    implicit none
    real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu,f1

    sigmaK= 1.0d0
    sigmae= 1.3d0
    Ce1=1.45d0
    Ce2=1.9d0
    Cmu=0.09d0
    f1=1.d0
    
    end

!!!***************************************************
!!!*						         	               *
!!!*       Nagano Takawa K - Epsilon Functions 	       	   *
!!!*								                   *
!!!***************************************************
    
subroutine  nagano_takawa_k_epsilon_functions(nut,f2,ny,y,kt,eps,Cmu)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: y(1:ny),Kt(1:ny),eps(1:ny),Cmu
    real*8, intent(out) :: nut(1:ny),f2(1:ny)
    real*8 Ret(1:ny), fmu(1:ny), Ret_min, eps_min
    integer j

    eps_min = 1.d-12

    do j=1,ny
        if (eps(j) <= eps_min) then
            Ret(j)= dabs(Kt(j)*Kt(j)/eps_min)
        else
            Ret(j)= dabs(Kt(j)*Kt(j)/eps(j))
        endif
    enddo

    Ret_min = 1.d-12
    do j=1,ny/2
        if (Ret(j) <= Ret_min) then
            fmu(j)= (1.d0 +4.1d0/Ret_min**0.75d0)*(1.d0 -dexp(-y(j)/26.d0))**2.d0
        else
            fmu(j)= (1.d0 +4.1d0/Ret(j)**0.75d0)*(1.d0 -dexp(-y(j)/26.d0))**2.d0
        endif
    enddo

    do j=ny/2+1,ny
        if (Ret(j) <= Ret_min) then
            fmu(j)= (1.d0 +4.1d0/Ret_min**0.75d0)*(1.d0 -dexp(-y(ny+1-j)/26.d0))**2.d0
        else
            fmu(j)= (1.d0 +4.1d0/Ret(j)**0.75d0)*(1.d0 -dexp(-y(ny+1-j)/26.d0))**2.d0
        endif
    enddo

    do j=1,ny
        nuT(j)= Cmu*fmu(j)*Ret(j)
    enddo

    do j=1,ny/2
        f2(j)= (1.d0 -0.3d0*dexp(-(Ret(j)/6.5d0)**2.d0))*(1.d0 -dexp(-y(j)/6.d0))**2.d0
    enddo

    do j=ny/2+1,ny
        f2(j)= (1.d0 -0.3d0*dexp(-(Ret(j)/6.5d0)**2.d0))*(1.d0 -dexp(-y(ny+1-j)/6.d0))**2.d0
    enddo
    
    end

!!!***************************************************
!!!*						         	               *
!!!*                U coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  u_coefficients(aU_w,aU_e,sU,nut,dnutdy,deta,Re_tau,d2etady2,detady)
    implicit none
    real*8, intent(in) :: nut,dnutdy,deta,Re_tau,d2etady2,detady
    real*8, intent(out) :: aU_w,aU_e,sU
    real*8 dev

    dev = deta*( (1.d0+nut)*d2etady2 + dnutdy*detady )/(4.d0*(1.d0+nut)*(detady)**2.d0)

    aU_w = 5.d-1 - dev
    aU_e = 5.d-1 + dev
    sU = (deta*deta)/(2.d0*Re_tau*(1.d0+nut)*(detady)**2.d0)

    end

!!!***************************************************
!!!*						         	               *
!!!*                K coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  K_coefficients(aK_w,aK_e,sK,eps,nut,dnutdy,dUdy,deta,sigmak,d2etady2,detady)
    implicit none
    real*8, intent(in) :: eps,nut,dnutdy,dUdy,deta,sigmak,d2etady2,detady
    real*8, intent(out) :: aK_w,aK_e,sK
    real*8 dev

    dev = deta*( (sigmak+nut)*d2etady2 + dnutdy*detady )/(4.d0*(sigmak+nut)*detady**2.d0)

    aK_w = 5.d-1 - dev
    aK_e = 5.d-1 + dev
    sK = (nut*dUdy*dUdy - eps)*(deta*deta)/(2.d0*(1.d0+nut/sigmak)*detady**2.d0)

    end

!!!***************************************************
!!!*						         	               *
!!!*                E coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  E_coefficients(aE_w,aE_e,sE,eps,Kt,nut,dnutdy,dUdy,deta,sigmae,Ce1,f1,Ce2,f2,d2etady2,detady)
    implicit none
    real*8, intent(in) :: eps,Kt,nut,dnutdy,dUdy,deta,sigmae,Ce1,f1,Ce2,f2,d2etady2,detady
    real*8, intent(out) :: aE_w,aE_e,sE
    real*8 K_min, Kb, dev
    logical method1

    dev = deta*( (sigmae+nut)*d2etady2 + dnutdy*detady )/(4.d0*(sigmae+nut)*detady**2.d0)

    K_min = 1.d-12

    Kb = (Ce1*f1*nut*dUdy*dUdy -Ce2*f2*eps)*(deta*deta/(2.d0*(1.d0+nut/sigmae)*detady**2.d0))

    method1 = .true.
    if (method1) then
        aE_w = 5.d-1 - dev
        aE_e = 5.d-1 + dev
        if (Kt<=K_min) then
            sE = 0.d0*Kb*eps/K_min
        else
            sE = Kb*eps/Kt
        endif
    else
        aE_w = (5.d-1 - dev)/(1.d0 - Kb/Kt)
        aE_e = (5.d-1 + dev)/(1.d0 - Kb/Kt)
        sE = 0.d0 
    endif
    
    end

!!!***************************************************
!!!*						         	             *
!!!*                T coefficients	       	         *
!!!*								                 *
!!!***************************************************

subroutine  T_coefficients(aT_w,aT_e,sT,nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTh2dy,d2Th2dy2,dTdy, &
    Pr,sigmaT,deta,d2etady2,detady)
    implicit none
    real*8, intent(in) :: nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTh2dy,d2Th2dy2,dTdy,Pr,sigmaT,deta,d2etady2,detady
    real*8, intent(out) :: aT_w,aT_e,sT
    real*8 A1,A2,A3, dev

    A1 = lambda / Pr + nut / sigmaT
    A2 = ( dlambdadT * dTdy + d2lambdadT2 * dTh2dy ) / Pr + dnutdy / sigmaT
    A3 =  dlambdadT * d2Th2dy2 / Pr

    !print*, ' A1 = ',A1,' A2 = ',A2,' A3 = ',A3

    dev = deta * ( d2etady2 / detady + A2 / A1 ) / ( 4.d0 * detady )

    aT_w = ( 5.d-1 - dev )
    aT_e = ( 5.d-1 + dev )
    sT = ( ( deta * deta ) / ( 2.d0 * detady**2.d0 ) ) * A3 / A1

    end

!!!***************************************************
!!!*						         	             *
!!!*                Th2 coefficients	       	     *
!!!*								                 *
!!!***************************************************

subroutine  Th2_coefficients(aTh2_w,aTh2_e,sTh2,nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTdy,d2Tdy2, &
    Pr,sigmaTh2,deta,d2etady2,detady,eps,Kt)
    implicit none
    real*8, intent(in) :: nut,dnutdy,lambda,dlambdadT,d2lambdadT2,dTdy,d2Tdy2,Pr,sigmaTh2,deta,d2etady2,detady
    real*8, intent(in) :: eps,Kt
    real*8, intent(out) :: aTh2_w,aTh2_e,sTh2
    real*8 A1,A2,A3,A4,dev, den

    A1 = lambda / Pr + nut / sigmaTh2
    A2 = 2.d0 * dlambdadT * dTdy / Pr + dnutdy / sigmaTh2
    A3 =  ( nut / sigmaTh2 ) * dTdy**2.d0
    A4 = ( ( eps / Kt ) * lambda - 2.d0 * dlambdadT * d2Tdy2 - 2.d0 *d2lambdadT2 * dTdy**2.d0 ) / Pr 

    dev = deta * ( d2etady2 / detady + A2 / A1 ) / ( 4.d0 * detady )

    den = ( 2.d0 * detady**2.d0 + deta * deta * ( A4 / A1 ) )

    aTh2_w = ( 5.d-1  - dev ) * ( 2.d0 * detady**2.d0 ) / den
    aTh2_e = ( 5.d-1  + dev ) * ( 2.d0 * detady**2.d0 ) / den
    sTh2 = deta * deta * ( A3 / A1 ) / den

    end

!!!*************************************************
!!!*						         	           *
!!!*          output fields k-eps                  *
!!!*								               *
!!!*************************************************

subroutine  output_fields_k_eps(ny,U,Kt,eps,nut,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps,epseps, &
    detady, d2etady2)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: deta,sigmaK,sigmaE,Ce1,Ce2,f1,f2(1:ny),Kt(1:ny),eps(1:ny),nut(1:ny),U(1:ny),detady(1:ny),d2etady2(1:ny)
    real*8, INTENT(OUT) :: tau_mu(1:ny),tau_R(1:ny),Pk(1:ny),Tk(1:ny),Dk(1:ny)
    real*8, INTENT(OUT) :: Peps(1:ny),Teps(1:ny),Deps(1:ny),epseps(1:ny)
    real*8 dUdy(1:ny),d2Ktdeta2(1:ny),d2epsdeta2(1:ny),dKtdeta(1:ny),depsdeta(1:ny),dnutdy(1:ny),dUdeta(1:ny),dnutdeta(1:ny)

    call ddeta(ny,U,dUdeta,deta)
    dUdy = dUdeta*detady
    call ddeta(ny,nut,dnutdeta,deta)
    dnutdy = dnutdeta*detady

    tau_mu = dUdy
    tau_R = nut*dUdy
    Pk = tau_R*dUdy

    call d2deta2(ny,Kt,d2Ktdeta2,deta)
    call ddeta(ny,Kt,dKtdeta,deta)

    Dk = d2Ktdeta2*detady**2.d0 + dKtdeta*d2etady2
    Tk = (nut/sigmaK)*Dk + (dKtdeta*detady/sigmaK)*dnutdy

    Peps(2:ny) = f1*Ce1*(eps(2:ny)/Kt(2:ny))*Pk(2:ny)
    Peps(1) = Peps(2)

    call d2deta2(ny,eps,D2epsdeta2,deta)
    call ddeta(ny,eps,depsdeta,deta)

    Deps = d2epsdeta2*detady**2.d0 + depsdeta*d2etady2
    Teps = (nut/sigmaE)*Deps + (depsdeta*detady/sigmaE)*dnutdy
    epsEps(2:ny-1) = -f2(2:ny-1)*Ce2*(eps(2:ny-1)/Kt(2:ny-1))*eps(2:ny-1)
    epsEps(1) = epsEps(2)
    epsEps(ny) = epsEps(ny-1)
    
    end

!!!*************************************************
!!!*						         	           *
!!!*          output fields thermal                  *
!!!*								               *
!!!*************************************************

subroutine  output_fields_thermal(ny,T,Th2,lambda,dlambdadT,d2lambdadT2,Kt,eps,nut,detady,d2etady2,deta,sigmaT,sigmaTh2,Pr, &
    q_lam,q_R,q_new,P_Th2,eps_Th2,T_Th2,D_Th2,H_Th2)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: deta,Kt(1:ny),eps(1:ny),nut(1:ny),detady(1:ny),d2etady2(1:ny)
    real*8, intent(in) :: T(1:ny),Th2(1:ny),lambda(1:ny),dlambdadT(1:ny),d2lambdadT2(1:ny),sigmaT,sigmaTh2,Pr
    real*8, INTENT(OUT) :: q_lam(1:ny),q_R(1:ny),q_new(1:ny),P_Th2(1:ny),eps_Th2(1:ny),T_Th2(1:ny),D_Th2(1:ny),H_Th2(1:ny)
    real*8 dTdy(1:ny),dTdeta(1:ny),d2Tdy2(1:ny),d2Tdeta2(1:ny),dTh2dy(1:ny),dTh2deta(1:ny),d2Th2dy2(1:ny),d2Th2deta2(1:ny)
    real*8 dnutdy(1:ny),dnutdeta(1:ny)

    call ddeta(ny,T,dTdeta,deta)
    dTdy = dTdeta*detady
    call ddeta(ny,Th2,dTh2deta,deta)
    dTh2dy = dTh2deta*detady
    call ddeta(ny,nut,dnutdeta,deta)
    dnutdy = dnutdeta*detady
    call d2deta2(ny,T,d2Tdeta2,deta)
    d2Tdy2 = d2Tdeta2*detady**2.d0 + dTdeta*d2etady2
    call d2deta2(ny,Th2,d2Th2deta2,deta)
    d2Th2dy2 = d2Th2deta2*detady**2.d0 + dTh2deta*d2etady2

    q_lam = ( lambda / Pr ) * dTdy
    q_R = ( nut / sigmaT ) * dTdy
    q_new = ( dlambdadT / Pr ) * dTh2dy

    P_Th2 = ( nut / sigmaT ) * dTdy**2.d0
    eps_Th2(2:ny-1) = ( eps(2:ny-1) / Kt(2:ny-1) ) * ( lambda(2:ny-1) / Pr ) * Th2(2:ny-1)
    eps_Th2(1) = eps_Th2(2)
    eps_Th2(ny) = eps_Th2(ny-1)

    T_Th2 = ( nut / sigmaTh2 )*d2Th2dy2 + ( dTh2dy / sigmaTh2 ) * dnutdy
    D_Th2 = ( lambda / Pr )*d2Th2dy2 + ( dTh2dy / Pr ) * dlambdadT * dTdy

    H_Th2 = ( 2.d0 / Pr ) * d2lambdadT2 * Th2 * dTdy**2.d0 + ( 2.d0 / Pr ) * dlambdadT * Th2 * d2Tdy2 &
    + ( 1.d0 / Pr ) * dlambdadT * dTh2dy * dTdy
    
    end


!!!*************************************************
!!!*						         	           *
!!!*               residuals                   *
!!!*								               *
!!!*************************************************

subroutine  residuals(ny,A,A0,resA)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: A(1:ny),A0(1:ny)
    real*8, INTENT(OUT) :: resA
    real*8 sumN, sumD
    integer j

    sumN = 0
    do j=1,ny
        sumN = sumN +dabs(A(j)- A0(j))
    enddo

    sumD = 0
    do j=1,ny
        sumD = sumD + dabs(A0(j))
    enddo

    resA = sumN/sumD

    end

!!!*************************************************
!!!*						         	             *
!!!*            grid                       *
!!!*								                 *
!!!*************************************************

subroutine grid(ny,dy_min,Re_tau,y,detady,d2etady2,deta)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(in) :: Re_tau,dy_min
    real*8, intent(out) :: y(1:ny),detady(1:ny),d2etady2(1:ny),deta
    integer j
    real*8 a, b, c, d, e, eta, y_max

    y_max = 2*Re_tau

    deta = y_max/(ny-1)

    a = (deta - dy_min)/(Re_tau*deta - deta*deta)
    b = (dy_min*Re_tau - deta*deta)/(Re_tau*deta - deta*deta)

    do j=1,ny/2
        eta = deta*(j-1)
        y(j) = a*eta**2.d0 + b*eta
        detady(j) = 1.d0/(2.d0*a*eta + b)
        d2etady2(j) = -2.d0*a/(2.d0*a*eta + b)**3.d0
    enddo

    c = -a
    d = 4.d0*Re_tau*a + b
    e = 2.d0*Re_tau*( 1.d0 - b -2.d0*a*Re_tau )

    do j=ny/2+1,ny
        eta = deta*(j-1)
        y(j) = c*eta**2.d0 + d*eta + e
        detady(j) = 1.d0/(2.d0*c*eta + d)
        d2etady2(j) = -2.d0*c/(2.d0*c*eta + d)**3.d0
    enddo

    print*, ' dy_max =', y(ny/2)-y(ny/2-1), ' dy_min =', y(2)-y(1), ' ratio =', (y(ny/2)-y(ny/2-1))/(y(2)-y(1))

    end

!!!*************************************************
!!!*						         	           *
!!!*            Thermal Diffusivity                *
!!!*								               *
!!!*************************************************
    
subroutine  thernal_diffusivity(lambda,dlambdadT,d2lambdadT2,T,Bp,Cp)
    implicit none
    real*8, intent(in) :: T,Bp,Cp
    real*8, intent(out) :: lambda,dlambdadT,d2lambdadT2
    integer j

    lambda = 1.d0 + Bp*T + Cp*T**2.d0
    dlambdadT = Bp + 2.d0*Cp*T
    d2lambdadT2 = 2.d0*Cp

    end

!!!*************************************************
!!!*						         	           *
!!!*            Thomas Algorithm                *
!!!*								               *
!!!*************************************************
    
subroutine  thomas_algorithm(var,ny,a_e,a_w,a_p,S)
    implicit none
    integer, intent(in) :: ny
    real*8, intent(inout) :: var(1:ny),a_e,a_w,a_p,S
    real*8 A(1:ny),C_apex(1:ny)
    integer j

    A(1) = 0.d0
    A(ny) = 0.d0
    C_apex(1) = var(1)
    C_apex(ny) = var(ny)

    do j=2,ny-1
        A(j) = a_e / ( a_p - a_w * A(j-1) )
        C_apex(j) = ( a_w * C_apex(j-1) + S ) / ( a_p - a_w * A(j-1) )
        var(j) = A(j) * var(j+1) + C_apex(j)
    enddo

    end

!!!*************************************************
!!!*						         	           *
!!!*            Correct Residuals               *
!!!*								               *
!!!*************************************************
    
subroutine  correct_residuals(alpha,res,res_old)
    implicit none
    real*8, intent(in) :: res,res_old
    real*8, intent(inout) :: alpha
    real*8 increment

    increment = 1.d-06

    if(res < res_old) then
        alpha = alpha * ( 1.d0 + increment )
    elseif(res > res_old) then
        alpha = alpha * ( 1.d0 + increment )
    else
        alpha = alpha
    endif

    alpha = min(1.d0,max(alpha,increment))

    end


