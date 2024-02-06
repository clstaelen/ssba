Program ssba
  ! Staelen & Huré, A&A, 2024 "Approaching the structure of rotating bodies from dimension reduction"
  ! compile and run with
  ! gfortran -Wall -unused ssba.f90; ./a.out
  Implicit None
  Integer,Parameter::AP=Kind(1.00D+00)
  Real(kind=AP),Parameter::pi=Atan(1._AP)*4,threshold=1.E-14_AP
  Integer::i,j,iter,n,itermax
  Real(kind=AP)::error,ebar_surf,npoly,int1,int2,mass,inemom,angmom,Om2
  Real(kind=AP),Allocatable::ebar(:),ell2(:),a(:),rho(:),hent(:),drho(:),debar(:),de2da(:),kernel(:),omega2(:)
  Logical::notfinished
  ! statements
  ! input
  n=2**8
  Allocate(ebar(0:n),ell2(0:n),a(0:n),rho(0:n),hent(0:n),drho(0:n),debar(0:n),de2da(0:n),kernel(0:n),omega2(0:n))
  itermax=100
  ! initialisation and starting guess (quadratic profiles)
  ebar_surf=0.9_AP
  npoly=3._AP
  print '("SSB-approximation (Staelen & Huré, A&A, 2024) - see Tab. 1")'
  print '("Nodes, surface axis ratio, poly. index",1X,i4,2(1X,1PE13.6))',n,ebar_surf,npoly
  a(:)=1._AP*(/(i,i=0,n)/)/n
  ebar(n)=ebar_surf
  ebar(0:n-1)=1._AP-(1._AP-ebar_surf)*a(0:n-1)**2
  ell2(:)=1._AP-ebar(:)**2
  hent(0)=1._AP
  hent(n)=0._AP
  hent(1:n-1)=1._AP-a(1:n-1)**2
  rho(0)=1._AP
  rho(n)=0._AP
  rho(1:n-1)=hent(1:n-1)**npoly
  de2da(0)=0._AP
  iter=0
  notfinished=.True.
  print '("SCF-conv. threshold",1X,1PE13.6)',threshold
  ! SCF-loop
  Do While (notfinished)
     iter=iter+1
     debar(0:n) = ebar(0:n)
     drho(0:n)=rho(0:n)
     ! computes de2/da (a>0)
     Do i=1,n
        kernel(0:n)=(/(mu(j,i),j=0,n)/)
        Call trapezoidal(n+1,rho(0:n),kernel(0:n),int1)
        kernel(0:n)=(/(chi(j,i),j=0,n)/)*2
        Call trapezoidal(n+1,rho(0:n),kernel(0:n),int2)
        de2da(i)=int2/int1
     End Do
     ! computes e2(a) and ebar(a)
     Do i=0,n-1 
        Call trapezoidal(n-i+1,a(i:n),de2da(i:n),int1)
        ell2(i)=1._AP-ebar_surf**2-int1
     End Do
     ebar(0:n)=Sqrt(1._AP-ell2(0:n))
     debar(0:n) = abs(debar(0:n)-ebar(0:n))
     ! computes the enthalpy H(a)
     Do i=0,n
        kernel(0:n)=(/(eta(j,n)-eta(j,i),j=0,n)/)*pi*2
        Call trapezoidal(n+1,rho(0:n),kernel(0:n),hent(i))
     End Do
     ! computes the mass density rho(a)
     rho(1:n-1)=(hent(1:n-1)/hent(0))**npoly
     drho(0:n)=abs(drho(0:n)-rho(0:n))
     error=max(maxval(debar),maxval(drho))
     ! Print '(i4,1X,1PE13.6)',iter,error
     If ((error<threshold).Or.(iter>itermax)) notfinished=.False.
  End Do
  Print ('("Converged ? ",l)'),error<threshold
  Print ('("Iterations, error",1X,i4,2(1X,1PE13.6))'),iter,error
  ! computes the rotation rate squared, omega2(a)
  Do i=0,n 
     kernel(0:n)=(/(kappa(j,i),j=0,n)/)
     Call trapezoidal(n+1,rho(0:n),kernel(0:n),int1)
     omega2(i)=-int1*pi*2
  End Do
  ! mass
  kernel(:)=a(:)**2*rho(:)/ebar(:)*(1._AP-ell2(:)-a(:)/6*de2da(:))
  Call trapezoidal(n+1,a(0:n),kernel(0:n),mass)
  mass=mass*pi*4
  ! moment of inertia
  kernel(:)=a(:)**4*rho(:)/ebar(:)*(1._AP-ell2(:)-a(:)/10*de2da(:))
  Call trapezoidal(n+1,a(0:n),kernel(0:n),inemom)
  inemom=inemom*pi*8/3
  ! angular momentum
  kernel(:)=rho(:)*a(:)**4/ebar(:)*(1._AP-ell2(:)-a(:)/10*de2da(:))*Sqrt(omega2(:))
  Call trapezoidal(n+1,a(0:n),kernel(0:n),angmom)
  angmom=angmom*pi*8/3
  Om2=(angmom/inemom)**2
  Print '("Moment of inertia I and ang. momentum J",2(1X,1PE13.6))',inemom,angmom
  Print '("Mass and mean rotation rate J/I squared",2(1X,1PE13.6))',mass,Om2
  Open(unit=1,file="ssba.dat")
  Do i=0,n 
     Write(1,'(i4,6(1X,1PE13.6))') i,a(i),de2da(i),ell2(i),rho(i),omega2(i),drho(i)
  End Do
  Print '("End. See output file ssba.dat")'
  Close(1)
  Deallocate(ebar,ell2,a,rho,hent,drho,debar,de2da,kernel,omega2)

Contains

  Subroutine trapezoidal(len,x,y,ans) 
    Implicit None
    Integer,Intent(in)::len
    Integer :: k
    Real(kind=AP),Dimension(1:len),Intent(in)::x,y 
    Real(kind=AP),Intent(out)::ans 
    ans=0._AP
    Do k=2,len
       ans=ans+(x(k)-x(k-1))*(y(k)+y(k-1))/2
    End Do
  End Subroutine trapezoidal

  Function ML_A(x)
    Implicit None
    Real(kind=AP)::ML_A,x
    ML_A=1._AP
    If (x<1._AP) ML_A=Asin(Sqrt(1._AP-x**2))/Sqrt(1._AP-x**2)
    If (x>1._AP) ML_A=Log(x+Sqrt(x**2-1._AP))/Sqrt(x**2-1._AP)
  End Function ML_A

  function kappa(j,i)
    implicit none ; real(kind=AP)::kappa,cji ; integer :: i,j
    if (j<i) then
      kappa = a(j)/a(i)/ebar(i)*4/3*(1._AP-ebar(i))
      if (ebar(j)/=1._AP) then
        cji = a(j)**2*ell2(j)/a(i)**2-ell2(i)
        kappa = a(j)*ebar(j)/a(i)/ell2(j)*((1._AP-a(j)**2*ell2(j)/a(i)**2*2)&
          &*ML_A(sqrt(1._AP-a(j)**2*ell2(j)/a(i)**2))&
          &+sqrt(1.+cji)*ML_A(sqrt(1._AP-a(j)**2*ell2(j)/a(i)**2/(1.+cji)))*2&
          &-ebar(i)*2-sqrt(1._AP-a(j)**2*ell2(j)/a(i)**2))
      end if
    else
      kappa = 1._AP-(3._AP-ell2(i)*2)/3 
      if (ebar(j)/=1._AP) then
        kappa = 1._AP+(3._AP-ell2(i)*2)/ell2(j)*(ebar(j)*ML_A(ebar(j))-1._AP)
      end if
    end if
  end function

  function chi(j,i)
    implicit none ; real(kind=AP)::chi,cji ; integer :: i,j
    chi = 0._AP
    if (j<i) then
      chi = a(j)**3/a(i)**4*(ebar(i)-1._AP)/ebar(i)
      if (ebar(j)/=1._AP) then
        cji = a(j)**2*ell2(j)/a(i)**2-ell2(i)
        chi = a(j)**3*ebar(j)/a(i)**4*(ML_A(sqrt(1._AP-a(j)**2*ell2(j)/a(i)**2))&
        &-ML_A(sqrt(1._AP-a(j)**2*ell2(j)/a(i)**2/(1._AP+cji)))/sqrt(1._AP+cji))
      end if
    end if
  end function

  function mu(j,i)
    implicit none ; real(kind=AP)::mu,cji ; integer :: i,j
    if (j<i) then 
      mu = -a(j)**3/a(i)**3/3/ebar(i)**3
      if (ebar(j)/=1._AP) then
        cji = a(j)**2*ell2(j)/a(i)**2-ell2(i)
        mu = a(j)*ebar(j)/a(i)/ell2(j)*&
          &(ML_A(sqrt(1.-a(j)**2*ell2(j)/a(i)**2/(1._AP+cji)))/sqrt(1._AP+cji)-1._AP/ebar(i))
    end if
    else
      mu = -1._AP/3
      if (ebar(j)/=1._AP) then
        mu = (ebar(j)*ML_A(ebar(j))-1._AP)/ell2(j)
    end if
    end if
  end function

  function eta(j,i)
    implicit none ; real(kind=AP)::eta ; integer :: i,j 
    if (j<i) then
      eta = a(j)**3/a(i)/ebar(i)*2/3
      if (ebar(j)/=1._AP) then
      eta=-a(j)*ebar(j)/ell2(j)*(a(i)*ebar(i)&
        &-sqrt(a(i)**2*ebar(i)**2+a(j)**2*ell2(j))&
        &*ML_A(sqrt(1._AP-a(j)**2*ell2(j)/(a(i)**2*ebar(i)**2+a(j)**2*ell2(j)))))
      end if
    else
      eta = a(j)**2-a(i)**2*ebar(i)**2/3
      if (ebar(j)/=1._AP) then
        eta = (a(i)**2*ebar(i)**2+a(j)**2*ell2(j))*ebar(j)/ell2(j)&
          &*ML_A(ebar(j))-a(i)**2*ebar(i)**2/ell2(j)
      end if
    end if
  end function
End Program ssba
