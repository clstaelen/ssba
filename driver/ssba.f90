Program ssba
  ! Staelen & Huré, A&A, 2024 "Approaching the structure of rotating bodies from dimension reduction"
  ! compile and run with
  ! gfortran -O3 -Wall -unused ssba.f90; ./a.out
  Implicit None
  Integer,Parameter::AP=Kind(1.00D+00)
  Real(kind=AP),Parameter::pi=Atan(1._AP)*4,threshold=1.E-14_AP
  Integer::i,j,iter,n,itermax
  Real(kind=AP)::error,ebar_surf,npoly,int1,int2,C1old,mass,inemom,angmom,Om2
  Real(kind=AP),Allocatable::ebar(:),ell2(:),a(:),c(:,:),rho(:),hent(:),dhent(:),de2da(:),kernel(:),omega2(:)
  Logical::notfinished
  ! statements
  ! input
  n=2**8
  Allocate(ebar(0:n),ell2(0:n),a(0:n),c(0:n,0:n),rho(0:n),hent(0:n),dhent(0:n),de2da(0:n),kernel(0:n),omega2(0:n))
  itermax=100
  ! initialisation and starting guess (linear profiles)
  ebar_surf=0.9_AP
  npoly=3._AP
  print '("SSB-approximation (Staelen & Huré, A&A, 2024) - see Tab. 1")'
  print '("Nodes, surface axis ratio, poly. index",1X,i4,2(1X,1PE13.6))',n,ebar_surf,npoly
  a(:)=1._AP*(/(i,i=0,n)/)/n
  ebar(n)=ebar_surf
  ebar(0:n-1)=1._AP-(1._AP-ebar_surf)*a(0:n-1)
  ell2(:)=1._AP-ebar(:)**2
  hent(0)=1._AP
  hent(1)=0._AP
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
     C1old=hent(0)
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
     ebar(:)=Sqrt(1._AP-ell2(:))
     ! confocal parameters
     Do i=0,n
        c(0:n,i)=(/(a(j)**2*ell2(j)/a(i)**2-ell2(i),j=0,n)/)
     End Do
     ! computes the enthalpy H(a)
     dhent=hent
     Do i=0,n
        kernel(0:n)=(/(eta(j,n)-eta(j,i),j=0,n)/)*pi*2
        Call trapezoidal(n+1,rho(0:n),kernel(0:n),hent(i))
     End Do
     ! computes the mass density rho(a)
     rho(1:n-1)=(hent(1:n-1)/hent(0))**npoly
     dhent=dhent-hent
     error=Sum(Abs(dhent(0:n)))/(n+1)
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
     Write(1,'(i4,6(1X,1PE13.6))') i,a(i),de2da(i),ell2(i),rho(i),omega2(i),dhent(i)
  End Do
  Print '("End. See output file ssba.dat")'
  Close(1)
  Deallocate(ebar,ell2,a,c,rho,hent,dhent,de2da,kernel,omega2)

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

  Function kappa(j,i)
    Implicit None
    Real(kind=AP)::kappa
    Integer :: i,j
    If (j<i) Then
       kappa=a(j)/a(i)/ebar(i)*4/3*(1._AP-ebar(i))
       If (ebar(j)/=1._AP)  kappa=a(j)*ebar(j)/a(i)/ell2(j)*((1._AP-a(j)**2*ell2(j)/a(i)**2*2)&
            &*ML_A(Sqrt(1._AP-a(j)**2*ell2(j)/a(i)**2))&
            &+Sqrt(1.+c(j,i))*ML_A(Sqrt(1._AP-a(j)**2*ell2(j)/a(i)**2/(1.+c(j,i))))*2&
            &-ebar(i)*2-Sqrt(1._AP-a(j)**2*ell2(j)/a(i)**2))
    Else
       kappa=1._AP-(3._AP-ell2(i)*2)/3 
       If (ebar(j)/=1._AP) kappa=1._AP+(3._AP-ell2(i)*2)/ell2(j)*(ebar(j)*ML_A(ebar(j))-1._AP)
    End If
  End Function kappa

  Function chi(j,i)
    Implicit None
    Real(kind=AP)::chi
    Integer :: i,j 
    chi=0._AP
    If (j<i) Then
       chi=a(j)**3/a(i)**4*(ebar(i)-1._AP)/ebar(i)
       If (ebar(j)/=1._AP) chi=a(j)**3*ebar(j)/a(i)**4*(ML_A(Sqrt(1._AP-a(j)**2*ell2(j)/a(i)**2))&
            &-ML_A(Sqrt(1._AP-a(j)**2*ell2(j)/a(i)**2/(1._AP+c(j,i))))/Sqrt(1._AP+c(j,i)))
    End If
  End Function chi

  Function mu(j,i)
    Implicit None
    Real(kind=AP)::mu
    Integer :: i,j
    If (j<i) Then 
       mu=-a(j)**3/a(i)**3/3/ebar(i)**3
       If (ebar(j)/=1._AP) mu=a(j)*ebar(j)/a(i)/ell2(j)&
            &*(ML_A(Sqrt(1.-a(j)**2*ell2(j)/a(i)**2/(1._AP+c(j,i))))/Sqrt(1._AP+c(j,i))-1._AP/ebar(i))
    Else
       mu=-1._AP/3
       If (ebar(j)/=1._AP) mu=(ebar(j)*ML_A(ebar(j))-1._AP)/ell2(j)
    End If
  End Function mu

  Function eta(j,i)
    Implicit None
    Real(kind=AP)::eta
    Integer :: i,j 
    If (j<i) Then
       eta=a(j)**3/a(i)/ebar(i)*2/3
       If (ebar(j)/=1._AP) eta=-a(j)*ebar(j)/ell2(j)*(a(i)*ebar(i)&
            &-Sqrt(a(i)**2*ebar(i)**2+a(j)**2*ell2(j))&
            &*ML_A(Sqrt(1._AP-a(j)**2*ell2(j)/(a(i)**2*ebar(i)**2+a(j)**2*ell2(j)))))
    Else
       eta=a(j)**2-a(i)**2*ebar(i)**2/3
       If (ebar(j)/=1._AP) eta = (a(i)**2*ebar(i)**2+a(j)**2*ell2(j))*ebar(j)/ell2(j)&
            &*ML_A(ebar(j))-a(i)**2*ebar(i)**2/ell2(j)
    End If
  End Function eta

End Program ssba
