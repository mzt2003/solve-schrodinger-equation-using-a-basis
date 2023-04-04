ccccccc      
      module parameter    !全局参数，也就是各个（物理）常数
            implicit none
            real*8,parameter :: hbarc=197.3269718d0   !hbar      ! NIST Ref 02.12.2014   ! MeV.fm           
            real*8,parameter :: finec=137.03599d0
            real*8,parameter :: amu=931.49432d0      !MeV
            real*8,parameter :: e2=1.43997d0         ! MeV.fm
            real*8,parameter :: PI=3.14159265        ! 圆周率pi
      end module
ccccccc
      module mesh              !mesh：积分的一些常量
            integer:: n_int    !积分格点数
            integer:: n_diff   !微分格点数
            integer ::n_basis  !基的个数
            real*8 :: hcm      !微分格点的步长
      end module
ccccccc
      module system   !系统的一些量，这里是一些np的两体信息
            implicit none
            real*8 :: z_1,z_2
            real*8 :: mass_1,mass_2
      end module
ccccccc
      module func      !一些用到的函数,把它们module在一个func里面
            use parameter
            contains
ccccccc
      real*8 function factorial(n) !阶乘函数，用实数表示整数，
            implicit none          !方便和实函数一起算，其输入为整型数
            integer::i,n           
            if(n==0) then
                  factorial=1
            else 
                  factorial=1
            do i=1,n
                  factorial=factorial*i
            end do
            end if
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc    
      real*8 function N_n(n)    !归一化系数，输入值为整数
            use parameter     !用到PI
            implicit none
            integer::n        !调用阶乘时，阶乘函数定义的输入是整型数
            real*8::alpha
            alpha=1.0
            N_n=sqrt(alpha/(sqrt(PI)*2**n*factorial(n))) 
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc 
      real*8 function Hermitian(n,z)!Hermitian多项式,n从0开始取   
            implicit none
            integer::k,n            !k是计数用的，n是输入的整数
            real*8::z               !z是输入的整数
            Hermitian=0             !初始化为0，后面要求和
            if(mod(n,2)==0) then    !偶数除以2余数为0
                  do k=0,n/2
                  Hermitian=Hermitian+(-1)**k*factorial(n)    
     &            *(2*z)**(n-2*k)/(factorial(k)*factorial(n-2*k))
                  end do
            else                    !else里面是奇数情况
                  do k=0,(n-1)/2
                  Hermitian=Hermitian+(-1)**k*factorial(n)*(
     &            2*z)**(n-2*k)/(factorial(k)*factorial(n-2*k))
                  end do
            end if
      end function       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!在调用的时候要注意   
!这几个函数的z(也就是径向距离)的定义是用real型定义的，
!所以要注意z不能输入1，2之类的整型数,
!而要写作2.0000之类的形式
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function HObasis(n,z)!HObasis:第n个能级的波函数
           implicit none        !n是整数，z是实数
           integer::n
           real*8::z
           HObasis=N_n(n)*Hermitian(n,z)*exp(-z*z/2)
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Gaussian 势，v0井深（使用时记得带负号），
!r0中心位置通常为0，a描述宽度。输出r处的势。
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
      function gausspot(r,v0,r0,a)
      implicit none
       real*8 r,v0,r0,gausspot,a
         if (a.gt.1e-6) then
           gausspot=V0*exp(-(r-r0)**2/a**2)
             else
               write(*,*)'a too small in gausspot!'
               stop
         endif
         return
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
! interpolation function for uniform grids 
!Y：求第Y个格点处的函数值，一般是小数,
!F：函数值数列, N：有N个格点
!该函数是直接copy的
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function FFR4(Y,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 F(N),P,P1,P2,Q,X,FFR4
      REAL*8 Y
      PARAMETER(X=.16666666666667)
      P=Y
      I=P
ccccccc
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFR4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFR4=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFR4=F(N)
      RETURN
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!高斯-勒让德积分。 N：积分格点数； x1,x2：积分下上限    
!输出x,w：格点位置和权重。
!gauleg相当于只对自变量区间作用，做一个重新离散化
!gauleg的核心是将这个新离散化的自变量数组给一组新的权重
!用这个新给的权重计算更稳定，这是重排的原因
!这个权重相当于dx(n)
!该函数是copy的
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gauleg(N,x1,x2,X,W)
        IMPLICIT NONE
        INTEGER N
        REAL*8 x1,x2,X(N),W(N)
        REAL*8 z1,z,xm,xl,pp,p3,p2,p1,pi,tol
        INTEGER m,i,j

        pi=acos(-1.0)
        tol=1.E-12

        m=(n+1)/2
        xm=0.5*(x2+x1)
        xl=0.5*(x2-x1)

        DO 10 i=1,m
         z=cos(pi*(i-0.25)/(N+0.5))

 20      CONTINUE
         p1=1.0E0
         p2=0.0E0
         DO 30 j=1,N
          p3=p2
          p2=p1
          p1=((2*j-1)*z*p2-(j-1)*p3)/j
 30      CONTINUE
         pp=N*(z*p1-p2)/(z*z-1.0E0)
         z1=z
         z=z1-p1/pp
         IF( abs(z1-z) .GT. tol) GOTO 20 ! Scheifenende

         X(i) = xm - xl*z
         X(n+1-i) = xm + xl*z
         W(i) = 2.E0*xl/((1.0-z*z)*pp*pp)
         W(n+1-i) = W(i)
 10     CONTINUE
        END SUBROUTINE gauleg 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 下面算二阶导数用了一个简单的三点差分
! 形式选取了subroutine的形式
! 其输入值依次为：函数（值）数组y，
! 以及计划输出的二阶导数数组d2y，
! 数组的size，n；以及步长dx（这里是均匀步长）
! 其中只有d2y是一个空的数组用来往里面填写并输出  
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine second_derivative(y,d2y,n,dx)
      implicit none
      integer::n,i
      real*8::dx
      real*8,dimension(0:n)::y
      real*8,dimension(0:n)::d2y
      do i=1,n-1
            d2y(i)=((y(i+1)-2*y(i)+y(i-1)))/(dx**2)
      end do
      d2y(0)=d2y(1)
      d2y(n)=d2y(n-1)
      end subroutine second_derivative 
      end module
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 因为要使用gauss-legend积分，这种积分方法比较稳定，其过程中需要
! 将积分变量区间重排，即重新给定权重，所以需要将具体积分算矩阵元
! 的过程中用到的数组重写好，这就需要一个subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! gauleg，FFR4，以及后续的求积分的函数或subroutine的逻辑是：
! 先讲要算的区间去用gauleg重新离散化并给定权重，然后将这些离散后
! 新的自变量点用FFR4插值计算出对应新的积分点上的函数值
! 而这些新的值就要存储在新的array里面
! 这就是后续subroutine的作用
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc   
      subroutine arrays_for_integrate(psi_1,d2psi_1,vpot_1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!这里所有的下标1均表示最后的输出量，而未加下标的均为过程量
!而这些过程量的使用就是暂存波函数以及势在不同位置处的值(size为n_diff+1的数组)
!在subroutine的括号里面加rr，rrw是因为subroutine的过程中
!call了gauleg，已经把rr和rrw存下来了，可以直接加 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
      use parameter     !使用相关的物理常数
      use system        !两体系统的信息
      use func          !程序中用到的很多函数
      use mesh          !格点信息
      implicit none
      real*8,dimension(0:n_diff)::vpot  !  势能，在微分分格点上
      real*8,dimension(0:n_diff)::psi    ! 某个波函数基
      real*8,dimension(0:n_diff)::d2psi ! psi的二阶导
      real*8,dimension(0:n_int)::rr,rrw ! 高斯积分需要用的变量
      real*8::v0        !势阱深
      integer::i,j      !循环变量
      real*8::r         !距离
      !积分节点上的函数值
      real*8,dimension(0:n_basis,0:n_int)::psi_1,d2psi_1
      real*8,dimension(0:n_int)::vpot_1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!下使用gauleg进行重新离散化并给定权重，其中： 
!真实积分的格点数为n_int+1，然后gauleg把0到hcm*n_diff的r区间去重新
!离散化，得到了新的位置数组rr和对应的权重rrw
!这里面积分上下限的范围一般比实际积分的范围要多给一些，也就是n_diff
!要更大一些，这样从中取的点及其权重更精确
!总结：用gauleg离散更大的自变量区间，再从中选取与实际积分格点数
!相同个数的自变量值(n_int+1个)去放到积分的数组里面
!在此之后就用FFR4去插值算这n_int+1个（从0取到n_int）
!新的横坐标对应的函数值（这里可以补充一下，因为FFR4用的也是n_diff+1
!大小的数组，而插值的精度显然是和格点数正相关的，所以这也是n_diff+1
!往大了取的原因）
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call gauleg(n_int+1,0.0d0,hcm*n_diff,rr,rrw)
      v0=-72.17      !任取势的深度，其值为负数
      do i=0,n_diff  !势的用来做FFR4的数组，一开始的size取n_diff+1
      r=i*hcm 
      vpot(i) = gausspot(r,v0,0.0d0,1.484d0)   
      end do 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!下面是两层循环，相当于i行j列的矩阵，第i行表示第i个基
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do i=0,n_basis
            do j=0,n_diff   
            psi(j)=HObasis(i,j*hcm)
            end do
            call second_derivative(psi,d2psi,n_diff,hcm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!call完之后二阶导数就是d2psi，其size还是n_diff+1
!下面可以直接用它做插值
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            do j=0, n_int
            psi_1(i,j)=FFR4(rr(j)/hcm, psi,n_diff+1)
            d2psi_1(i,j)=FFR4(rr(j)/hcm,d2psi,n_diff+1)
            end do 
      end do
      do i=0,n_int
            vpot_1(i)=FFR4(rr(i)/hcm,vpot,n_diff+1)
      end do
ccccccc   
    
      end subroutine arrays_for_integrate
ccccccc
      program nptest
            use parameter     !使用相关的物理常数
            use system        !将两体系统的信息包括在内
            use func          !程序中用到的很多函数
            use mesh          !格点信息
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !     
      !
      !
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      !HObasis及其二阶导数
      real*8,allocatable,dimension(:,:)::psi,d2psi 
      real*8,allocatable,dimension(:)::vpot   !势
      real*8,allocatable,dimension(:,:)::H      !Hamiltonian矩阵
      real*8,allocatable,dimension(:)::rr,rrw   !积分点和权重   
      real*8::mu        !约化质量，因为开头use了system
                        !所以电荷以及质量已经可以直接赋值了
      integer::i,j,k    !计数量，在积分的时候用
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*),Hermitian(3,2.0d0)!看一下Hermitian多项式对不对
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!下面给mesh参数以及系统参数赋值
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      n_int=40          !积分格点数目
      n_diff=300        !微分格点数目
      n_basis=8         !basis的数量
      hcm=0.01d0        !fm
      mass_1=1.0078     
      mass_2=1.0086     
      mu=amu*(mass_2*mass_1)/(mass_1+mass_2)    !约化质量
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!给参数赋值后就可以allocate具体大小的数组了
!如波函数 Hamiltonian矩阵
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!波函数及其二阶导数
      allocate(psi(0:n_basis,0:n_int),d2psi(0:n_basis,0:n_int))
!势
      allocate(vpot(0:n_int))
!Hamiltonian矩阵 
      allocate(H(0:n_basis,0:n_basis))
!积分点及其权重
      allocate(rr(0:n_int),rrw(0:n_int))
!用前面的subroutine来获得basis及其二阶导数，还有势
!以及积分格点和权重
!这就是具体积分时用到的所有量 
      call gauleg(n_int+1,0.0d0,hcm*n_diff,rr,rrw)
      call arrays_for_integrate(psi,d2psi,vpot)
!下进行积分算矩阵元
          do i=0,n_basis
              do j=0,n_basis
                H(i,j)=0
                do k=0,n_int
                    H(i,j)=H(i,j)+ psi(i,k)* 
     &               (d2psi(j,k)*(-hbarc*hbarc/2/mu)
     &               +vpot(k)*psi(j,k))*rrw(k)
                end do
              end do
          end do
      write(22,*) H 
!deallocate所有allocate的内存
      deallocate(psi,d2psi,vpot,H,rr,rrw)
      
      end program
ccccccc


    
    
    
