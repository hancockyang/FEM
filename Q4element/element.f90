! module of elements
! class of element is
!  of features:
!   1. elem_no: element no
!   2. node_no: global node number of the element, dim of (4)
!   3. coordxy: x-y coords of the nodes, dim of (4,2)
!   4. coordrs: r-s coords of the nodes, dim of (4,2)
!   5. disp   : displacement of the nodes, dim of (4)
!
! class of element list
!  it is the box of pointers that point to class element
!  it also contains time_instance
        
      module element_mod
           implicit none
           ! Temperature Q4 element
           type Q4_temp
              ! elem_no   : element number
              ! node_no   : global node number, vector of length(4)
              ! coordxy   : global x-y coords (x,y), matrix of dim (4,2)
              ! coordrs   : global r-s coords (r,s), matrix of dim (4,2)
              ! Ke        : element stiffness matrix, matrix of dim (4,4)
              ! Me        : element mass matrix, matrix of dim (4,4)
              ! Ce        : element damping matrix, matrix of dim (4,4)
              ! temp      : temperature at four nodes, vector of length(4)
              ! Re        : element load vector, vector of length(4)
              ! a0,a1,a2  : material properties, scalar
              ! rho,c     : mass density, scalar
              ! f         : body force
              integer:: elem_no
              integer:: node_no(4)
              double precision,dimension(4,2)::coordxy,coordrs 
              double precision,dimension(4,4):: Ke, Me, Ce
              double precision,dimension(4)::temp,Re
              double precision::a0,a1,a2,rho,c,f
              contains
                 ! procedures
                 ! Ke_integrand   : integrand of Ke
                 ! Me_integrand   : integrand of Me
                 ! Ke_matrix      : calculate Ke matrix by numerical integration
                 ! Me_matrix      : calculate Me matrix by numerical integration
                 ! force_vector   : calculate body force vector
                 procedure:: initial
                 procedure:: Ke_integrand
                 procedure:: Me_integrand
                 procedure:: Ke_matrix
                 procedure:: Me_matrix
                 procedure:: force_vector
           end type Q4_temp
        ! Inheritance from Q4_temp class
        !  Q4_temp_bc_b: bottom boundary
           type, extends(Q4_temp)::Q4_temp_bc_b
              ! bc_node_no: boundary node number
              ! q1        : heat flux along bottom boundary
              ! q2        : heat flux along bottom boundary
              ! e12       : boundary force along edge12
              ! bc_temp   : boundary temperature
              integer::bc_node_no(2)
              double precision::q1,q2,e12(4)
              double precision,dimension(2)::bc_temp
              contains
                procedure:: edgforce_b
           end type Q4_temp_bc_b

        !  Q4_temp_bc_t: top boundary
           type, extends(Q4_temp)::Q4_temp_bc_t
              ! bc_node_no: boundary node number
              ! q3        : heat flux along top boundary
              ! q4        : heat flux along top boundary
              ! e34       : boundary force along edge34
              ! bc_temp   : boundary temperature
              integer::bc_node_no(2)
              double precision::q3,q4,e34(4)
              double precision,dimension(2)::bc_temp
              contains
                procedure:: edgforce_t
           end type Q4_temp_bc_t

           
        !  Q4_temp_bc_l: left boundary
           type, extends(Q4_temp)::Q4_temp_bc_l
              ! bc_node_no: boundary node number
              ! q1        : heat flux along left boundary
              ! q4        : heat flux along left boundary
              ! e41       : boundary force along edge41
              ! bc_temp   : boundary temperature
              integer::bc_node_no(2)
              double precision::q1,q4,e41(4)
              double precision,dimension(2)::bc_temp
              contains
                procedure:: edgforce_l
           end type Q4_temp_bc_l
       
        !  Q4_temp_bc_r: right boundary
           type, extends(Q4_temp)::Q4_temp_bc_r
              ! bc_node_no: boundary node number
              ! q2         : heat flux along top boundary
              ! q3         : heat flux along top boundary
              ! e23       : boundary force along edge23
              ! bc_temp   : boundary temperature
              integer::bc_node_no(2)
              double precision::q2,q3,e23(4)
              double precision,dimension(2)::bc_temp
              contains
                procedure:: edgforce_r
           end type Q4_temp_bc_r

           contains
              subroutine initial(Q4)
                implicit none
                class(Q4_temp),intent(inout)::Q4
                Q4%elem_no=0
                Q4%node_no=0
                Q4%coordrs=0.0d0
                Q4%coordxy=0.0d0
                Q4%Ke=0.0d0
                Q4%Me=0.0d0
                Q4%Ce=0.0d0
                Q4%temp=0.0d0
                Q4%Re=0.0d0
                Q4%a0=0.0d0
                Q4%a1=0.0d0
                Q4%a2=0.0d0
                Q4%rho=0.0d0
                Q4%c=0.0d0
                Q4%f=0.0d0

              end subroutine initial

              subroutine Ke_matrix(Q4)
                implicit none
                class(Q4_temp),intent(inout)::Q4
                double precision:: r, s
                double precision, dimension(4,4):: G
                double precision, dimension(2):: quad_point, weight
                integer::k,l
                quad_point = (/ 1.0d0/dsqrt(3.0d0),&
                       -1.0d0/dsqrt(3.0d0) /)
                weight  = (/ 1.0d0, 1.0d0 /)
                do k = 1,2
                   r = quad_point(k)
                   do l = 1,2
                      s = quad_point(l)
		              CALL Q4%Ke_integrand(r,s,G)
		              Q4%Ke = Q4%Ke + weight(l)*weight(k)*G
                   enddo   
                enddo
              end subroutine Ke_matrix
              
              subroutine Me_matrix(Q4)
                implicit none
                class(Q4_temp),intent(inout)::Q4
                double precision:: r, s
                double precision, dimension(4,4):: M
                double precision, dimension(2):: quad_point, weight
                integer::k,l
                quad_point = (/ 1.0d0/dsqrt(3.0d0),&
                       -1.0d0/dsqrt(3.0d0) /)
                weight  = (/ 1.0d0, 1.0d0 /)
                do k = 1,2
                   r = quad_point(k)
                   do l = 1,2
                      s = quad_point(l)
		              call Q4%Me_integrand(r,s,M)
                      Q4%Me = Q4%Me + weight(l)*weight(k)*M
                   enddo   
                enddo
              end subroutine Me_matrix
              
              subroutine Ke_integrand(Q4,r,s,G)
                implicit none
                class(Q4_temp),intent(in)::Q4
                double precision,intent(in)::r,s
                double precision:: xr,xs,yr,ys,detJ
                double precision,dimension(4):: Nr,Ns,shapefun
                double precision,dimension(2,2)::Jaccobi,Kappa,Jinv
                double precision,dimension(2,4)::Nrs,B
                double precision,dimension(4,4),intent(out)::G
                ! initialize 
                Kappa(1,1) = Q4%a1
                Kappa(1,2) = 0.0d0
                Kappa(2,1) = 0.0d0
                Kappa(2,2) = Q4%a2
        ! shape function:
                shapefun = .25d0*(/ (1.0d0 - r)*(1.0d0 - s),&
                        (1.0d0 + r)*(1.0d0 - s),&
                        (1.0d0 + s)*(1.0d0 + r),&
                        (1.0d0 - r)*(1.0d0 + s) /);
        ! derivative of shape function
                Nr = (/ -(1.0d0 - s)/4.0d0, (1.0d0 - s)/4.0d0,&
                        (1.0d0 + s)/4.0d0, ( -1.0d0 - s)/4.0d0 /)
                Ns = (/ -(1.0d0 - r)/4.0d0, (-1.0d0 - r)/4.0d0,&
                        (1.0d0 + r)/4.0d0, (1.0d0 - r)/4.0d0 /);
                Nrs(1,:) = Nr
                Nrs(2,:) = Ns
        ! xr, xs, yr, ys
                xr = DOT_PRODUCT(Nr,Q4%coordxy(:,1))
                xs = DOT_PRODUCT(Ns,Q4%coordxy(:,1))
                yr = DOT_PRODUCT(Nr,Q4%coordxy(:,2))
                ys = DOT_PRODUCT(Ns,Q4%coordxy(:,2))
        ! Jacobi         
                Jaccobi(1,1) = xr;
                Jaccobi(1,2) = yr;
                Jaccobi(2,1) = xs;
                Jaccobi(2,2) = ys;
        ! determinant of Jacobi and inverse of Jacobi
                detJ = xr*ys - xs*yr
                Jinv(1,1) =  ys;        
                Jinv(1,2) = -yr;
                Jinv(2,1) = -xs;
                Jinv(2,2) =  xr;
                Jinv = 1.0d0/detJ*Jinv
                B = MATMUL(Jinv,Nrs);
                G = MATMUL(MATMUL(TRANSPOSE(B),Kappa),B)*detJ &
                        + Q4%a0*spread(shapefun,dim = 2, ncopies = &
                        size(shapefun))*spread(shapefun,dim = 1,&
                        ncopies = size(shapefun))*detJ
              end subroutine Ke_integrand

              subroutine Me_integrand(Q4,r,s,M)
                implicit none
                class(Q4_temp),intent(in)::Q4
                double precision,intent(in)::r,s
                double precision:: xr,xs,yr,ys,detJ
                double precision,dimension(4):: Nr,Ns,shapefun
                double precision,dimension(2,4)::Nrs
                double precision,dimension(4,4),intent(out)::M
        ! shape function:
                shapefun = .25d0*(/ (1.0d0 - r)*(1.0d0 - s),&
                        (1.0d0 + r)*(1.0d0 - s),&
                        (1.0d0 + s)*(1.0d0 + r),&
                        (1.0d0 - r)*(1.0d0 + s) /);
        ! derivative of shape function
                Nr = (/ -(1.0d0 - s)/4.0d0, (1.0d0 - s)/4.0d0,&
                        (1.0d0 + s)/4.0d0, ( -1.0d0 - s)/4.0d0 /)
                Ns = (/ -(1.0d0 - r)/4.0d0, (-1.0d0 - r)/4.0d0,&
                        (1.0d0 + r)/4.0d0, (1.0d0 - r)/4.0d0 /);
                Nrs(1,:) = Nr
                Nrs(2,:) = Ns
        ! xr, xs, yr, ys
                xr = DOT_PRODUCT(Nr,Q4%coordxy(:,1))
                xs = DOT_PRODUCT(Ns,Q4%coordxy(:,1))
                yr = DOT_PRODUCT(Nr,Q4%coordxy(:,2))
                ys = DOT_PRODUCT(Ns,Q4%coordxy(:,2))
        ! determinant of Jacobi and inverse of Jacobi
                detJ = xr*ys - xs*yr
        ! mass integrant
                M = Q4%rho*Q4%c*spread(shapefun,dim = 2, ncopies = &
                        size(shapefun))*spread(shapefun,dim = 1,&
                        ncopies = size(shapefun))*detJ
              end subroutine Me_integrand

              subroutine force_vector(Q4)
                Implicit none
                class(Q4_temp),intent(inout)::Q4
                double precision:: f, f1, f2, f3, f4
                double precision, dimension(4):: x, y, rf
                f = Q4%f
                x = Q4%coordxy(:,1)
                y = Q4%coordxy(:,2)
                f1 = x(4)*(2.0d0*y(1) - y(2) - y(3)) + &
                        (2.0d0*x(1) - x(3))*(y(2) - y(4))+&
                        x(2)*(-2.0d0*y(1) + y(3) + y(4))
                f2 = x(4)*y(1) + 2.0d0*x(1)*y(2) - x(1)*y(3) - &
                        x(4)*y(3) + 2.0d0*x(2)*(y(3) - y(1)) - &
                        x(1)*y(4) + x(3)*(y(1) - 2.0d0*y(2) + y(4))
                f3 = x(4)*(y(1) + y(2) -2.0d0*y(3)) +  &
                        (x(1) - 2.0d0*x(3))*(y(2) - y(4)) - &
                        x(2)*(y(1) - 2.0d0*y(3) + y(4))
                f4 = 2.0d0*x(4)*(y(1) - y(3)) + x(2)*(y(3) - y(1)) - &
                        x(3)*(y(1) + y(2) - 2.0d0*y(4)) + &
                        x(1)*(y(2) + y(3) - 2.0d0*y(4))
                rf = (/ f1, f2, f3, f4/)
                Q4%Re = f*1.0d0/12.0d0*rf
              end

              subroutine edgforce_t(Q4_t)
                use dispmodule
                implicit none
                class(Q4_temp_bc_t),intent(inout)::Q4_t
                double precision:: L34,q3,q4
                double precision, dimension(4)::x,y,q
                double precision, dimension(4,4):: B
                q3 = Q4_t%q3
                q4 = Q4_t%q4
                q=(/0.0d0,0.0d0,q3,q4/)
                x = Q4_t%coordxy(:,1)
                y = Q4_t%coordxy(:,2)
                L34 = sqrt((x(4)-x(3))**2+(y(4)-y(3))**2)
                Q4_t%e34=0.0d0
                B(:,1)=(/0.d0,0.0d0,0.0d0,0.0d0/)
                B(:,2)=(/0.d0,0.0d0,0.0d0,0.0d0/)
                B(:,3)=(/0.d0,0.0d0,2.0d0/3.0d0,1.0d0/3.0d0/)
                B(:,4)=(/0.d0,0.0d0,1.0d0/3.0d0,2.0d0/3.0d0/)
                Q4_t%e34=MATMUL(B,q)*L34*0.5d0
              end subroutine edgforce_t

              subroutine edgforce_b(Q4_b)
                implicit none
                class(Q4_temp_bc_b),intent(inout)::Q4_b
                double precision:: L12,q1,q2
                double precision, dimension(4)::x,y,q
                double precision, dimension(4,4)::B
                q1=Q4_b%q1
                q2=Q4_b%q2
                q=(/q1,q2,0.0d0,0.0d0/)
                x = Q4_b%coordxy(:,1)
                y = Q4_b%coordxy(:,2)
                L12 = sqrt((x(2)-x(1))**2 + (y(2)-y(1))**2)
                Q4_b%e12=0.0d0
                B(:,1)=(/2.0d0/3.0d0,1.0d0/3.0d0,0.0d0,0.0d0/)
                B(:,2)=(/1.0d0/3.0d0,2.0d0/3.0d0,0.0d0,0.0d0/)
                B(:,3)=(/0.0d0,0.0d0,0.0d0,0.0d0/)
                B(:,4)=(/0.0d0,0.0d0,0.0d0,0.0d0/)
                Q4_b%e12=MATMUL(B,q)*L12*0.5d0
              end subroutine edgforce_b

              subroutine edgforce_l(Q4_l)
                implicit none
                class(Q4_temp_bc_l),intent(inout)::Q4_l
                double precision:: L41,q4,q1
                double precision,dimension(4)::x,y,q
                double precision,dimension(4,4)::B
                q1=Q4_l%q1
                q4=Q4_l%q4
                q=(/q1,0.0d0,0.0d0,q4/)
                x = Q4_l%coordxy(:,1)
                y = Q4_l%coordxy(:,2)
                L41 = sqrt((x(4)-x(1))**2 + (y(4)-y(1))**2)
                Q4_l%e41=0.0d0
                B(:,1)=(/2.0d0/3.0d0,0.0d0,0.0d0,1.0d0/3.0d0/)
                B(:,2)=(/0.0d0,0.0d0,0.0d0,0.0d00/)
                B(:,3)=(/0.0d0,0.0d0,0.0d0,0.0d00/)
                B(:,4)=(/1.0d0/3.0d0,0.0d0,0.0d0,2.0d0/3.0d0/)
                Q4_l%e41=MATMUL(B,q)*L41*0.5d0
              end subroutine edgforce_l

              subroutine edgforce_r(Q4_r)
                implicit none
                class(Q4_temp_bc_r),intent(inout)::Q4_r
                double precision:: q2,q3,L23
                double precision,dimension(4)::x,y,q
                double precision,dimension(4,4)::B
                q2=Q4_r%q2
                q3=Q4_r%q3
                q=(/0.0d0,q2,q3,0.0d0/)
                x = Q4_r%coordxy(:,1)
                y = Q4_r%coordxy(:,2)
                L23 = sqrt((x(3)-x(2))**2 + (y(3)-y(2))**2)
                Q4_r%e23=0.0d0
                B(:,1)=(/0.0d0,0.0d0,0.0d0,0.0d0/)
                B(:,2)=(/0.0d0,2.0d0/3.0d0,1.0d0/3.0d0,0.0d0/)
                B(:,3)=(/0.0d0,1.0d0/3.0d0,2.0d0/3.0d0,0.0d0/)
                B(:,4)=(/0.0d0,0.0d0,0.0d0,0.0d0/)
                Q4_r%e23=MATMUL(B,q)*L23*0.5d0
              end subroutine edgforce_r
           !type element_list
           !   type(element), dimension(:),pointer:: elements
           !   double precision:: time_instance
           !end type element_list
      end module element_mod
