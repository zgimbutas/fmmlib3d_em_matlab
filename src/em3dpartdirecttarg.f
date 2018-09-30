cc Copyright (C) 2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This is the end of the debugging code and the beginning 
c        of the Maxwell particle direct evaluation routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine em3dpartdirecttarg(zk,nsource,
     $     source,ifcjvec,cjvec,ifcmvec,cmvec,ifevec,evec,ifhvec,hvec,
     $     ntarget,target,ifevectarg,evectarg,ifhvectarg,hvectarg)
        implicit real *8 (a-h,o-z)
c
c       Maxwell interactions in R^3: evaluate all pairwise dipole
c       interactions (ignoring self-interaction) 
c       and interactions with targets via direct O(N^2) algorithm.
c
c       We use green3e and green3m for the Green's functions,
c       without the (1/4 pi) scaling.  Self-interactions are not-included.
c   
c       INPUT PARAMETERS:
c
c       zk: complex *16: frequency parameter
c       nsource: integer:  number of sources
c       source: real *8 (3,nsource):  source locations
c       ifcjvec:  electric dipole computation flag
c                  ifcjvec = 1   =>  include electric dipole contribution
c                                     otherwise do not
c       cjvec: complex *16 (3,nsource): electric dipole strengths
c       ifcmvec:  magnetic dipole computation flag
c                  ifcmvec = 1   =>  include magnetic dipole contribution
c                                     otherwise do not
c       cmvec: complex *16 (3,nsource): magnetic dipole strengths
c
c       ifevec:  E-field flag 
c                   (1=compute E-field, otherwise no)
c       ifhvec:  H-field flag 
c                   (1=compute H-field, otherwise no)
c
c       ntarget: integer:  number of targets
c       target: real *8 (3,ntarget):  target locations
c       ifevectarg:  target E-field flag 
c                   (1=compute E-field, otherwise no)
c       ifhvectarg:  target H-field flag 
c                   (1=compute H-field, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       evec: complex *16 (3,nsource): E field at source locations
c       hvec: complex *16 (3,nsource): H field at source locations
c       evectarg: complex *16 (3,ntarget): E field at target locations
c       hvectarg: complex *16 (3,ntarget): H field at target locations
c
        dimension source(3,1)
        complex *16 cjvec(3,1),cmvec(3,1)
        complex *16 evec(3,1),hvec(3,1)
        dimension target(3,1)
        complex *16 evectarg(3,1),hvectarg(3,1)
c
        complex *16 evectemp(3),hvectemp(3)
        complex *16 zk
c
        do i=1,nsource
        if( ifevec .eq. 1) then
           evec(1,i)=0
           evec(2,i)=0
           evec(3,i)=0
        endif
        if( ifhvec .eq. 1) then
           hvec(1,i)=0
           hvec(2,i)=0
           hvec(3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifevectarg .eq. 1) then
           evectarg(1,i)=0
           evectarg(2,i)=0
           evectarg(3,i)=0
        endif
        if( ifhvectarg .eq. 1) then
           hvectarg(1,i)=0
           hvectarg(2,i)=0
           hvectarg(3,i)=0
        endif
        enddo
c
        if( ifcjvec .eq. 0 .and. ifcmvec .eq. 0 ) return
c
c       ... sources
c
        if( ifevec .eq. 1 .or. ifhvec .eq. 1 ) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,evectemp,hvectemp) 
        do 6160 j=1,nsource
        do 6150 i=1,nsource

        if( i .eq. j ) goto 6150

        if( ifcjvec .eq. 1 .and. ifcmvec .eq. 1 ) then
        call dipole3emt(zk,source(1,i),source(1,j),
     $     cjvec(1,i),cmvec(1,i),evectemp,hvectemp)
        else
        if( ifcjvec .eq. 1 ) then
        call dipole3et(zk,source(1,i),source(1,j),
     $     cjvec(1,i),evectemp,hvectemp)
        endif
        if( ifcmvec .eq. 1 ) then
        call dipole3mt(zk,source(1,i),source(1,j),
     $     cmvec(1,i),evectemp,hvectemp)
        endif
        endif
c
        if( ifevec .eq. 1 ) then
        evec(1,j)=evec(1,j)+evectemp(1)
        evec(2,j)=evec(2,j)+evectemp(2)
        evec(3,j)=evec(3,j)+evectemp(3)
        endif
        if( ifhvec .eq. 1 ) then
        hvec(1,j)=hvec(1,j)+hvectemp(1)
        hvec(2,j)=hvec(2,j)+hvectemp(2)
        hvec(3,j)=hvec(3,j)+hvectemp(3)
        endif

 6150   continue
 6160   continue
C$OMP END PARALLEL DO
        endif
c
c       ... targets
c
        if( ifevectarg .eq. 1 .or. ifhvectarg .eq. 1 ) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,evectemp,hvectemp) 
        do j=1,ntarget
        do i=1,nsource

        if( ifcjvec .eq. 1 .and. ifcmvec .eq. 1 ) then
        call dipole3emt(zk,source(1,i),target(1,j),
     $     cjvec(1,i),cmvec(1,i),evectemp,hvectemp)
        else
        if( ifcjvec .eq. 1 ) then
        call dipole3et(zk,source(1,i),target(1,j),
     $     cjvec(1,i),evectemp,hvectemp)
        endif
        if( ifcmvec .eq. 1 ) then
        call dipole3mt(zk,source(1,i),target(1,j),
     $     cmvec(1,i),evectemp,hvectemp)
        endif
        endif
c
        if( ifevectarg .eq. 1 ) then
        evectarg(1,j)=evectarg(1,j)+evectemp(1)
        evectarg(2,j)=evectarg(2,j)+evectemp(2)
        evectarg(3,j)=evectarg(3,j)+evectemp(3)
        endif
        if( ifhvectarg .eq. 1 ) then
        hvectarg(1,j)=hvectarg(1,j)+hvectemp(1)
        hvectarg(2,j)=hvectarg(2,j)+hvectemp(2)
        hvectarg(3,j)=hvectarg(3,j)+hvectemp(3)
        endif

        enddo
        enddo
C$OMP END PARALLEL DO
        endif
c
        return
        end
c
c
c
c
c
        subroutine dipole3erhot(rk,source,target,cjvec,rho_e,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the electric pair of (cjvec, rho_e) located at the
c       arbitrary source location
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric current   
c       rho_e (complex *16) - the strength of the electric charge
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3),source(3),target(3)
        complex *16 cjvec(3),rho_e,evec(3),hvec(3),fout,rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
c       ... vector potential
c
        fout=exp(ima*rk*cd)/cd * (ima*rk)
c
        evec(1)=+fout*cjvec(1)
        evec(2)=+fout*cjvec(2)
        evec(3)=+fout*cjvec(3)
c
c       ... minus the gradient of scalar potential
c
        fout=rho_e*(1-ima*rk*cd)*exp(ima*rk*cd)/cd**3 
c
        evec(1)=evec(1)+dx*fout
        evec(2)=evec(2)+dy*fout
        evec(3)=evec(3)+dz*fout
c
c
        call green3m(rk,xyz,cjvec,hvec)
        hvec(1)=hvec(1) 
        hvec(2)=hvec(2) 
        hvec(3)=hvec(3) 
c
        return
        end
c
c
c
c
c
        subroutine dipole3mrhot(rk,source,target,cmvec,rho_m,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the magnetic pair of (cmvec, rho_m) located at the
c       arbitrary source location
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cmvec (complex *16) - the strength of the magnetic current   
c       rho_m (complex *16) - the strength of the magnetic charge
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3),source(3),target(3)
        complex *16 cmvec(3),rho_m,evec(3),hvec(3),fout,rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
c       ... vector potential
c
        fout=exp(ima*rk*cd)/cd * (ima*rk)
c
        hvec(1)=+fout*cmvec(1)
        hvec(2)=+fout*cmvec(2)
        hvec(3)=+fout*cmvec(3)
c
c       ... minus the gradient of scalar potential
c
        fout=rho_m*(1-ima*rk*cd)*exp(ima*rk*cd)/cd**3
c
        hvec(1)=hvec(1)+dx*fout
        hvec(2)=hvec(2)+dy*fout
        hvec(3)=hvec(3)+dz*fout
c
        hvec(1)=-hvec(1)
        hvec(2)=-hvec(2)
        hvec(3)=-hvec(3)
c
c
        call green3m(rk,xyz,cmvec,evec)
        evec(1)=evec(1) 
        evec(2)=evec(2) 
        evec(3)=evec(3) 
c
        return
        end
c
c
c
c
c
        subroutine dipole3emrhot(rk,source,target,
     $     cjvec,rho_e,cmvec,rho_m,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the electric pair of (cjvec, rho_e) and the magnetic pair of
c       (cmvec, rho_m) located at the arbitrary source location
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric current   
c       rho_e (complex *16) - the strength of the electric charge
c       cmvec (complex *16) - the strength of the magnetic current   
c       rho_m (complex *16) - the strength of the magnetic charge
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        complex *16 rk
        dimension xyz(3)
        complex *16 cjvec(3),cmvec(3)
        complex *16 rho_e,rho_m
        complex *16 evec(3),hvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
c
        call dipole3erhot(rk,source,target,cjvec,rho_e,evec1,hvec1)
        call dipole3mrhot(rk,source,target,cmvec,rho_m,evec2,hvec2)
c
        evec(1)=evec1(1)+evec2(1)
        evec(2)=evec1(2)+evec2(2)
        evec(3)=evec1(3)+evec2(3)
c
        hvec(1)=hvec1(1)+hvec2(1)
        hvec(2)=hvec1(2)+hvec2(2)
        hvec(3)=hvec1(3)+hvec2(3)
c
        return
        end
c
c
c
c
c
        subroutine em3dccpartdirecttarg(zk,nsource,
     $     source,ifcjvec,cjvec,rho_e,ifcmvec,cmvec,rho_m,
     $     ifevec,evec,ifhvec,hvec,
     $     ntarget,target,ifevectarg,evectarg,ifhvectarg,hvectarg)
        implicit real *8 (a-h,o-z)
c
c       Maxwell interactions in R^3: evaluate all pairwise dipole and charge
c       interactions (ignoring self-interaction) 
c       and interactions with targets via direct O(N^2) algorithm.
c
c       We use green3e and green3m for the Green's functions,
c       without the (1/4 pi) scaling.  Self-interactions are not-included.
c   
c       INPUT PARAMETERS:
c
c       zk: complex *16: frequency parameter
c       nsource: integer:  number of sources
c       source: real *8 (3,nsource):  source locations
c       ifcjvec:  electric dipole computation flag
c                  ifcjvec = 1   =>  include electric dipole contribution
c                                     otherwise do not
c       cjvec: complex *16 (3,nsource): electric dipole strengths
c       rho_e: complex *16 (nsource): electric charge strengths
c       ifmjvec:  magnetic dipole computation flag
c                  ifmjvec = 1   =>  include magnetic dipole contribution
c                                     otherwise do not
c       cmvec: complex *16 (3,nsource): magnetic dipole strengths
c       rho_m: complex *16 (nsource): magnetic charge strengths
c
c       ifevec:  E-field flag 
c                   (1=compute E-field, otherwise no)
c       ifhvec:  H-field flag 
c                   (1=compute H-field, otherwise no)
c
c       ntarget: integer:  number of targets
c       target: real *8 (3,ntarget):  target locations
c       ifevectarg:  target E-field flag 
c                   (1=compute E-field, otherwise no)
c       ifhvectarg:  target H-field flag 
c                   (1=compute H-field, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       evec: complex *16 (3,nsource): E field at source locations
c       hvec: complex *16 (3,nsource): H field at source locations
c       evectarg: complex *16 (3,ntarget): E field at target locations
c       hvectarg: complex *16 (3,ntarget): H field at target locations
c
        dimension source(3,1)
        complex *16 cjvec(3,1),cmvec(3,1)
        complex *16 rho_e(1),rho_m(1)
        complex *16 evec(3,1),hvec(3,1)
        dimension target(3,1)
        complex *16 evectarg(3,1),hvectarg(3,1)
c
        complex *16 evectemp(3),hvectemp(3)
        complex *16 zk
c
        do i=1,nsource
        if( ifevec .eq. 1) then
           evec(1,i)=0
           evec(2,i)=0
           evec(3,i)=0
        endif
        if( ifhvec .eq. 1) then
           hvec(1,i)=0
           hvec(2,i)=0
           hvec(3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifevectarg .eq. 1) then
           evectarg(1,i)=0
           evectarg(2,i)=0
           evectarg(3,i)=0
        endif
        if( ifhvectarg .eq. 1) then
           hvectarg(1,i)=0
           hvectarg(2,i)=0
           hvectarg(3,i)=0
        endif
        enddo
c
        if( ifcjvec .eq. 0 .and. ifcmvec .eq. 0 ) return
c
c       ... sources
c
        if( ifevec .eq. 1 .or. ifhvec .eq. 1 ) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,evectemp,hvectemp) 
        do 6160 j=1,nsource
        do 6150 i=1,nsource

        if( i .eq. j ) goto 6150

        if( ifcjvec .eq. 1 .and. ifcmvec .eq. 1 ) then
        call dipole3emrhot(zk,source(1,i),source(1,j),
     $     cjvec(1,i),rho_e(i),cmvec(1,i),rho_m(i),evectemp,hvectemp)
        else
        if( ifcjvec .eq. 1 ) then
        call dipole3erhot(zk,source(1,i),source(1,j),
     $     cjvec(1,i),rho_e(i),evectemp,hvectemp)
        endif
        if( ifcmvec .eq. 1 ) then
        call dipole3mrhot(zk,source(1,i),source(1,j),
     $     cmvec(1,i),rho_m(i),evectemp,hvectemp)
        endif
        endif
c
        if( ifevec .eq. 1 ) then
        evec(1,j)=evec(1,j)+evectemp(1)
        evec(2,j)=evec(2,j)+evectemp(2)
        evec(3,j)=evec(3,j)+evectemp(3)
        endif
        if( ifhvec .eq. 1 ) then
        hvec(1,j)=hvec(1,j)+hvectemp(1)
        hvec(2,j)=hvec(2,j)+hvectemp(2)
        hvec(3,j)=hvec(3,j)+hvectemp(3)
        endif

 6150   continue
 6160   continue
C$OMP END PARALLEL DO
        endif
c
c       ... targets
c
        if( ifevectarg .eq. 1 .or. ifhvectarg .eq. 1 ) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,evectemp,hvectemp) 
        do j=1,ntarget
        do i=1,nsource

        if( ifcjvec .eq. 1 .and. ifcmvec .eq. 1 ) then
        call dipole3emrhot(zk,source(1,i),target(1,j),
     $     cjvec(1,i),rho_e(i),cmvec(1,i),rho_m(i),evectemp,hvectemp)
        else
        if( ifcjvec .eq. 1 ) then
        call dipole3erhot(zk,source(1,i),target(1,j),
     $     cjvec(1,i),rho_e(i),evectemp,hvectemp)
        endif
        if( ifcmvec .eq. 1 ) then
        call dipole3mrhot(zk,source(1,i),target(1,j),
     $     cmvec(1,i),rho_m(i),evectemp,hvectemp)
        endif
        endif
c
        if( ifevectarg .eq. 1 ) then
        evectarg(1,j)=evectarg(1,j)+evectemp(1)
        evectarg(2,j)=evectarg(2,j)+evectemp(2)
        evectarg(3,j)=evectarg(3,j)+evectemp(3)
        endif
        if( ifhvectarg .eq. 1 ) then
        hvectarg(1,j)=hvectarg(1,j)+hvectemp(1)
        hvectarg(2,j)=hvectarg(2,j)+hvectemp(2)
        hvectarg(3,j)=hvectarg(3,j)+hvectemp(3)
        endif

        enddo
        enddo
C$OMP END PARALLEL DO
        endif
c
        return
        end
c
c
c
c
c
