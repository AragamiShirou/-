!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!                       Read .dat
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
subroutine dataget(varname,u)
parameter(ix=33,iy=18,iz=11,iit=12)
real hgt(ix,iy,iz,iit),u(ix,iy,iz,iit)
character(1) varname

open(1,file="./micaps/"//varname//".grd",form="unformatted",&
    access="direct",status="old",recl=ix*iy*iz*iit*4)
read(1,rec=1) u(:,:,:,:)
close(1)
!u(:,:,:,:)=hgt(:,:,:,itb:ite)
end subroutine dataget

!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!           Generate longitudes and latitudes
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
subroutine getlatlon(lat,latb,dlat,iy)
parameter(pi=3.14159)
real lat(iy),latb,dlat
integer iy
do i=1,iy
    lat(i)=(latb+(i-1)*dlat)*pi/180
end do
end subroutine getlatlon

program main
parameter(ix=33,iy=18,iz=11,it=12)!,itb=415,ite=453)
parameter(pi=3.14159,a=6.371e6)
real h(ix,iy,iz,it),u(ix,iy,iz,it),v(ix,iy,iz,it),div(ix-2,iy-2,iz,it),&
    zeta(ix-2,iy-2,iz,it),adv(ix-4,iy-4,iz,it),omega(ix-2,iy-2,iz,it),&
    omegap(ix-2,iy-2,iz,it)
real dx,dy,lat(iy),lon(ix)
integer hpa(iz)
data hpa/1000,925,850,700,500,400,300,250,200,150,100/!,70,50,30,20,10/

call dataget("u",u)!,itb,ite,it)
call dataget("v",v)!,itb,ite,it)

call getlatlon(lat,12.,4.,iy)
call getlatlon(lon,80.,-4.,ix)

dtheta=4.*pi/180
do l=1,it
    do k=1,iz
        do i=2,ix-1
            do j=2,iy-1
                dx=a*cos(lat(j))*dtheta
                dy=a*dtheta
                div(i-1,j-1,k,l)=(u(i+1,j,k,l)-u(i-1,j,k,l))/(2*dx)+&
                    (v(i,j+1,k,l)-v(i,j-1,k,l))/(2*dy)-v(i,j,k,l)*tan(lat(j))/a
                zeta(i-1,j-1,k,l)=(v(i+1,j,k,l)-v(i-1,j,k,l))/(2*dx)-&
                    (u(i,j+1,k,l)-u(i,j-1,k,l))/(2*dy)+u(i,j,k,l)*tan(lat(j))/a
            end do
        end do
    end do
end do
!print *, div(2,3,7,16)
!print *, zeta(2,3,9,16)

do l=1,it
    do k=1,iz
        do j=2,iy-1
            do i=2,ix-1
                dx=a*cos(lat(j))*dtheta
                dy=a*dtheta
                adv(i,j,k,l)=u(i,j,k,l)*(zeta(i+1,j,k,l)-zeta(i-1,j,k,l))/(2*dx)+&
                    v(i,j,k,l)*(zeta(i,j+1,k,l)-zeta(i,j-1,k,l))/(2*dy)
            end do
        end do
    end do
end do

print *, adv(:,:,:,6)


open(111,file="div.grd",form="unformatted",access="direct",&
    status="replace",recl=(ix-2)*(iy-2)*iz*4)
open(112,file="zeta.grd",form="unformatted",access="direct",&
    status="replace",recl=(ix-2)*(iy-2)*iz*4)
open(113,file="uvadv.grd",form="unformatted",access="direct",&
    status="replace",recl=(ix-4)*(iy-4)*iz*4)
do i=1,it
    write(111,rec=i) div(:,:,:,i)
    write(112,rec=i) zeta(:,:,:,i)
    write(113,rec=i) adv(:,:,:,i)
end do
close(111)
close(112)
close(113)

!8888888888888888888888888888888888888888888888888888888888
!                   omega
!8888888888888888888888888888888888888888888888888888888888
omega(:,:,1,:)=0
do l=1,it
    do k=2,iz
        dp=abs(hpa(k)-hpa(k-1))
        do j=1,iy-2
            do i=1,ix-2
                omega(i,j,k,l)=omega(i,j,k-1,l)+0.5*(div(i,j,k-1,l)+div(i,j,k,l))*dp
            end do
        end do
    end do
end do

mm=iz*(iz+1)*0.5
omegap(:,:,iz,:)=0
do l=1,it
    do k=1,iz-1
        do j=1,iy-2
            do i=1,ix-2
                omegap(i,j,k,l)=omega(i,j,k,l)-(k+1)*k*(omega(i,j,iz,l)-0)*0.5/mm   
            end do
        end do
    end do
end do


open(1111,file="omega.grd",access="direct",form="unformatted",&
    status="replace",recl=(ix-2)*(iy-2)*iz*4)
do i=1,it
    write(1111,rec=i) omegap(:,:,:,i)
end do
close(1111)
end
