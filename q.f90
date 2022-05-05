!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!       Linear Interpolation to calculate a&b
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
subroutine calculate_ab(xa,xb,t)
    parameter(xa1=17.26939,xb1=35.86,xa2=21.87456,xb2=7.66)
    real xa,xb
    real,intent(in):: t 
    
    !print *, t
    xak=(xa2-xa1)/(-40+15)
    xbk=(xb2-xb1)/(-40+15)
    xam=xa1-xak*(-15)
    xbm=xb1-xbk*(-15)
    xa=xak*t+xam
    xb=xbk*t+xbm
    !print *, xa
    !print *, xb
    !pause
    if(t>-15) then
        xa=xa1
        xb=xb1
    else if(t<-40) then
        xa=xa2
        xb=xb2
    end if
end subroutine calculate_ab

!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!                   read micaps data
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
subroutine dataget(varname,t,u)
    parameter(ix=33,iy=18,iz=11,it=12)
    real t(ix,iy,iz,it),dtd(ix,iy,iz,it),u(ix,iy,iz,it)
    integer hpa(iz),ddhh0(it)
    character(4) ddhh1,p1
    character(3) p
    data hpa/1000,925,850,700,500,400,300,250,200,150,100/
    data ddhh0/1620,1708,1720,1808,1820,1908,1920,2008,2020,2108,2120,2208/
    character(len=*),intent(in):: varname

    do k=1,iz
        write(p,'(i3)') hpa(k)
        if (hpa(k).eq.1000) write(p1,'(i4)') hpa(k)
        do i=1,it
            write(ddhh1,'(i4)') ddhh0(i)
            if (hpa(k).eq.1000) then
                open(1,file="../micaps/"//varname//"/"//p1//"/0904"//ddhh1//".000",&
                    status="old")
            else
                open(1,file="../micaps/"//varname//"/"//p//"/0904"//ddhh1//".000",&
                    status="old")
            end if
            iii=4
            if(varname.eq."uv") iii=3
            do ii=1,iii
                read(1,*)
            end do
            read(1,*) t(:,:,k,i)
            if(varname.eq."uv") read(1,*) u(:,:,k,i)
            close(1)
        end do
    end do
end subroutine dataget

!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!           Generate longitudes and latitudes 
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
subroutine getlatlon(lat,latb,dlat,iy)
    parameter(pi=3.14159)
    integer,intent(in):: iy
    real,intent(out):: lat(iy)
    real,intent(in):: latb,dlat

    !print *, iy,latb,dlat
    do i=1,iy
        lat(i)=(latb+dlat*(i-1))*pi/180
    end do
    print *, lat
end subroutine getlatlon

program main
parameter(ix=33,iy=18,iz=11,it=12)
parameter(pi=3.14159,g=9.8,a=6371e3)
real t(ix,iy,iz,it),dtd(ix,iy,iz,it),td(ix,iy,iz,it),e(ix,iy,iz,it),&
    q(ix,iy,iz,it),u(ix,iy,iz,it),v(ix,iy,iz,it),uq(ix,iy,iz,it),vq(ix,iy,iz,it),&
    vqdiv(ix-2,iy-2,iz,it),adv(ix-2,iy-2,iz,it),&
    lat0(iy),lat(iy),lon(ix)
real xa,xb
integer hpa(iz)
data hpa/1000,925,850,700,500,400,300,250,200,150,100/

call dataget("temper",t,u)
call dataget("t-td",dtd,v)
do l=1,it
    do k=1,iz
        do j=1,iy
            do i=1,ix
                call calculate_ab(xa,xb,t(i,j,k,l))
                td(i,j,k,l)=t(i,j,k,l)-dtd(i,j,k,l)
                e(i,j,k,l)=6.1078*exp(xa*td(i,j,k,l)/(273.16+td(i,j,k,l)-xb))
                q(i,j,k,l)=622*e(i,j,k,l)/(hpa(k)-0.378*e(i,j,k,l))
            end do
        end do
    end do
end do

open(1234,file="q.grd",form="unformatted",access="direct",status="replace",recl=ix*iy*iz*4)
do i=1,it
    write(1234,rec=i) q(:,:,:,i)
end do
close(1234)

call dataget("uv",u,v)

do l=1,it
    do k=1,iz
        do j=1,iy
            do i=1,ix
                uq(i,j,k,l)=u(i,j,k,l)*q(i,j,k,l)/g
                vq(i,j,k,l)=v(i,j,k,l)*q(i,j,k,l)/g
            end do
        end do
    end do
end do

open(111,file="uvq.grd",form="unformatted",access="direct",&
    status="replace",recl=ix*iy*iz*4)
irec=0
do i=1,it
    irec=irec+1
    write(111,rec=irec) uq(:,:,:,i)
end do
do i=1,it
    irec=irec+1
    write(111,rec=irec) vq(:,:,:,i)
end do
close(111)
print *, "uvq"
pause
call getlatlon(lat0,12.,4.,iy)
call getlatlon(lon,32.,4.,ix)

do i=1,iy
    lat(iy-i+1)=lat0(i)
end do
print *, lat

dtheta=4.*pi/180
do l=1,it
    do k=1,iz
        do j=2,iy
            dx=a*cos(lat(j))*dtheta
            dy=a*dtheta
            do i=2,ix
                vqdiv(i-1,j-1,k,l)=0.5*(uq(i+1,j,k,l)-uq(i-1,j,k,l))/dx+&
                    0.5*(vq(i+1,j,k,l)-vq(i-1,j,k,l))/dy
                adv(i-1,j-1,k,l)=0.5*u(i,j,k,l)*(q(i+1,j,k,l)-q(i-1,j,k,l))/dx+&
                    0.5*v(i,j,k,l)*(q(i,j+1,k,l)-q(i,j-1,k,l))/dy
            end do
        end do
    end do
end do
print *, adv 
pause
print *, vqdiv
pause

open(123,file="vqdiv.grd",form="unformatted",access="direct",&
    status="replace",recl=(ix-2)*(iy-2)*iz*4)
open(222,file="qadv.grd",form="unformatted",access="direct",&
    status="replace",recl=(ix-2)*(iy-2)*iz*4)
do i=1,it
    write(123,rec=i) vqdiv(:,:,:,i)
    write(222,rec=i) adv(:,:,:,i)
end do
close(123)
close(222)
end
