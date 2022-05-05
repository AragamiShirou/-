!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!                       Read .dat
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
subroutine dataget(datafile,hgt)
    parameter(ix=33,iy=18,iz=11,iit=12)
    real hgt(ix,iy,iz,iit)!,u(ix,iy,iz,nt)
    character(len=200) datafile

    open(1,file="./micaps/"//trim(datafile),form="unformatted",&
        access="direct",status="old",recl=ix*iy*iz*iit*4)
    read(1,rec=1) hgt 
    close(1)
end subroutine dataget

!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!           Generate longitudes and latitudes
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
subroutine getlatlon(lat,latb,dlat,iy)
    parameter(pi=3.14159)
    real lat(iy),latb,late
    integer iy
    do i=1,iy
        lat(i)=(latb+(i-1)*dlat)*pi/180
    end do
end subroutine getlatlon

program main
parameter(ix=33,iy=18,iz=11,it=12)!,itb=417,ite=453)
parameter(pi=3.14159,a=6371e3,gf=0.67197e5)
real hgt(ix,iy,iz,it),zetag(ix-2,iy-2,iz,it)
real dx,dy,lat(iy),lon(ix)
integer hpa(iz)
data hpa/1000,925,850,700,500,400,300,250,200,150,100/

call dataget("hgt.grd",hgt)!,itb,ite,it)
call getlatlon(lon,32.,4,ix)
call getlatlon(lat,-80.,-4,iy)

dtheta=4*pi/180
do l=1,it
    do k=1,iz
        do j=2,iy-1
            do i=2,ix-1
                dx=a*cos(lat(j))*dtheta
                dy=a*dtheta
                m=gf/sin(lat(j))
                ddhddx=(hgt(i+1,j,k,l)+hgt(i-1,j,k,l)-2*hgt(i,j,k,l))/(dx**2)
                ddhddy=(hgt(i,j+1,k,l)+hgt(i,j-1,k,l)-2*hgt(i,j,k,l))/(dy**2)
                dhdy=(hgt(i,j+1,k,l)-hgt(i,j-1,k,l))/(2*dy)
                zetag(i-1,j-1,k,l)=m*(ddhddx+ddhddy-dhdy*tan(lat(j))/a)
            end do
        end do
    end do
end do
print *, zetag
open(111,file="zetag.grd",form="unformatted",access="direct",&
    status="replace",recl=(ix-2)*(iy-2)*iz*4)
do i=1,it
    write(111,rec=i) zetag(:,:,:,i)
end do
close(111)
end

