module precision_rounding_mod
  use constants_mod
  use sp_domain_str_mod
  use sp_constraints_str_mod
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::precision_thermostat_bound,precision_scale_bound,precision_state_ratio,precision_constraint_ratios
contains
  subroutine precision_thermostat_bound(domain,constraints,x,force,half,vel,oracle,bound)
    type(s_domain),intent(in)::domain
    type(s_constraints),intent(in)::constraints
    real(wip),intent(in)::x(:,:,:),force(:,:,:),half(:,:,:),vel(:,:,:)
    real(8),intent(out)::oracle(9),bound(9)
    real(8)::u,g2,total,dm,p(3,3),mag(3,3),q(3,3),dq(3,3),m,energy(3),err(3),r(3),f(3),term(3)
    real(8)::magnitude(9),gd
    integer::cell,atom,group,entry,n,ids(8),h,j
    u=0.5d0*real(epsilon(1.0_wip),8);g2=gamma(2,u)
    oracle=0;bound=0;magnitude=0
    do cell=1,size(x,3)
      do atom=1,constraints%No_HGr(cell)
        m=real(domain%mass(atom,cell),8)
        energy=m*real(half(:,atom,cell),8)**2
        oracle(1:3)=oracle(1:3)+energy;bound(1:3)=bound(1:3)+g2*abs(energy)
        magnitude(1:3)=magnitude(1:3)+abs(energy)
        energy=m*real(vel(:,atom,cell),8)**2
        oracle(4:6)=oracle(4:6)+energy;bound(4:6)=bound(4:6)+g2*abs(energy)
        magnitude(4:6)=magnitude(4:6)+abs(energy)
      end do
      do group=0,constraints%connect
        if(group==0)then
          n=domain%num_water(cell)
        else
          n=constraints%HGr_local(group,cell)
        end if
        do entry=1,n
          if(group==0)then
            j=3;ids(1:j)=domain%water_list(1:j,entry,cell)
          else
            j=group+1;ids(1:j)=constraints%HGr_bond_list(1:j,entry,group,cell)
          end if
          total=0;p=0;mag=0
          do h=1,j
            atom=ids(h);m=real(domain%mass(atom,cell),8);total=total+m
            term=m*real(half(:,atom,cell),8);p(:,1)=p(:,1)+term;mag(:,1)=mag(:,1)+abs(term)
            term=m*real(vel(:,atom,cell),8);p(:,2)=p(:,2)+term;mag(:,2)=mag(:,2)+abs(term)
            term=m*real(x(:,atom,cell),8);p(:,3)=p(:,3)+term;mag(:,3)=mag(:,3)+abs(term)
          end do
          dm=gamma(j,u)*total
          if(total<=dm.or.total<=0)error stop 'Invalid mass in thermostat rounding oracle'
          q=p/total
          dq=(gamma(2*j,u)*mag+abs(q)*dm)/(total-dm)
          dq=dq+u/(1-u)*(abs(q)+dq)
          do h=1,2
            energy=total*q(:,h)**2
            err=dm*(abs(q(:,h))+dq(:,h))**2+total*dq(:,h)*(2*abs(q(:,h))+dq(:,h)) &
                +g2*(total+dm)*(abs(q(:,h))+dq(:,h))**2
            oracle(3*h-2:3*h)=oracle(3*h-2:3*h)+energy
            bound(3*h-2:3*h)=bound(3*h-2:3*h)+err
            magnitude(3*h-2:3*h)=magnitude(3*h-2:3*h)+abs(energy)
          end do
          do h=1,j
            atom=ids(h);r=real(x(:,atom,cell),8)-q(:,3);f=real(force(:,atom,cell),8)
            oracle(7:9)=oracle(7:9)-r*f
            bound(7:9)=bound(7:9)+abs(f)*(dq(:,3)+g2*(abs(r)+dq(:,3)))
            magnitude(7:9)=magnitude(7:9)+abs(r*f)
          end do
        end do
      end do
    end do
    ! Cover DP evaluation and accumulation in both the oracle and producers.
    gd=gamma(32*sum(domain%num_atom(1:size(x,3)))+1024,0.5d0*epsilon(1.0d0))
    bound=bound+gd*(magnitude+bound+1d0)
  contains
    real(8) function gamma(n,u)
      integer,intent(in)::n
      real(8),intent(in)::u
      gamma=real(n,8)*u/(1-real(n,8)*u)
    end function
  end subroutine
  subroutine precision_scale_bound(domain,constraints,vel,factor,oracle,bound)
    type(s_domain),intent(in)::domain
    type(s_constraints),intent(in)::constraints
    real(wip),intent(in)::vel(:,:,:),factor(3)
    real(8),intent(out)::oracle(:,:,:),bound(:,:,:)
    real(8)::u,alpha,total,dm,p(3),mag(3),q(3),dq(3),m,g(3),dg(3),v(3),change(3),error(3)
    integer::cell,atom,group,entry,n,j,h,ids(8)
    u=0.5d0*real(epsilon(1.0_wip),8);alpha=u/(1-u)
    oracle=0;bound=0
    g=real(factor,8)-1d0;dg=alpha*abs(g)
    do cell=1,size(vel,3)
      do atom=1,constraints%No_HGr(cell)
        oracle(:,atom,cell)=real(vel(:,atom,cell),8)*real(factor,8)
        bound(:,atom,cell)=alpha*abs(oracle(:,atom,cell))
      end do
      do group=0,constraints%connect
        if(group==0)then
          n=domain%num_water(cell)
        else
          n=constraints%HGr_local(group,cell)
        end if
        do entry=1,n
          if(group==0)then
            j=3;ids(1:j)=domain%water_list(1:j,entry,cell)
          else
            j=group+1;ids(1:j)=constraints%HGr_bond_list(1:j,entry,group,cell)
          end if
          total=0;p=0;mag=0
          do h=1,j
            atom=ids(h);m=real(domain%mass(atom,cell),8);total=total+m
            v=m*real(vel(:,atom,cell),8);p=p+v;mag=mag+abs(v)
          end do
          dm=real(j,8)*u/(1-real(j,8)*u)*total
          if(total<=dm.or.total<=0)error stop 'Invalid mass in scaling rounding oracle'
          q=p/total
          dq=(real(2*j,8)*u/(1-real(2*j,8)*u)*mag+abs(q)*dm)/(total-dm)
          dq=dq+alpha*(abs(q)+dq)
          change=g*q
          error=abs(g)*dq+(abs(q)+dq)*dg+alpha*(abs(g)+dg)*(abs(q)+dq)
          do h=1,j
            atom=ids(h);v=real(vel(:,atom,cell),8)+change
            oracle(:,atom,cell)=v
            bound(:,atom,cell)=error+alpha*(abs(v)+error)
            bound(:,atom,cell)=bound(:,atom,cell)+ &
              64d0*epsilon(1d0)*(abs(real(vel(:,atom,cell),8))+abs(change)+1d0)
          end do
        end do
      end do
    end do
  end subroutine
  subroutine precision_state_ratio(domain,actual,reference,dt,position,ratio,worst,av,rv)
    type(s_domain),intent(in)::domain
    real(wip),intent(in)::actual(:,:,:,:),reference(:,:,:,:)
    real(8),intent(in)::dt
    logical,intent(in)::position
    real(8),intent(out)::ratio,av,rv
    integer,intent(out)::worst(4)
    real(8)::u,atol,rtol,value,expected,limit,coordinate_scale,vector_scale,error
    integer::field,k,atom,cell
    if(dt<=0)error stop 'Invalid timestep for native state policy'
    u=0.5d0*real(epsilon(1.0_wip),8);ratio=0;worst=0;av=0;rv=0
    do field=1,6
      do cell=1,size(actual,3)
        do atom=1,domain%num_atom(cell)
          do k=1,3
            value=real(actual(k,atom,cell,field),8);expected=real(reference(k,atom,cell,field),8)
            if(.not.ieee_is_finite(value).or..not.ieee_is_finite(expected))error stop 'Nonfinite native state'
            atol=1d-10;rtol=1d-12
            if(field==5)then
              atol=1d-8;rtol=1d-10
            end if
            limit=atol+rtol*abs(expected)
            vector_scale=sqrt(sum(real(reference(:,atom,cell,field),8)**2))
            if(field==2)limit=limit+8*u*(1+vector_scale)
            if(field==1.or.(.not.position.and.(field==3.or.field==4))) &
              limit=limit+16*u*(1+vector_scale)
            if(position.and.field==1)then
              coordinate_scale=max(sqrt(sum(real(actual(:,atom,cell,2),8)**2)), &
                sqrt(sum(real(reference(:,atom,cell,2),8)**2)))
              limit=limit+16*u*(1+coordinate_scale)/dt
            end if
            error=abs(value-expected)/limit
            if(error>ratio)then
              ratio=error;worst=[field,k,atom,cell];av=value;rv=expected
            end if
          end do
        end do
      end do
    end do
  end subroutine

  subroutine precision_constraint_ratios(domain,constraints,array,ratios,raw)
    type(s_domain),intent(in)::domain
    type(s_constraints),intent(in)::constraints
    real(wip),intent(in)::array(:,:,:,:)
    real(8),intent(out)::ratios(2),raw(2)
    real(8)::u,alpha,tol,target
    integer::cell,group,entry,h,first,second,atom,other
    u=0.5d0*real(epsilon(1.0_wip),8);alpha=u/(1-u)
    tol=real(constraints%shake_tolerance,8);ratios=0;raw=0
    do cell=1,size(array,3)
      do group=1,constraints%connect
        do entry=1,constraints%HGr_local(group,cell)
          first=constraints%HGr_bond_list(1,entry,group,cell)
          do h=1,group
            second=constraints%HGr_bond_list(h+1,entry,group,cell)
            target=real(constraints%HGr_bond_dist(h+1,entry,group,cell),8)
            call pair(first,second,target,.false.)
          end do
        end do
      end do
      do entry=1,domain%num_water(cell)
        do first=1,2
          do second=first+1,3
            atom=domain%water_list(first,entry,cell);other=domain%water_list(second,entry,cell)
            target=real(constraints%water_rOH,8)
            if(first==2)target=real(constraints%water_rHH,8)
            call pair(atom,other,target,.true.)
          end do
        end do
      end do
    end do
  contains
    subroutine pair(a,b,target,water)
      integer,intent(in)::a,b
      real(8),intent(in)::target
      logical,intent(in)::water
      real(8)::x(3),y(3),va(3),vb(3),r(3),v(3),er(3),ev(3),error(2),limit(2),v0a(3),v0b(3)
      x=real(array(:,a,cell,2),8);y=real(array(:,b,cell,2),8)
      va=real(array(:,a,cell,4),8);vb=real(array(:,b,cell,4),8);r=x-y;v=va-vb
      er=alpha*(abs(x)+abs(y));ev=alpha*(abs(va)+abs(vb))
      if(water)then
        v0a=real(array(:,a,cell,3),8);v0b=real(array(:,b,cell,3),8)
        ev=ev+alpha*(abs(va-v0a)+abs(vb-v0b)+alpha*(abs(va)+abs(vb)))
      end if
      error=[abs(sum(r*r)-target*target),abs(sum(r*v))]
      limit(1)=2*tol*target+1d-12+2*sum(abs(r)*er)+sum(er*er)
      limit(2)=1.01d0*tol+1d-12+sum(abs(r)*ev+abs(v)*er+er*ev)
      if(.not.all(ieee_is_finite(error)).or..not.all(ieee_is_finite(limit))) &
        error stop 'Nonfinite native constraint comparison'
      ratios=max(ratios,error/limit);raw=max(raw,error)
    end subroutine
  end subroutine
end module
