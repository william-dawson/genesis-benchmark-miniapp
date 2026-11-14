#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
module genesis_thermo_mod
  use cudafor
  use constants_mod
  use sp_domain_str_mod
  use sp_constraints_str_mod
  use sp_constraints_mod,only:compute_kin_group,compute_virial_group,update_vel_group_3d
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::thermo_enabled,thermo_verify,thermo_statistics,thermo_output
  public::thermo_scale_reference,thermo_scale_check
  real(8),device,allocatable::total_force(:,:,:),half_velocity(:,:,:),partial(:,:)
  real(8),allocatable,pinned::host_partial(:,:),reference(:,:,:),actual(:,:,:)
contains
  logical function thermo_enabled()
    integer,save::selected=-1
    integer::status,code
    character(16)::value
    if(selected<0) then
      selected=0
      call get_environment_variable('GENESIS_THERMO_GPU',value,status=status)
      if(status==0) then
        read(value,*,iostat=code) selected
        if(code/=0.or.selected<0.or.selected>1) error stop 'GENESIS_THERMO_GPU must be 0 or 1'
      end if
    end if
    thermo_enabled=selected==1
  end function

  logical function thermo_verify()
    integer::status
    character(16)::value
    call get_environment_variable('GENESIS_THERMO_VERIFY',value,status=status)
    thermo_verify=status==0.and.trim(value)=='1'
  end function

  subroutine checked(status)
    integer,intent(in)::status
    if(status/=cudaSuccess) error stop 'GENESIS thermostat CUDA failure'
  end subroutine

  subroutine thermo_statistics(domain,constraints,v,x,fl,fs,m,im,na,ns,nw,wl,hl,li,half,values)
    type(s_domain),intent(in)::domain
    type(s_constraints),intent(in)::constraints
    real(8),device,intent(in)::v(:,:,:),x(:,:,:),fl(:,:,:),fs(:,:,:),m(:,:),im(:,:)
    integer,device,intent(in)::na(:),ns(:),nw(:),wl(:,:,:),hl(:,:),li(:,:,:,:)
    real(8),intent(in)::half
    real(8),intent(out)::values(9)
    integer::width,count,cell,atom,n
    real(8)::expected(9),virial(3,3),kin(3),ekin,ratio
    real(8),allocatable::host_x(:,:,:),host_force(:,:,:),host_half(:,:,:)
    width=size(v,2);count=size(v,3)
    if(.not.allocated(total_force)) then
      allocate(total_force(3,width,count),half_velocity(3,width,count))
      allocate(partial(9,count),host_partial(9,count))
    end if
    if(.not.allocated(reference)) allocate(reference(3,width,count),actual(3,width,count))
    if(any(shape(total_force)/=shape(v))) error stop 'GPU thermostat shape changed'
    call force_half<<<(width*count+255)/256,256>>>(width,count,na,im,fl,fs,half,total_force,half_velocity)
    call checked(cudaGetLastError())
    call group_statistics<<<count,128>>>(width,count,size(wl,1),size(wl,2),size(li,1),size(li,2), &
       size(hl,1),v,x,total_force,half_velocity,m,ns,nw,wl,hl,li,partial)
    call checked(cudaGetLastError())
    host_partial=partial;values=0
    do cell=1,count
      values=values+host_partial(:,cell)
    end do
    if(.not.all(ieee_is_finite(values))) error stop 'Nonfinite thermostat statistics'
    if(thermo_verify()) then
      allocate(host_x(3,width,count),host_force(3,width,count),host_half(3,width,count))
      actual=v;host_x=x;host_force=fl;host_half=fs
      host_force=host_force+host_half;reference=total_force;ratio=0
      do cell=1,count
        n=domain%num_atom(cell)
        ratio=max(ratio,maxval(abs(host_force(:,1:n,cell)-reference(:,1:n,cell))/ &
          (1.0d-8+1.0d-10*abs(host_force(:,1:n,cell)))))
        do atom=1,n
          host_half(:,atom,cell)=(half*domain%inv_mass(atom,cell))*host_force(:,atom,cell)
        end do
      end do
      reference=half_velocity
      do cell=1,count
        n=domain%num_atom(cell)
        if(.not.all(ieee_is_finite(reference(:,1:n,cell)))) error stop 'Nonfinite GPU thermostat half velocity'
        ratio=max(ratio,maxval(abs(host_half(:,1:n,cell)-reference(:,1:n,cell))/ &
          (1.0d-10+1.0d-12*abs(host_half(:,1:n,cell)))))
      end do
      if(ratio>1) error stop 'GPU thermostat force/half velocity differs from CPU'
      write(*,'(a,es20.12)') 'GENESIS_THERMO_VERIFY: force and half velocity pass ',ratio
      call compute_kin_group(constraints,count,domain%num_water,domain%water_list,domain%mass,host_half,kin,ekin)
      expected(1:3)=kin
      call compute_kin_group(constraints,count,domain%num_water,domain%water_list,domain%mass,actual,kin,ekin)
      expected(4:6)=kin;virial=0
      call compute_virial_group(constraints,count,domain%num_water,domain%water_list,domain%mass,host_x,host_force,virial)
      expected(7:9)=[virial(1,1),virial(2,2),virial(3,3)]
      ratio=maxval(abs(values-expected)/(1.0d-8+1.0d-10*abs(expected)))
      if(.not.all(ieee_is_finite(expected)).or.ratio>1) then
        print *, 'GENESIS_THERMO statistics actual/reference ',values,expected
        error stop 'GPU thermostat statistics differ from CPU'
      end if
      write(*,'(a,es20.12)') 'GENESIS_THERMO_VERIFY: statistics pass ',ratio
    end if
  end subroutine

  subroutine thermo_output(domain)
    type(s_domain),intent(inout)::domain
    integer::cell,n
    if(.not.allocated(total_force)) error stop 'GPU thermostat output before statistics'
    reference=total_force;actual=half_velocity
    do cell=1,size(reference,3)
      n=domain%num_atom(cell)
      domain%force(:,1:n,cell)=reference(:,1:n,cell)
      domain%velocity_half(:,1:n,cell)=actual(:,1:n,cell)
    end do
  end subroutine

  subroutine thermo_scale_reference(domain,constraints,v,factor)
    type(s_domain),intent(in)::domain
    type(s_constraints),intent(in)::constraints
    real(8),device,intent(in)::v(:,:,:)
    real(8),intent(in)::factor(3)
    if(.not.thermo_verify()) return
    if(.not.allocated(reference)) then
      allocate(reference,mold=domain%velocity)
      allocate(actual,mold=domain%velocity)
    end if
    reference=v
    call update_vel_group_3d(constraints,size(v,3),domain%num_water,domain%water_list,factor,domain%mass,reference)
  end subroutine

  subroutine thermo_scale_check(domain,v)
    type(s_domain),intent(in)::domain
    real(8),device,intent(in)::v(:,:,:)
    integer::cell,n
    real(8)::ratio
    if(.not.thermo_verify()) return
    actual=v;ratio=0
    do cell=1,size(v,3)
      n=domain%num_atom(cell)
      if(.not.all(ieee_is_finite(actual(:,1:n,cell)))) error stop 'Nonfinite GPU thermostat velocity'
      ratio=max(ratio,maxval(abs(actual(:,1:n,cell)-reference(:,1:n,cell))/ &
        (1.0d-10+1.0d-12*abs(reference(:,1:n,cell)))))
    end do
    if(ratio>1) error stop 'GPU thermostat group scaling differs from CPU'
    write(*,'(a,es20.12)') 'GENESIS_THERMO_VERIFY: scaling pass ',ratio
  end subroutine

  attributes(global) subroutine force_half(width,count,na,im,fl,fs,half,force,velocity)
    integer,value::width,count
    integer,device,intent(in)::na(count)
    real(8),device,intent(in)::im(width,count),fl(3,width,count),fs(3,width,count)
    real(8),value::half
    real(8),device,intent(out)::force(3,width,count),velocity(3,width,count)
    integer::index,cell,atom,k
    real(8)::factor
    index=(blockIdx%x-1)*blockDim%x+threadIdx%x
    cell=(index-1)/width+1;atom=mod(index-1,width)+1
    if(cell>count) return
    if(atom>na(cell)) return
    factor=half*im(atom,cell)
    do k=1,3
      force(k,atom,cell)=fl(k,atom,cell)+fs(k,atom,cell)
      velocity(k,atom,cell)=factor*force(k,atom,cell)
    end do
  end subroutine

  attributes(global) subroutine group_statistics(width,count,water_height,water_max,list_height,group_max,connect, &
      v,x,force,half,m,ns,nw,wl,hl,li,result)
    integer,value::width,count,water_height,water_max,list_height,group_max,connect
    real(8),device,intent(in)::v(3,width,count),x(3,width,count),force(3,width,count),half(3,width,count),m(width,count)
    integer,device,intent(in)::ns(count),nw(count),wl(water_height,water_max,count),hl(connect,count)
    integer,device,intent(in)::li(list_height,group_max,connect,count)
    real(8),device,intent(out)::result(9,count)
    real(8),shared::sums(9,128)
    real(8)::acc(9),momentum(3),half_momentum(3),center(3),total,mass,vc(3),vh(3)
    integer::thread,cell,entry,group,n,ids(8),h,k,stride,atom
    thread=threadIdx%x;cell=blockIdx%x;acc=0
    do atom=thread,ns(cell),128
      mass=m(atom,cell)
      acc(1:3)=acc(1:3)+mass*half(:,atom,cell)*half(:,atom,cell)
      acc(4:6)=acc(4:6)+mass*v(:,atom,cell)*v(:,atom,cell)
    end do
    do group=0,connect
      if(group==0) then
        n=3;k=nw(cell)
      else
        n=group+1;k=hl(group,cell)
      end if
      do entry=thread,k,128
        if(group==0) then
          ids(1:n)=wl(1:n,entry,cell)
        else
          ids(1:n)=li(1:n,entry,group,cell)
        end if
        total=0;momentum=0;half_momentum=0;center=0
        do h=1,n
          atom=ids(h);mass=m(atom,cell);total=total+mass
          momentum=momentum+mass*v(:,atom,cell)
          half_momentum=half_momentum+mass*half(:,atom,cell)
          center=center+mass*x(:,atom,cell)
        end do
        vc=momentum/total;vh=half_momentum/total;center=center/total
        acc(1:3)=acc(1:3)+total*vh*vh
        acc(4:6)=acc(4:6)+total*vc*vc
        do h=1,n
          atom=ids(h)
          acc(7:9)=acc(7:9)-(x(:,atom,cell)-center)*force(:,atom,cell)
        end do
      end do
    end do
    sums(:,thread)=acc
    call syncthreads()
    stride=64
    do while(stride>0)
      if(thread<=stride) sums(:,thread)=sums(:,thread)+sums(:,thread+stride)
      call syncthreads()
      stride=stride/2
    end do
    if(thread==1) result(:,cell)=sums(:,1)
  end subroutine
end module
