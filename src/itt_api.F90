module itt_api
  use,intrinsic :: iso_c_binding
  implicit none

  interface
     !__itt_domain* __itt_domain_create( const __itt_char *name ); 
     function itt_domain_create(name) bind(c)
       import
       type(C_PTR),value :: name
       type(C_PTR) :: itt_domain_create
     end function itt_domain_create

     !void __itt_frame_begin_v3(const __itt_domain *domain, __itt_id *id); 
     subroutine itt_frame_begin(domain,id) bind(c)
       import
       type(C_PTR),value :: domain
       type(C_PTR),value :: id
     end subroutine itt_frame_begin

     !void __itt_frame_end_v3(const __itt_domain *domain, __itt_id *id);
     subroutine itt_frame_end(domain,id) bind(c)
       import
       type(C_PTR),value :: domain
       type(C_PTR),value :: id
     end subroutine itt_frame_end
  end interface
end module itt_api
