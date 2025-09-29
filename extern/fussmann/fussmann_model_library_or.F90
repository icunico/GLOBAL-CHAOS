module fussmann_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use fussmann_P1 
   use fussmann_C1
   use fussmann_X
   use fussmann_Y
   use fussmann_Z

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: fussmann_model_factory

contains

   subroutine create(self, name, model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('P1'); allocate(type_fussmann_P1::model)
         case ('C1'); allocate(type_fussmann_C1::model)
         case ('X'); allocate(type_fussmann_X::model)
         case ('Y'); allocate(type_fussmann_Y::model)
         case ('Z'); allocate(type_fussmann_Z::model)
      end select
   end subroutine create

end module fussmann_model_library
