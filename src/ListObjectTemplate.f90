!#define LISTNAME(a) CONCAT(a,List)

!> Defines a List type for handling singly linked lists. This type is already specialized for the usage for 
!! GENE. It is used there for intermediate storage of entries of the gyromatrix (in gyro_average_df), 
!! before allocating the sparse banded matrix.
MODULE MODULELISTNAME
  USE MODULENODENAME
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: push, pop, LISTNAME, initialize, finalize, isEmpty, &
       &isAtEnd, get_next, init_iterator, length, output

  !> This type contains the data for one node of the list.
  TYPE Node
     type(TYPENAME) :: DATA
     TYPE(Node), pointer :: next !< pointer to the next Node in the List
  END TYPE Node

  !> To be able to have several lists, each list contains its own Head and Tail.
  TYPE LISTNAME
     private
     TYPE(Node), POINTER :: Head !< First Node in the List
     TYPE(Node), POINTER :: Tail !< Last Node in the List
     TYPE(Node), pointer :: iterator_node
     INTEGER :: NumberOfEntries
  END TYPE LISTNAME

  !> Initializes a List. This interface has to be called before the usage of the List.
  INTERFACE initialize
     module procedure lo_initialize
  END INTERFACE

  !> Finalizes a List, including deallocation of occupied memory. It has to be called
  !! after the last usage of the List.
  INTERFACE finalize
     MODULE PROCEDURE lo_finalize
  END INTERFACE

  !> Add a node to the end of the List.
  INTERFACE push
     MODULE PROCEDURE push_values
  END INTERFACE

  !> Get the first Node from the List and remove it from the List.
  INTERFACE pop
     module procedure pop_node
  END INTERFACE

  INTERFACE isAtEnd
     module procedure lo_isAtEnd
  END INTERFACE

  INTERFACE isEmpty
     module procedure lo_isEmpty
  END INTERFACE

  INTERFACE Length
     module procedure lo_length
  END INTERFACE

  INTERFACE get_next
     MODULE procedure lo_get_next
  END INTERFACE

  INTERFACE init_iterator
     module procedure lo_init_iterator
  END INTERFACE

  INTERFACE output
     module procedure lo_output_list
  END INTERFACE

CONTAINS
  SUBROUTINE lo_initialize(this)
    TYPE(LISTNAME) :: this

    ALLOCATE(this%Head)
    this%Head%next => null()
    this%Tail => this%Head
    this%NumberOfEntries = 0
  END SUBROUTINE lo_initialize
  
  SUBROUTINE lo_finalize(this)
    TYPE(LISTNAME) :: this
    
    ! Local variables
    type(TYPENAME) :: temp_node

    DO WHILE (ASSOCIATED(this%Head%next))
       CALL pop(this, temp_node)
    END DO
    IF (ASSOCIATED(this%Head)) DEALLOCATE(this%Head)
  END SUBROUTINE lo_finalize

  !> A Node type is added to the List at the end.
  SUBROUTINE push_node(this, ai_node)
    TYPE(LISTNAME) :: this
    TYPE(Node) :: ai_node !< contains the Node which is to be appended

    ALLOCATE(this%Tail%next)
    this%Tail => this%Tail%next
    this%Tail = ai_node
    this%Tail%next  => NULL()
    this%NumberOfEntries = this%NumberOfEntries + 1
  END SUBROUTINE push_node

  !> One can also append the data of a Node directly with this routine.
  SUBROUTINE push_values(this, nodedata)
    TYPE(LISTNAME) :: this
    TYPE(TYPENAME) :: nodedata

    ALLOCATE(this%Tail%next)
    this%Tail => this%Tail%next
    this%Tail%data = nodedata
    this%Tail%next  => NULL()
    this%NumberOfEntries = this%NumberOfEntries + 1
  END SUBROUTINE push_values

  !> returns the first node from the List 
  SUBROUTINE pop_node(this, ao_node)
    type(LISTNAME) :: this
    TYPE(TYPENAME), INTENT(OUT) :: ao_node !< return value, contains the first node of the List

    ! Local variables
    TYPE(Node), pointer :: node_ptr

    IF (ASSOCIATED(this%Head%next)) THEN
       ao_node = this%Head%next%data
       node_ptr => this%Head%next%next
       IF (ASSOCIATED(this%Tail,this%Head%next)) THEN
          ! last node is to be deallocated, so we have to reassign Tail
          this%Tail => this%Head
       END IF
       deallocate(this%Head%next)
       this%Head%next => node_ptr
       this%NumberOfEntries = this%NumberOfEntries - 1
    END IF
  END SUBROUTINE pop_node
    
  !> Test if the List has any content. Returns true if the List is empty, otherwise false.
  FUNCTION lo_isEmpty(this) 
    TYPE(LISTNAME) :: this
    logical :: lo_isEmpty

    !lo_isEmpty = (this%NumberOfEntries.EQ.0)
    lo_isEmpty = ASSOCIATED(this%Tail,this%Head)
  END FUNCTION lo_isEmpty

  FUNCTION lo_length(this)
    TYPE(LISTNAME) :: this
    integer :: lo_length

    lo_length = this%NumberOfEntries
  END FUNCTION lo_length

  !> returns a pointer to the Head Node
  FUNCTION get_ptr_to_head(this) RESULT(ptr_to_head)
    TYPE(LISTNAME) :: this
    TYPE(Node), pointer :: ptr_to_head

    ptr_to_head => this%Head
  END FUNCTION get_ptr_to_head

  FUNCTION lo_get_next(this) RESULT(nodedata)
    type(LISTNAME) :: this
    TYPE(TYPENAME),POINTER :: nodedata

    IF (ASSOCIATED(this%iterator_node%next)) THEN
       this%iterator_node=>this%iterator_node%next
       nodedata => this%iterator_node%data
    ELSE
       nodedata => null()
    END IF
  END FUNCTION lo_get_next

  FUNCTION lo_isAtEnd(this) 
    TYPE(LISTNAME) :: this
    logical :: lo_isAtEnd

    lo_isAtEnd = ASSOCIATED(this%iterator_node,this%Tail)
  END FUNCTION lo_isAtEnd

  SUBROUTINE lo_init_Iterator(this)
    TYPE(LISTNAME) :: this

    this%iterator_node => this%Head
  END SUBROUTINE lo_init_Iterator

  SUBROUTINE lo_output_list(this)
    TYPE(LISTNAME) :: this

    TYPE(Node), pointer :: ptr

    ptr => this%Head%next
    DO WHILE (ASSOCIATED(ptr))
       CALL output(ptr%data)
       ptr => ptr%next
    END DO
  END SUBROUTINE lo_output_list
END MODULE MODULELISTNAME
