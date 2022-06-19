    MODULE PYGDMPRECISION
! Define precision of variables
!if you want to run in single-precision set  DP=KIND(0.E0)
      INTEGER,PARAMETER :: DP=KIND(0.E0)

!if you want to run in double-precision set  DP=KIND(0.D0)
!       INTEGER,PARAMETER :: DP=KIND(0.D0)

    END MODULE PYGDMPRECISION