        !COMPILER-GENERATED INTERFACE MODULE: Fri Oct 09 19:36:16 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INTERP__genmod
          INTERFACE 
            SUBROUTINE INTERP(BUF,T,NCF,NCM,NA,IFL,PV)
              INTEGER(KIND=4) :: NCM
              INTEGER(KIND=4) :: NCF
              REAL(KIND=8) :: BUF(NCF,NCM,*)
              REAL(KIND=8) :: T(2)
              INTEGER(KIND=4) :: NA
              INTEGER(KIND=4) :: IFL
              REAL(KIND=8) :: PV(NCM,*)
            END SUBROUTINE INTERP
          END INTERFACE 
        END MODULE INTERP__genmod
