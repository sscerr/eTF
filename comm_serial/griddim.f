          SUBROUTINE get_grid_dims(nproc,nprow,npcol)
! ...     This subroutine factorizes the number of processors (NPROC)
! ...     into NPROW and NPCOL,  that are the sizes of the 2D processors mesh.
            integer nproc,nprow,npcol
            integer sqrtnp,i
            sqrtnp = INT( SQRT( DBLE(nproc) ) + 1 )
            DO i=1,sqrtnp
              IF(MOD(nproc,i).EQ.0) nprow = i
            END DO
            npcol = nproc/nprow
            RETURN
          END SUBROUTINE

