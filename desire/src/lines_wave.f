c lines_wave
c transfrom LTE data form grid file (SIR format) to a RH format file (.grid --> .wave)
c ntl: Number of lines (each one can have several blends, not considered in the total number ntl)
c nli: Number of wavelengths in the wavelength grid
c nble(iln) number of blends of each line
c
c 10/10/20 brc: First version.
c          epm: Save the wavetable in memory instead of file.
c_____________________________________________________________________________

        subroutine lines_wave( )

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        character*2 atom_nlte(kl)
        integer nlin(kl),npas(kl),nble(kl)
        integer ntl,ntau,ist(4)
        integer iln,ixx,ilte,ikk0,ntot,jj
        integer nlow_i(kl), nup_i(kl),linea_nlte(kl)
        real*4  dlongd(kld),dlamda0(kl)
        real*8  wavedbl_arr(kl),paso,waveini,w_air(kld),w_vac(kld)

c       Salva los valores cuando esta subrutina termine.
c       (Es innecesario si se compila con -fno-automatic)
        save w_vac

        common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/niveles/nlow_i,nup_i,atom_nlte,linea_nlte
        common/wavearrdble/wavedbl_arr
        common/ldeo/dlongd,dlamda0

        ixx=0
        ikk0=0  !contador sobre pasos
        ntot=0
        do iln=1,ntl
           ilte=0
           ix0=ixx
           do ible=1,nble(iln)
              ixx=ixx+1
              if(linea_nlte(ixx) .eq. 0)then
                 ilte=ilte+1
              end if
           end do

           if(ilte .ge. 1)then
              paso=(dlongd(ikk0+2)-dlongd(ikk0+1))*1.d-3  !en A
              waveini=wavedbl_arr(ix0+1)+dlongd(ikk0+1)*1.d-3
              do ipaso=-1,npas(iln)+1
                 ntot=ntot+1
                 w_air(ntot)=waveini+ipaso*paso
              end do
           endif
           ikk0=ikk0+npas(iln)
        end do
        call air_to_vacuum(ntot,w_air,w_vac)  !w_air en A, w_vac en nm

c       Pasa el array a RH donde se registra en una estructura C global.
        call sirwave(ntot,w_vac)

        return
        end
