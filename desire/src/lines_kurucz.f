c lines_kurucz
c transfrom an atomic parameters file (SIR FORMAT) to a Kurucz format file
c
c 07/07/20 brc: First version.
c          epm: Save Kurucz data in memory instead of file.
c_____________________________________________________________________________

        subroutine lines_kurucz( )

        implicit real*4 (a-h,o-z)
        include 'PARAMETER'

        integer*4 ist(4),nlin(kl),npas(kl),nble(kl)
        integer*4 ntau,ntl,nlines

        real*8    wavenm,wlengt_air,wlengt_vac,wavedbl_arr(kl)
        real*8    aloggamma_rad,aloggamma_4,aloggamma_6
        real*8    energy_low,energy_up
        real*8    eVtocm_1
        real*8    w_vac_fn  !funcion en departures.f
        parameter (eVtocm_1=8065.5443)

c       Otras variables de los commons.
        real*8    elemcode_arr(kl)
        real*8    loggf8_arr(kl)
        real*4    gf_arr(kl),energy_arr(kl)
        integer*4 nlow_i(kl), nup_i(kl),linea_nlte(kl)
        character atom_nlte(kl)*2

c       Variables a pasar a RH.
        real*8    field1(kl),field2(kl),field3(kl),field4(kl)
        real*8    field5(kl),field6(kl),field7(kl),field8(kl)

c       Salva los valores cuando esta subrutina termine.
c       (Es innecesario si se compila con -fno-automatic)
        save field1,field2,field3,field4,field5,field6,field7,field8

        common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/wavearrdble/wavedbl_arr
        common/elemcode/elemcode_arr
        common/loggfarr/gf_arr,energy_arr
        common/loggfarr8/loggf8_arr
        common/niveles/nlow_i,nup_i,linea_nlte,atom_nlte

        ixx=0
        nlines=0
        do iln=1,ntl  !numero total de lineas
           do ible=1,nble(iln)  !numero de blends de cada linea
              ixx=ixx+1

              if(nlow_i(ixx).eq.0 .and. nup_i(ixx).eq.0)then  !LTE case
                 nlines=nlines+1

                 wlengt_air=wavedbl_arr(ixx)
                 wavenm=wlengt_air/10.0d0
                 wlengt_vac=w_vac_fn(wlengt_air)*1.0d-8  !vacumm ldo in cm

                 energy_low=dble(energy_arr(ixx))*eVtocm_1
                 energy_up=energy_low+1.0d0/wlengt_vac

                 aloggamma_rad=dlog(0.22233d0/wlengt_vac)
                 aloggamma_4=-6.00d0
                 aloggamma_6=-7.60d0

                 field1(nlines) = wavenm
                 field2(nlines) = loggf8_arr(ixx)
                 field3(nlines) = elemcode_arr(ixx)
                 field4(nlines) = energy_low
                 field5(nlines) = energy_up
                 field6(nlines) = aloggamma_rad
                 field7(nlines) = aloggamma_4
                 field8(nlines) = aloggamma_6

              end if
           end do
        end do

        call sirkurucz(nlines,field1,field2,field3,field4,
     &                        field5,field6,field7,field8)

        return
        end
