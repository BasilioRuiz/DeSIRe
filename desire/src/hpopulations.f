c calulate H populations in LTE

        subroutine hpopulations(T,pe,popH)

        implicit real*4 (a-h,o-z)
        real*4 popH(*)

        integer i,ncontr
        parameter (ncontr=28)
        real*4 u0(ncontr),u1(ncontr),u2(ncontr)
        real*4 du0,du1,du2
        real*4 pp(10),p(28),pe,pg,nH,kT,kboltzT,chi,chi10,u00

        kT=8.6173325e-5*T      !eV
        kboltzT=1.38054e-16*T  !ergs
        call gasc(T,pe,pg,pp)
        !pp(1)                     !p(h)/p(h')
        !pp(6)                     !p(h+)/p(h')
        !pp(7)                     !p(h2)/p(h')
        !pp(8)                     !pe/p(h')
        !pp(9)=pe/(1.38054e-16*t)  !n(e)=pe/kt

        nH=pp(1)/pp(8)*pe/kboltzT !n(HI)

        call nelfctb(1,T,u0(1),u1(1),u2(1),du0,du1,du2)

        chi10=13.59843449  !eV
        u00=u0(1)
        popH(1)=nH*2./u00
        do i=2,5
           chi=chi10*(1.-1./real(i*i))
           popH(i)=nH*2.*i*i*exp(-chi/kT)/u00
        end do
        popH(6)=pp(6)/pp(8)*pe/kboltzT

        return
        end
