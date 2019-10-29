c  writeabun_RH writes the abundances in the RH format


	subroutine write_abun_RH(abu,filename)  
	implicit real*4 (a-h,o-z)
	parameter (na=92)
	character*(*) filename
        real abu(*)

        CHARACTER*2 ATM(na)/'H','HE','LI','BE','B','C','N','O','F','NE',
     *'NA','MG','AL','SI','P','S','CL','AR','K','CA','SC','TI','V','CR',
     *'MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR',
     *'RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN',
     *'SN','SB','TE','I','XE','CS','BA','LA','CE','PR','ND','PM',
     *'SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA','W',
     *'RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN',
     *'FR','RA','AC','TH','PA','U'/
       CHARACTER*2 ATM2(7)/'NP','PU','AM','CM','BK','CF','ES'/
        
        ican=54
        open(ican,file=filename)
        do i=1,92 
	  write(ican,*)ATM(i),abu(i)
        end do
        do i=1,7
          write(ican,*)ATM2(i),-7.96
        end do

        close(ican)
        return
        end


