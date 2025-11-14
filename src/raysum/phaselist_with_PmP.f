c####&

c Subroutines to list phases of different types.


c Direct phases, excluding PSP conversions (low amplitude):
c (Modified to make unconverted phase always appear first)
      subroutine ph_direct(phaselist,nseg,numph,nlay,iphase)
      
        implicit none
        include 'params.h'
        
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph
        integer nlay,iseg,iph,itemp,iphase,tempphase,conv_count
        integer phase_order(3)
        
        do iph=1,3
          phase_order(iph)=iphase+iph-1
          if (phase_order(iph) .gt. 3) then
            phase_order(iph)=phase_order(iph)-3
          end if
          write(*,*) 'Phase order and iph: ', phase_order(iph), ' ', iph
        end do

        do iph=1,(3**(nlay-1))
          itemp=numph+1
c          write(*,*) 'Current iph: itemp=',itemp
          nseg(itemp)=nlay
          phaselist(1,1,itemp)=nlay
          phaselist(1,2,itemp)=iphase
c          write(*,*) 'Start phaselist:',phaselist(1,1,itemp),
c     &               phaselist(1,2,itemp)
          tempphase=iph-1
          do iseg=nlay,2,-1          
            phaselist(iseg,1,itemp)=nlay-iseg+1
            phaselist(iseg,2,itemp)=tempphase-(tempphase/3)*3+1
            phaselist(iseg,2,itemp)=phase_order(phaselist(iseg,2,itemp))
            tempphase=tempphase/3
c            write(*,*) 'Currrent phaselist:',phaselist(iseg,1,itemp),
c     &                phaselist(iseg,2,itemp)
          end do
c Check if current phase involves multiple conversions
          if (conv_count(phaselist(1,1,itemp),nseg(itemp))
     &               .le. 1) then
            numph=numph+1
          end if
        end do        

      end
      
c --------------------------------------
      
c Direct phases, including all possible conversions
      subroutine ph_direct_all(phaselist,nseg,numph,nlay,iphase)       
        
        implicit none
        include 'params.h'
        
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph
        integer nlay,iseg,iph,itemp,iphase,tempphase
        
        do iph=1,(3**(nlay-1))
          itemp=numph+1
          nseg(itemp)=nlay
          phaselist(1,1,itemp)=nlay
          phaselist(1,2,itemp)=iphase
          tempphase=iph-1
          do iseg=nlay,2,-1          
            phaselist(iseg,1,itemp)=nlay-iseg+1
            phaselist(iseg,2,itemp)=tempphase-(tempphase/3)*3+1
            tempphase=tempphase/3
          end do
          numph=numph+1
        end do
                
      end

c ----------------------------------

c Get surface multiples for a given interface. Conversions are
c assumed to happen only at the surface and the bottoming layer.
c blay is the bottoming layer for the multiples; iphase is the
c incident phase.
      subroutine ph_fsmults(phaselist,nseg,numph,nlay,blay,iphase)
      
        implicit none
        include 'params.h'
        
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph
        integer nlay,blay,iphase

        if (iphase .eq. 1) then
          call ph_fsmults_P(phaselist,nseg,numph,nlay,blay)
        else
          call ph_fsmults_S(phaselist,nseg,numph,nlay,blay,iphase)
        end if
        
      end  
      
      
c -------------------------------------------

c Incident-P case for ph_fsmults
      subroutine ph_fsmults_P(phaselist,nseg,numph,nlay,blay)
      
        implicit none
        include 'params.h'
        
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph
        integer nlay,blay,iph,tempphase,iseg
        
c        Unconverted multiple
        numph=numph+1
c        write(*,*) 'Current phase number unconverted:', numph
        nseg(numph)=nlay+2*blay
        do iseg=1,nlay
          phaselist(iseg,1,numph)=nlay-iseg+1
          phaselist(iseg,2,numph)=1
        end do
        do iseg=(nlay+1),(nlay+blay)
          phaselist(iseg,1,numph)=iseg-nlay
          phaselist(iseg,2,numph)=4
        end do
        do iseg=(nlay+blay+1),nseg(numph)
          phaselist(iseg,1,numph)=nseg(numph)-iseg+1
          phaselist(iseg,2,numph)=1
        end do

c        Conversion at blay bounce. Consider all possible qS
c        polarizations
        do iph=1,2**blay
          numph=numph+1
c          write(*,*) 'Current phase number blay bounce:', numph
          nseg(numph)=nlay+2*blay
          do iseg=1,nlay
            phaselist(iseg,1,numph)=nlay-iseg+1
            phaselist(iseg,2,numph)=1          
          end do
          do iseg=(nlay+1),(nlay+blay)
            phaselist(iseg,1,numph)=iseg-nlay
            phaselist(iseg,2,numph)=4
          end do 
          tempphase=iph-1
          do iseg=nseg(numph),(nlay+blay+1),-1
            phaselist(iseg,1,numph)=nseg(numph)-iseg+1
            phaselist(iseg,2,numph)=tempphase-(tempphase/2)*2+2
            tempphase=tempphase/2
          end do
        end do

c        Conversion at surface
        do iph=1,2**(blay*2)
          numph=numph+1
          nseg(numph)=nlay+2*blay
          do iseg=1,nlay
            phaselist(iseg,1,numph)=nlay-iseg+1
            phaselist(iseg,2,numph)=1          
          end do
          tempphase=iph-1
          do iseg=nseg(numph),(nlay+1),-1
            phaselist(iseg,2,numph)=tempphase-(tempphase/2)*2+2
            tempphase=tempphase/2
            if (iseg .le. (nlay+blay)) then
              phaselist(iseg,2,numph)=phaselist(iseg,2,numph)+3
              phaselist(iseg,1,numph)=iseg-nlay
            else
              phaselist(iseg,1,numph)=nseg(numph)-iseg+1
            end if
          end do     
        end do

c        Conversion at blay transmission -- left out for now.      
        
      end


       
c Incident-S case for ph_fsmults
      subroutine ph_fsmults_S(phaselist,nseg,numph,nlay,blay,iphase)
      
        implicit none
        include 'params.h'
        
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph
        integer nlay,blay,iph,tempphase,iseg,iphase
        
c        Unconverted multiples (S variants)
        do iph=1,2**(nlay-1+2*blay)
          numph=numph+1
          nseg(numph)=nlay+2*blay
          phaselist(1,1,numph)=nlay
          phaselist(1,2,numph)=iphase
          tempphase=iph-1
          do iseg=nseg(numph),(nseg(numph)-blay+1),-1
            phaselist(iseg,1,numph)=nseg(numph)-iseg+1
            phaselist(iseg,2,numph)=tempphase-(tempphase/2)*2+2
            tempphase=tempphase/2
          end do
          do iseg=(nseg(numph)-blay),(nseg(numph)-2*blay+1),-1
            phaselist(iseg,1,numph)=iseg-nseg(numph)+2*blay
            phaselist(iseg,2,numph)=tempphase-(tempphase/2)*2+5
            tempphase=tempphase/2
          end do
          do iseg=(nseg(numph)-2*blay),2,-1
            phaselist(iseg,1,numph)=nlay-iseg+1
            phaselist(iseg,2,numph)=tempphase-(tempphase/2)*2+2
            tempphase=tempphase/2
          end do
        end do

c        Conversion at blay bounce.
        do iph=1,2**(nlay-1+blay)
          numph=numph+1
          nseg(numph)=nlay+2*blay
          phaselist(1,1,numph)=nlay
          phaselist(1,2,numph)=iphase
          do iseg=nseg(numph),(nseg(numph)-blay+1),-1
            phaselist(iseg,1,numph)=nseg(numph)-iseg+1
            phaselist(iseg,2,numph)=1
          end do
          tempphase=iph-1
          do iseg=(nseg(numph)-blay),(nseg(numph)-2*blay+1),-1
            phaselist(iseg,1,numph)=iseg-nseg(numph)+2*blay
            phaselist(iseg,2,numph)=tempphase-(tempphase/2)*2+5
            tempphase=tempphase/2
          end do
          do iseg=(nseg(numph)-2*blay),2,-1
            phaselist(iseg,1,numph)=nlay-iseg+1
            phaselist(iseg,2,numph)=tempphase-(tempphase/2)*2+2
            tempphase=tempphase/2
          end do
        end do

c        Conversion at surface reflection
        do iph=1,2**(nlay-1)
          numph=numph+1
          nseg(numph)=nlay+2*blay
          phaselist(1,1,numph)=nlay
          phaselist(1,2,numph)=iphase
          do iseg=nseg(numph),(nseg(numph)-blay+1),-1
            phaselist(iseg,1,numph)=nseg(numph)-iseg+1
            phaselist(iseg,2,numph)=1
          end do
          do iseg=(nseg(numph)-blay),(nseg(numph)-2*blay+1),-1
            phaselist(iseg,1,numph)=iseg-nseg(numph)+2*blay
            phaselist(iseg,2,numph)=1
          end do
          tempphase=iph-1
          do iseg=(nseg(numph)-2*blay),2,-1
            phaselist(iseg,1,numph)=nlay-iseg+1
            phaselist(iseg,2,numph)=tempphase-(tempphase/2)*2+2
            tempphase=tempphase/2
          end do
        end do


c        Conversion at blay transmission -- left out for now.      
        
      end      
c -----------------------------------
c Sobroutine to create phases important for receiver function calculation 
c These phases include (in additon to phases crated from ph_direct subroutine)
c the folowing phases: PpPp, PpPs, PpSs + PsPs and PsSs.
      subroutine ph_rfmults(phaselist,nseg,numph,nlay,blay,iphase)
      
        implicit none
        include 'params.h'
        
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph
        integer nlay,blay,iphase

        if (iphase .eq. 1) then
          call ph_rfmults_P(phaselist,nseg,numph,nlay,blay)
c        else
c          call ph_rfmults_S(phaselist,nseg,numph,nlay,blay,iphase)
        end if
        
      end 

c -------------------------------------------

c Incident-P case for ph_rfmults
      subroutine ph_rfmults_P(phaselist,nseg,numph,nlay,blay)
      
        implicit none
        include 'params.h'
        
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph
        integer nlay,blay,iph,tempphase,iseg
        
        
c         This part has been commented out as we do not need 
c         unconverted P waves in RFs, since in the process of 
c         deconvolution they are removed anyway.
c        Unconverted multiple
        numph=numph+1
c        write(*,*) 'Current phase number unconverted:', numph
        nseg(numph)=nlay+2*blay
        do iseg=1,nlay
          phaselist(iseg,1,numph)=nlay-iseg+1
          phaselist(iseg,2,numph)=1
        end do
        do iseg=(nlay+1),(nlay+blay)
          phaselist(iseg,1,numph)=iseg-nlay
          phaselist(iseg,2,numph)=4
        end do
        do iseg=(nlay+blay+1),nseg(numph)
          phaselist(iseg,1,numph)=nseg(numph)-iseg+1
          phaselist(iseg,2,numph)=1
        end do

c        Conversion at blay bounce. Consider all possible qS
c        polarizations and create phases
        do iph=1,2
          numph=numph+1
c          write(*,*) 'Current phase number blay bounce:', numph
          nseg(numph)=nlay+2*blay
          do iseg=1,nlay
            phaselist(iseg,1,numph)=nlay-iseg+1
            phaselist(iseg,2,numph)=1          
          end do
          do iseg=(nlay+1),(nlay+blay)
            phaselist(iseg,1,numph)=iseg-nlay
            phaselist(iseg,2,numph)=4
          end do 
          do iseg=nseg(numph),(nlay+blay+1),-1
            phaselist(iseg,1,numph)=nseg(numph)-iseg+1
            phaselist(iseg,2,numph)=iph+1
          end do
        end do

c        Conversion at surface
        do iph=1,2**2
          numph=numph+1
          nseg(numph)=nlay+2*blay
          do iseg=1,nlay
            phaselist(iseg,1,numph)=nlay-iseg+1
            phaselist(iseg,2,numph)=1          
          end do
          do iseg=(nlay+1),(nlay+blay)
            phaselist(iseg,1,numph)=iseg-nlay
c           This will adjust phases going down to 5 and 6
c           that is sign for Sv and Sh          
            phaselist(iseg,2,numph)=(iph-1)-((iph-1)/2)*2+5
          end do
          do iseg=(nlay+blay+1),nseg(numph)
            phaselist(iseg,1,numph)=nseg(numph)-iseg+1
            phaselist(iseg,2,numph)=((iph-1)/2)+2
          end do  
        end do

c        Conversion at blay transmission, namley PsSs and PsPs.
        do iph=1,2**2
          numph=numph+1
          nseg(numph)=nlay+2*blay
c         P phase in the halh-space
          phaselist(1,1,numph)=nlay
          phaselist(1,2,numph)=1
          do iseg=2,nlay
            phaselist(iseg,1,numph)=nlay-iseg+1
            phaselist(iseg,2,numph)=(iph-1)-((iph-1)/2)*2+2          
          end do
c         Part going down
          do iseg=(nlay+1),(nlay+blay)
            phaselist(iseg,1,numph)=iseg-nlay
            if (iph .le. 2) then
              phaselist(iseg,2,numph)=(iph-1)-((iph-1)/2)*2+5
            else
              phaselist(iseg,2,numph)=4
            end if
          end do
c         Last leg 
          do iseg=nseg(numph),(nlay+blay+1),-1
            phaselist(iseg,1,numph)=nseg(numph)-iseg+1
            phaselist(iseg,2,numph)=(iph-1)-((iph-1)/2)*2+2 
          end do
        end do

      end

c RECEIVER FUNCTION PHASE LISTING
c Subroutines to list phases for receiver function modeling.

c      Direct phases, conversions only at Moho.
      subroutine ph_direct_rf(phaselist,nseg,numph,nlay,blay,iphase)
        
        implicit none
        include 'params.h'
        
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph
        integer nlay,blay,iseg,iph,iphase

        if (blay .eq. 1) then        
          do iph=1,3
            numph=numph+1
            nseg(numph)=nlay
            do iseg=1,nlay-blay    
              phaselist(iseg,1,numph)=nlay-iseg+1
              phaselist(iseg,2,numph)=iphase
            end do
            do iseg=nlay-blay+1,nlay   
              phaselist(iseg,1,numph)=nlay-iseg+1
              phaselist(iseg,2,numph)=(iph-1)-((iph-1)/3)*3+1
            end do       
          end do
        else 
          do iph=1,2
            numph=numph+1
            nseg(numph)=nlay
            do iseg=1,nlay-blay    
              phaselist(iseg,1,numph)=nlay-iseg+1
              phaselist(iseg,2,numph)=iphase
            end do
            do iseg=nlay-blay+1,nlay   
              phaselist(iseg,1,numph)=nlay-iseg+1
              phaselist(iseg,2,numph)=(iph-1)-((iph-1)/2)*2+2
            end do       
          end do
        end if 
                
      end

   
c -----------------------------------


c Count conversions in a given phase
      integer function conv_count(phase,nseg)
      
        implicit none
        include 'params.h'
        
        integer phase(maxseg,2),nseg,iseg,lastphase,thisphase,cc
        
        cc=0
        lastphase=phase(1,2)
        do iseg=2,nseg
          thisphase=phase(iseg,2)         
          if (((lastphase .eq. 1) .and. (thisphase .ne. 1)) .or.
     &        ((lastphase .ne. 1) .and. (thisphase .eq. 1))) then
            cc=cc+1
          end if
          lastphase=thisphase
        end do
        
        conv_count=cc
        
      end


c -----------------------------------

        
        
