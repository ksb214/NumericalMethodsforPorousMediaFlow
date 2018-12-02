c---------------------------------------------------------------- zfact
      function zfact(p,t,grav,fco2,fh2s)
c----------------------------------------------------------------------
c  Calculates gas compressibility (z-factor).
c  INPUT Variables:
c   p, psia - pressure
c   t, Rankine - temperature
c   grav - gas specific gravity (to air)
c   fco2, fraction - mole fraction of the CO2 in the gas
c   fh2s, fraction - mole fraction of the H2S in the gas
c
c  OUTPUT
c   zfact - gas compressibility factor
c
      implicit real*8(a-h,o-z)
c........................................... check the pressure
      if(p .lt. 14.0) then
	write(*,85) p
 85   format(5x,'Pressure is too low: Set Z-factor = 1.0, P =',f10.2)
        z = 1.0
	goto 40
      end if	
c
c.............. pseudo-critical temp. & press. by Sutton (SPE,1985)
      pc = 756.8-131.0*grav-3.6*grav*grav
      tc = 169.2+349.5*grav-74.0*grav*grav
c..........................effect of non-hydrocarbon ipurities on Z
      a = fco2 + fh2s
      b = fh2s
      err = 120.0*(a**0.9 - a**1.6) + 15.*(b**0.5 - b**4)      
      pc = pc*(tc-err)/(tc+b*(1.-b)*err)
      tc = tc - err
c.................................... pseudo-reduced temp. & press.
      tr = t/tc
      pr = p/pc
c............................. Dranchuk & Abou-Kassem Eq.(JCP,1975)
      a1 = 0.3265
      a2 =-1.07
      a3 =-0.5339
      a4 = 0.01569
      a5 =-0.05165
      a6 = 0.5475
      a7 =-0.7361
      a8 = 0.1844
      a9 = 0.1056
      a10= 0.6134
      a11= 0.7210
c
      c1 = a1+a2/tr+a3/tr**3+a4/tr**4+a5/tr**5
      c2 = a6+a7/tr+a8/tr**2
      c3 = a9*(a7/tr+a8/tr**2)
c
      crit = 1.0E-6
      itcon=0
 20   z = 1.0
      do 30 iter = 1,20
         dr = 0.27*pr/(z*tr)
         c4 = a10*(1+a11*dr**2)*(dr**2/tr**3)*exp(-a11*dr**2)
         dc4dr = (2.d0*a10*dr/(tr**3))
     &           * (1.d0 + a11*(dr**2) - (a11*(dr**2))**2)
     &           * dexp(-a11*(dr**2))
         dzdr = c1 + 2.d0*c2*dr
     &           - 5.d0*c3*(dr**4) + dc4dr
c.................................... function statement for DAK Eq.
         fun = z - (1.d0 + c1*dr + c2*(dr**2)
     &             - c3*(dr**5) + c4)
         dfun = 1.d0+(c1+2.d0*c2*dr-5.d0*c3*(dr**4)+dc4dr)*dr/z
         deL = -(fun/dfun)
         z = z + deL
         if(dabs(del) .lt. crit) goto 40
 30   continue
      itcon = itcon + 1
      crit = crit*10.
      if(itcon .le. 2) goto 20      
      write(6,*) 'Warning: zfact did not converge pressure at.',p
      write(*,*) 'Warning: zfact did not converge pressure at.',p
      z = 1.0
 40   continue
      zfact=z
c
      return
      end
c
            