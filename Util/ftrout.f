      subroutine dradb2 (ido,l1,cc,ch,wa1)
      double precision cc(ido,2,l1), ch(ido,l1,2), wa1(1), ti2, tr2
c
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
            tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
            ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
            ti2 = cc(i,1,k)+cc(ic,2,k)
            ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
            ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
  103    continue
  104 continue
c
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
         ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
  106 continue
c
  107 return
      end
      subroutine dradb3 (ido,l1,cc,ch,wa1,wa2)
      double precision cc(ido,3,l1), ch(ido,l1,3), wa1(1), wa2(1),
     1  ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui, taur, ti2, tr2
      data taur / -0.5 d0 /
      data taui  /  0.8660254037 8443864676 3723170752 93618d0/
c
      do 101 k=1,l1
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ci3 = taui*(cc(1,3,k)+cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
  101 continue
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
  102    continue
  103 continue
c
      return
      end
      subroutine dradb4 (ido,l1,cc,ch,wa1,wa2,wa3)
      double precision cc(ido,4,l1), ch(ido,l1,4), wa1(1), wa2(1),
     1  wa3(1), ci2, ci3, ci4, cr2, cr3, cr4, sqrt2, ti1, ti2, ti3, ti4,
     2  tr1, tr2, tr3, tr4
      data sqrt2 /  1.414213562 3730950488 0168872420 970 d0 /
c
      do 101 k=1,l1
         tr1 = cc(1,1,k)-cc(ido,4,k)
         tr2 = cc(1,1,k)+cc(ido,4,k)
         tr3 = cc(ido,2,k)+cc(ido,2,k)
         tr4 = cc(1,3,k)+cc(1,3,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,2) = tr1-tr4
         ch(1,k,3) = tr2-tr3
         ch(1,k,4) = tr1+tr4
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ti1 = cc(i,1,k)+cc(ic,4,k)
            ti2 = cc(i,1,k)-cc(ic,4,k)
            ti3 = cc(i,3,k)-cc(ic,2,k)
            tr4 = cc(i,3,k)+cc(ic,2,k)
            tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
            tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
            ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1-tr4
            cr4 = tr1+tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
            ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
            ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
            ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
            ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
            ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
c
  105 continue
      do 106 k=1,l1
         ti1 = cc(1,2,k)+cc(1,4,k)
         ti2 = cc(1,4,k)-cc(1,2,k)
         tr1 = cc(ido,1,k)-cc(ido,3,k)
         tr2 = cc(ido,1,k)+cc(ido,3,k)
         ch(ido,k,1) = tr2+tr2
         ch(ido,k,2) = sqrt2*(tr1-ti1)
         ch(ido,k,3) = ti2+ti2
         ch(ido,k,4) = -sqrt2*(tr1+ti1)
  106 continue
c
  107 return
      end
      subroutine dradb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      double precision cc(ido,5,l1), ch(ido,l1,5), wa1(1), wa2(1),
     1  wa3(1), wa4(1), ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5,
     2  di2, di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3,
     3  ti4, ti5, tr11, tr12, tr2, tr3, tr4, tr5
      data tr11  /  0.3090169943 7494742410 2293417182 81906d0/
      data ti11  /  0.9510565162 9515357211 6439333379 38214d0/
      data tr12  / -0.8090169943 7494742410 2293417182 81906d0/
      data ti12  /  0.5877852522 9247312916 8705954639 07277d0/
c
      do 101 k=1,l1
         ti5 = cc(1,3,k)+cc(1,3,k)
         ti4 = cc(1,5,k)+cc(1,5,k)
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         tr3 = cc(ido,4,k)+cc(ido,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci5 = ti11*ti5+ti12*ti4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(1,k,5) = cr2+ci5
  101 continue
      if (ido .eq. 1) return
c
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            ti5 = cc(i,3,k)+cc(ic,2,k)
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ti4 = cc(i,5,k)+cc(ic,4,k)
            ti3 = cc(i,5,k)-cc(ic,4,k)
            tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
            tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
  102    continue
  103 continue
c
      return
      end
      subroutine dradbg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      double precision cc(ido,ip,l1), c1(ido,l1,ip), c2(idl1,ip),
     1  ch(ido,l1,ip), ch2(idl1,ip), wa(1), ai1, ai2, ar1, ar1h, ar2,
     2  ar2h, arg, dc2, dcp, ds2, dsp, tpi
      data tpi   /  6.2831853071 7958647692 5286766559 00577d0/
c
      arg = tpi/dfloat(ip)
      dcp = dcos(arg)
      dsp = dsin(arg)
      idp2 = ido+2
      nbd = (ido-1)/2
      ipp2 = ip+2
      ipph = (ip+1)/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  101    continue
  102 continue
      go to 106
c
  103 do 105 i=1,ido
         do 104 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
c
  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 107 k=1,l1
            ch(1,k,j) = cc(ido,j2-2,k)+cc(ido,j2-2,k)
            ch(1,k,jc) = cc(1,j2-1,k)+cc(1,j2-1,k)
  107    continue
  108 continue
c
      if (ido .eq. 1) go to 116
      if (nbd .lt. l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k=1,l1
            do 109 i=3,ido,2
               ic = idp2-i
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  109       continue
  110    continue
  111 continue
      go to 116
c
  112 do 115 j=2,ipph
         jc = ipp2-j
         do 114 i=3,ido,2
            ic = idp2-i
            do 113 k=1,l1
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  113       continue
  114    continue
  115 continue
c
  116 ar1 = 1.
      ai1 = 0.
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 117 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+ar1*ch2(ik,2)
            c2(ik,lc) = ai1*ch2(ik,ip)
  117    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 118 ik=1,idl1
               c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc)
  118       continue
  119    continue
  120 continue
c
      do 122 j=2,ipph
         do 121 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  121    continue
  122 continue
c
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 k=1,l1
            ch(1,k,j) = c1(1,k,j)-c1(1,k,jc)
            ch(1,k,jc) = c1(1,k,j)+c1(1,k,jc)
  123    continue
  124 continue
c
      if (ido .eq. 1) go to 132
      if (nbd .lt. l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k=1,l1
            do 125 i=3,ido,2
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  125       continue
  126    continue
  127 continue
      go to 132
c
  128 do 131 j=2,ipph
         jc = ipp2-j
         do 130 i=3,ido,2
            do 129 k=1,l1
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  129       continue
  130    continue
  131 continue
  132 continue
c
      if (ido .eq. 1) return
      do 133 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  133 continue
c
      do 135 j=2,ip
         do 134 k=1,l1
            c1(1,k,j) = ch(1,k,j)
  134    continue
  135 continue
c
      if (nbd .gt. l1) go to 139
      is = -ido
      do 138 j=2,ip
         is = is+ido
         idij = is
         do 137 i=3,ido,2
            idij = idij+2
            do 136 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  136       continue
  137    continue
  138 continue
      go to 143
c
  139 is = -ido
      do 142 j=2,ip
         is = is+ido
         do 141 k=1,l1
            idij = is
            do 140 i=3,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  140       continue
  141    continue
  142 continue
c
  143 return
      end
      subroutine dradf2 (ido,l1,cc,ch,wa1)
      double precision cc(ido,l1,2), ch(ido,2,l1), wa1(1), ti2, tr2
c
      do 101 k=1,l1
         ch(1,1,k) = cc(1,k,1)+cc(1,k,2)
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,2)
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k) = cc(i,k,1)+ti2
            ch(ic,2,k) = ti2-cc(i,k,1)
            ch(i-1,1,k) = cc(i-1,k,1)+tr2
            ch(ic-1,2,k) = cc(i-1,k,1)-tr2
  103    continue
  104 continue
c
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(1,2,k) = -cc(ido,k,2)
         ch(ido,1,k) = cc(ido,k,1)
  106 continue
c
  107 return
      end
      subroutine dradf3 (ido,l1,cc,ch,wa1,wa2)
      double precision cc(ido,l1,3), ch(ido,3,l1), wa1(1), wa2(1),
     1  ci2, cr2, di2, di3, dr2, dr3, taui, taur, ti2, ti3, tr2, tr3
      data taur / -0.5 d0 /
      data taui  /  0.8660254037 8443864676 3723170752 93618d0/
c
      do 101 k=1,l1
         cr2 = cc(1,k,2)+cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2
         ch(1,3,k) = taui*(cc(1,k,3)-cc(1,k,2))
         ch(ido,2,k) = cc(1,k,1)+taur*cr2
  101 continue
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2+dr3
            ci2 = di2+di3
            ch(i-1,1,k) = cc(i-1,k,1)+cr2
            ch(i,1,k) = cc(i,k,1)+ci2
            tr2 = cc(i-1,k,1)+taur*cr2
            ti2 = cc(i,k,1)+taur*ci2
            tr3 = taui*(di2-di3)
            ti3 = taui*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(ic-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(ic,2,k) = ti3-ti2
  102    continue
  103 continue
c
      return
      end
      subroutine dradf4 (ido,l1,cc,ch,wa1,wa2,wa3)
      double precision cc(ido,l1,4), ch(ido,4,l1), wa1(1), wa2(1),
     1  wa3(1), ci2, ci3, ci4, cr2, cr3, cr4, hsqt2, ti1, ti2, ti3,
     2  ti4, tr1, tr2, tr3, tr4
      data hsqt2 /   .7071067811 8654752440 0844362104 85 d0 /
c
      do 101 k=1,l1
         tr1 = cc(1,k,2)+cc(1,k,4)
         tr2 = cc(1,k,1)+cc(1,k,3)
         ch(1,1,k) = tr1+tr2
         ch(ido,4,k) = tr2-tr1
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
         ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1 = cr2+cr4
            tr4 = cr4-cr2
            ti1 = ci2+ci4
            ti4 = ci2-ci4
            ti2 = cc(i,k,1)+ci3
            ti3 = cc(i,k,1)-ci3
            tr2 = cc(i-1,k,1)+cr3
            tr3 = cc(i-1,k,1)-cr3
            ch(i-1,1,k) = tr1+tr2
            ch(ic-1,4,k) = tr2-tr1
            ch(i,1,k) = ti1+ti2
            ch(ic,4,k) = ti1-ti2
            ch(i-1,3,k) = ti4+tr3
            ch(ic-1,2,k) = tr3-ti4
            ch(i,3,k) = tr4+ti3
            ch(ic,2,k) = tr4-ti3
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
c
      do 106 k=1,l1
         ti1 = -hsqt2*(cc(ido,k,2)+cc(ido,k,4))
         tr1 = hsqt2*(cc(ido,k,2)-cc(ido,k,4))
         ch(ido,1,k) = tr1+cc(ido,k,1)
         ch(ido,3,k) = cc(ido,k,1)-tr1
         ch(1,2,k) = ti1-cc(ido,k,3)
         ch(1,4,k) = ti1+cc(ido,k,3)
  106 continue
c
  107 return
      end
      subroutine dradf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      double precision cc(ido,l1,5), ch(ido,5,l1), wa1(1), wa2(1),
     1  wa3(1), wa4(1), ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2,
     2  di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3, ti4,
     3  ti5, tr11, tr12, tr2, tr3, tr4, tr5
      data tr11  /  0.3090169943 7494742410 2293417182 81906d0/
      data ti11  /  0.9510565162 9515357211 6439333379 38214d0/
      data tr12  / -0.8090169943 7494742410 2293417182 81906d0/
      data ti12  /  0.5877852522 9247312916 8705954639 07277d0/
c
      do 101 k=1,l1
         cr2 = cc(1,k,5)+cc(1,k,2)
         ci5 = cc(1,k,5)-cc(1,k,2)
         cr3 = cc(1,k,4)+cc(1,k,3)
         ci4 = cc(1,k,4)-cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2+cr3
         ch(ido,2,k) = cc(1,k,1)+tr11*cr2+tr12*cr3
         ch(1,3,k) = ti11*ci5+ti12*ci4
         ch(ido,4,k) = cc(1,k,1)+tr12*cr2+tr11*cr3
         ch(1,5,k) = ti12*ci5-ti11*ci4
  101 continue
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
            di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
            cr2 = dr2+dr5
            ci5 = dr5-dr2
            cr5 = di2-di5
            ci2 = di2+di5
            cr3 = dr3+dr4
            ci4 = dr4-dr3
            cr4 = di3-di4
            ci3 = di3+di4
            ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
            ch(i,1,k) = cc(i,k,1)+ci2+ci3
            tr2 = cc(i-1,k,1)+tr11*cr2+tr12*cr3
            ti2 = cc(i,k,1)+tr11*ci2+tr12*ci3
            tr3 = cc(i-1,k,1)+tr12*cr2+tr11*cr3
            ti3 = cc(i,k,1)+tr12*ci2+tr11*ci3
            tr5 = ti11*cr5+ti12*cr4
            ti5 = ti11*ci5+ti12*ci4
            tr4 = ti12*cr5-ti11*cr4
            ti4 = ti12*ci5-ti11*ci4
            ch(i-1,3,k) = tr2+tr5
            ch(ic-1,2,k) = tr2-tr5
            ch(i,3,k) = ti2+ti5
            ch(ic,2,k) = ti5-ti2
            ch(i-1,5,k) = tr3+tr4
            ch(ic-1,4,k) = tr3-tr4
            ch(i,5,k) = ti3+ti4
            ch(ic,4,k) = ti4-ti3
  102    continue
  103 continue
c
      return
      end
      subroutine dradfg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      double precision cc(ido,ip,l1), c1(ido,l1,ip), c2(idl1,ip),
     1  ch(ido,l1,ip), ch2(idl1,ip), wa(1), ai1, ai2, ar1, ar1h, ar2,
     2  ar2h, arg, dc2, dcp, ds2, dsp, tpi
      data tpi   /  6.2831853071 7958647692 5286766559 00577d0/
c
      arg = tpi/dfloat(ip)
      dcp = dcos(arg)
      dsp = dsin(arg)
      ipph = (ip+1)/2
      ipp2 = ip+2
      idp2 = ido+2
      nbd = (ido-1)/2
      if (ido .eq. 1) go to 119
      do 101 ik=1,idl1
         ch2(ik,1) = c2(ik,1)
  101 continue
      do 103 j=2,ip
         do 102 k=1,l1
            ch(1,k,j) = c1(1,k,j)
  102    continue
  103 continue
c
      if (nbd .gt. l1) go to 107
      is = -ido
      do 106 j=2,ip
         is = is+ido
         idij = is
         do 105 i=3,ido,2
            idij = idij+2
            do 104 k=1,l1
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  104       continue
  105    continue
  106 continue
      go to 111
c
  107 is = -ido
      do 110 j=2,ip
         is = is+ido
         do 109 k=1,l1
            idij = is
            do 108 i=3,ido,2
               idij = idij+2
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  108       continue
  109    continue
  110 continue
c
  111 if (nbd .lt. l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k=1,l1
            do 112 i=3,ido,2
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  112       continue
  113    continue
  114 continue
      go to 121
c
  115 do 118 j=2,ipph
         jc = ipp2-j
         do 117 i=3,ido,2
            do 116 k=1,l1
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  116       continue
  117    continue
  118 continue
      go to 121
c
  119 do 120 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  120 continue
c
  121 do 123 j=2,ipph
         jc = ipp2-j
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)+ch(1,k,jc)
            c1(1,k,jc) = ch(1,k,jc)-ch(1,k,j)
  122    continue
  123 continue
c
      ar1 = 1.d0
      ai1 = 0.d0
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 124 ik=1,idl1
            ch2(ik,l) = c2(ik,1)+ar1*c2(ik,2)
            ch2(ik,lc) = ai1*c2(ik,ip)
  124    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 125 ik=1,idl1
               ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j)
               ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc)
  125       continue
  126    continue
  127 continue
c
      do 129 j=2,ipph
         do 128 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+c2(ik,j)
  128    continue
  129 continue
c
      if (ido .lt. l1) go to 132
      do 131 k=1,l1
         do 130 i=1,ido
            cc(i,1,k) = ch(i,k,1)
  130    continue
  131 continue
      go to 135
c
  132 do 134 i=1,ido
         do 133 k=1,l1
            cc(i,1,k) = ch(i,k,1)
  133    continue
  134 continue
c
  135 do 137 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 136 k=1,l1
            cc(ido,j2-2,k) = ch(1,k,j)
            cc(1,j2-1,k) = ch(1,k,jc)
  136    continue
  137 continue
c
      if (ido .eq. 1) return
      if (nbd .lt. l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k=1,l1
            do 138 i=3,ido,2
               ic = idp2-i
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  138       continue
  139    continue
  140 continue
      return
c
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 143 i=3,ido,2
            ic = idp2-i
            do 142 k=1,l1
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  142       continue
  143    continue
  144 continue
c
      return
      end
      subroutine drfftb (n,r,wsave)
      double precision r(1), wsave(1)
c
      if (n .eq. 1) return
c
      call drftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
c
      return
      end
      subroutine drfftf (n,r,wsave)
      double precision r(1), wsave(1)

      if (n .eq. 1) return

      call drftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))

      return
      end
      subroutine drffti (n,wsave)
      double precision wsave(1)
c
      if (n .eq. 1) return
c
      call drfti1 (n,wsave(n+1),wsave(2*n+1))
c
      return
      end
      subroutine drftb1 (n,c,ch,wa,ifac)
      double precision c(1), ch(1), wa(1)
      integer ifac(1)
c
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call dradb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call dradb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
c
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call dradb2 (ido,l1,c,ch,wa(iw))
         go to 105
  104    call dradb2 (ido,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
c
  106    if (ip .ne. 3) go to 109
         ix2 = iw+ido
         if (na .ne. 0) go to 107
         call dradb3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call dradb3 (ido,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
c
  109    if (ip .ne. 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 110
         call dradb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call dradb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
c
  112    if (na .ne. 0) go to 113
         call dradbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call dradbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (ido .eq. 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue
c
      if (na .eq. 0) return
      do 117 i=1,n
         c(i) = ch(i)
  117 continue
c
      return
      end
      subroutine drftf1 (n,c,ch,wa,ifac)
      double precision c(1), ch(1), wa(1)
      integer ifac(1)
c
      nf = ifac(2)
      na = 1
      l2 = n
      iw = n
      do 111 k1=1,nf
         kh = nf-k1
         ip = ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip .ne. 4) go to 102
c
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call dradf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call dradf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
         go to 110
c
  102    if (ip .ne. 2) go to 104
         if (na .ne. 0) go to 103
         call dradf2 (ido,l1,c,ch,wa(iw))
         go to 110
  103    call dradf2 (ido,l1,ch,c,wa(iw))
         go to 110
c
  104    if (ip .ne. 3) go to 106
         ix2 = iw+ido
         if (na .ne. 0) go to 105
         call dradf3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 110
  105    call dradf3 (ido,l1,ch,c,wa(iw),wa(ix2))
         go to 110
c
  106    if (ip .ne. 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 107
         call dradf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call dradf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
c
  108    if (ido .eq. 1) na = 1-na
         if (na .ne. 0) go to 109
         call dradfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         na = 1
         go to 110
  109    call dradfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
         na = 0
c
  110    l2 = l1
  111 continue
c
      if (na .eq. 1) return
      do 112 i=1,n
         c(i) = ch(i)
  112 continue
c
      return
      end
      subroutine drfti1 (n,wa,ifac)
      double precision wa(1), arg, argh, argld, fi, tpi
      integer ifac(1), ntryh(4)
      data ntryh(1), ntryh(2), ntryh(3), ntryh(4) /4, 2, 3, 5/
      data tpi   /  6.2831853071 7958647692 5286766559 00577d0/
c
      nl = n
      nf = 0
      j = 0
c
  101 j = j+1
      if (j.le.4) ntry = ntryh(j)
      if (j.gt.4) ntry = ntry + 2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr.ne.0) go to 101
c
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
  107 if (nl .ne. 1) go to 104
      ifac(1) = n
      ifac(2) = nf
c
      argh = tpi/dfloat(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 110 k1=1,nfm1
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = dfloat(ld)*argh
            fi = 0.d0
            do 108 ii=3,ido,2
               i = i+2
               fi = fi+1.d0
               arg = fi*argld
               wa(i-1) = dcos(arg)
               wa(i) = dsin(arg)
  108       continue
            is = is+ido
  109    continue
c
         l1 = l2
  110 continue
c
      return
      end
