module conversion
  implicit none
  doubleprecision, parameter :: e_SI=1.6021766D-19, kB_SI=1.380649D-23,kB_au=3.166811563d0*(10**(-6)), KB_ev = 8.617333d-5, h_SI=6.62607015D-34, h_ev=4.135667696D-15,&
       hbar_SI=1.054571817D-34, hbar_ev=6.582119569D-16, epsilon0_SI= 8.8541878128D-12, me_SI=9.1093897D-31, c_SI = 2.99892D8, c_au = 137.d0, epsilon0_ev=55.26349406D6
  
contains
  
  
  doubleprecision function Jtoev(num)
    doubleprecision num
    jtoev=num*6.242D18
    return 
  end function Jtoev
  doubleprecision function evtoJ(num)
    doubleprecision num
    evtoj=num*1.6022D-19
    return 
  end function evtoJ
  doubleprecision function jtohartree(num)
    doubleprecision num
    jtohartree=num*2.294D17
    return 
  end function jtohartree
  doubleprecision function evtohartree(num)
    doubleprecision num
    evtohartree= num*0.0367493D0
    return 
  end function evtohartree
  doubleprecision function hartreetoev(num)
    doubleprecision num
    hartreetoev=num*27.2114D0
    return 
  end function hartreetoev
  doubleprecision function hartreetoj(num)
    doubleprecision num
    hartreetoj=num*4.3597D-18
      return 
    end function hartreetoj
  doubleprecision function Dbtoau(num)
    doubleprecision num
    Dbtoau=num*0.393456D0
    return 
  end function Dbtoau
  doubleprecision function autoDb(num)
    doubleprecision num
    autoDb=num/0.393456D0
    return 
  end function autoDb
  doubleprecision function cm1toev(num)
    doubleprecision num
    cm1toev=num/8065.5D0
    return 
  end function cm1toev
  doubleprecision function evtocm1(num)
    doubleprecision num
    evtocm1=num*8065.5D0
    return 
  end function evtocm1
  doubleprecision function DbtoCm(num)
    doubleprecision num
    DbtoCm=num*3.336D-30
    return
  end function DbtoCm
  doubleprecision function CmtoDb(num)
    doubleprecision num
    CmtoDb=num/(3.336D-30)
    return
  end function CmtoDb
end module conversion
