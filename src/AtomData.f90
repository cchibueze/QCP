!#######################################################################################################################!
    
    
    
    

!=======================================================================================================================!
!                                                                                                                       !
!                                                       ATOMIC DATA                                                     !
!                                                                                                                       !
!=======================================================================================================================!

    
    
    
!#######################################################################################################################!

    
    
    
    
    
    
    
    
    
!#######################################################################################################################!   
module atomic_data
implicit none
save
!=======================================================================================================================!
!                                                    MODULE VARIABLES DELCARATION                                       !
!=======================================================================================================================!

!to declare nucleus types with a string (list position in this case)
character (len=12) , dimension(1:100) :: atomname=(/'HYDROGEN    ', 'HELIUM      ', 'LITHIUM     ', 'BERYLLIUM   ', 'BORON       ', 'CARBON      ', 'NITROGEN    ', 'OXYGEN      ', 'FLUORINE    ', 'NEON        ', &            
                                                    'SODIUM      ', 'MAGNESIUM   ', 'ALUMINIUM   ', 'SILICON     ', 'PHOSPHORUS  ', 'SULFUR      ', 'CHLORINE    ', 'ARGON       ', 'POTASSIUM   ', 'CALCIUM     ', &
                                                    'SCANDIUM    ', 'TITANIUM    ', 'VANADIUM    ', 'CHROMIUM    ', 'MANGANESE   ', 'IRON        ', 'COBALT      ', 'NICKEL      ', 'COPPER      ', 'ZINC        ', & 
                                                    'GALLIUM     ', 'GERMANIUM   ', 'ARSENIC     ', 'SELENIUM    ', 'BROMINE     ', 'KRYPTON     ', 'RUBIDIUM    ', 'STRONTIUM   ', 'YTTRIUM     ', 'ZIRCONIUM   ', &
                                                    'NIOBIUM     ', 'MOLYBDENUM  ', 'TECHNETIUM  ', 'RUTHENIUM   ', 'RHODIUM     ', 'PALLADIUM   ', 'SILVER      ', 'CADMIUM     ', 'INDIUM      ', 'TIN         ', &
                                                    'ANTIMONY    ', 'TELLURIUM   ', 'IODINE      ', 'XENON       ', 'CESIUM      ', 'BARIUM      ', 'LANTHANUM   ', 'CERIUM      ', 'PRASEODYMIUM', 'NEODYMIUM   ', &
                                                    'PROMETHIUM  ', 'SAMARIUM    ', 'EUROPIUM    ', 'GADOLINIUM  ', 'TERBIUM     ', 'DYSPROSIUM  ', 'HOLMIUM     ', 'ERBIUM      ', 'THULIUM     ', 'YBERTIUM    ', &
                                                    'LUTETIUM    ', 'HAFNIUM     ', 'TANTALUM    ', 'TUNGSTEN    ', 'RHENIUM     ', 'OSMIUM      ', 'IRIDIUM     ', 'PLATINUM    ', 'GOLD        ', 'MERCURY     ', &
                                                    'THALLIUM    ', 'LEAD        ', 'BISMUTH     ', 'POLONIUM    ', 'ASTATINE    ', 'RADON       ', 'FRANCIUM    ', 'RADIUM      ', 'ACTINIUM    ', 'THORIUM     ', &  
                                                    'PROTACTINIUM', 'URANIUM     ', 'NEPTUNIUM   ', 'PLUTONIUM   ', 'AMERICIUM   ', 'CURIUM      ', 'BERKELIUM   ', 'CALIFORNIUM ', 'EINSTEINIUM ', 'FERMIUM     '  /)  

character (len=2) , dimension(1:100) :: atomsymbol=(/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
                                                      'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
                                                      'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
                                                      'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
                                                      'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
                                                      'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
                                                      'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
                                                      'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
                                                      'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
                                                      'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm'  /) !to declare nucleus types with a integer (list position in this case)

real (kind = 8) , dimension(1:100) :: atomic_mass=(/ &
    1.008, 4.002602, 6.94, 9.0121831, 10.81, 12.011, 14.007, 15.999, 18.998403163, 20.1797,             & ! 1–10
    22.98976928, 24.305, 26.9815385, 28.085, 30.973761998, 32.06, 35.45, 39.948, 39.0983, 40.078,       & ! 11–20
    44.955908, 47.867, 50.9415, 51.9961, 54.938044, 55.845, 58.933194, 58.6934, 63.546, 65.38,          & ! 21–30
    69.723, 72.63, 74.921595, 78.971, 79.904, 83.798, 85.4678, 87.62, 88.90584, 91.224,                 & ! 31–40
    92.90637, 95.95, 98.0, 101.07, 102.9055, 106.42, 107.8682, 112.414, 114.818, 118.71,                & ! 41–50
    121.76, 127.6, 126.90447, 131.293, 132.90545196, 137.327, 138.90547, 140.116, 140.90766, 144.242,   & ! 51–60
    145.0, 150.36, 151.964, 157.25, 158.92535, 162.5, 164.93033, 167.259, 168.93422, 173.054,           & ! 61–70
    174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084, 196.966569, 200.592,        & ! 71–80
    204.38, 207.2, 208.9804, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0377,                        & ! 81–90
    231.03588, 238.02891, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0 /)                       ! 91–100


real (kind = 8) , dimension(1:100) :: nuclear_charge !a.u.
integer , dimension(1:100) :: atomcount=0

contains
    subroutine set_nuclearcharge()
    integer :: i
    do i = 1,100
        nuclear_charge(i) = real(i)
    enddo
    end subroutine
end module atomic_data
!#######################################################################################################################!