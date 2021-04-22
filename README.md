# Models of cell-cell interactions in biological tissues

## Anita Youngblood - 160003606 - acry@st-andrews.ac.uk
## Supervisor: Jochen Kursawe

Code for MT5599 project on Delta/Notch interactions in biological tissues. 

Main models used are found in 
- TestDeltaNotchReporterLimiPrepattern.hpp 
(This uses Sprinzak et al.'s three-substance system of ODEs and incorporates the possibility of a prepattern. 
 Suggested initial conditions are included and can be commented in/out to choose the prepattern.)
- TestDeltaNotchReporterProtrusionLimi.hpp
(This uses Sprinzak et al. and Hunter et al.'s three-substance system of ODEs and incorporates the possibility of signalling through cellular protrusions.
 Protrusions are parameterised in DeltaNotchReporterProtrusionOdeSystemLimi.cpp by length, tip length, angle, and angular opening.)
- TestDeltaNotchReporterProtrusionLimiDivisionReporterDependent.hpp
(This uses Sprinzak et al.'s three-substance system of ODEs and incorporates the possibility of signalling through cellular protrusions as well as Notch-dependent cell
 division. Protrusions are parameterised in DeltaNotchReporterProtrusionOdeSystemLimi.cpp by length, tip length, angle, and angular opening.
 Division probability at each time-step is calculated based on Reporter expression using a Hill function with parameters specified in ReporterDependentCellCycleModel.cpp)

These can be implemented in the usual manner, for example using the Chaste docker (https://github.com/Chaste/chaste-docker)

