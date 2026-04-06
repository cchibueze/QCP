# QCP
Quantum Chemistry Package written by Chima Chibueze

Clone the repository and carry out the following compile steps:

1. Go into the QCP repository
2. Within the terminal: usr@pc % mkdir build
2. Within the terminal: usr@pc % cd ./build/
3. Within the terminal: usr@pc % cmake ../CMakeLists.txt
4. Within the terminal: usr@pc % make

Run a QCP calculation by:
1. Within the terminal: usr@pc % vi ./QCP/INPUT.dat        (make/alter input file)
2. Within the terminal: usr@pc % cd ./QCP/build/src/       (go to the program directory)
3. Within the terminal: usr@pc % ./my_exe                  (run the program)
4. Within the terminal: usr@pc % vi ../../OUTPUT.dat       (check the output)

What follows is the structure of INPUT.dat for an HF calculation in a minimal basis for H2
Of course, one can alter the xyz coordinates and the keywords for different calculation
One must copy past the lines between (thus not including) the quotation marks: " "

"

basis sto-3g

NGC 0
GOWF hf
maxiter 100
gradeps 0.0001

NR  0
maxiter 100
gradeps 0.001

SP  1

MP2 0

CID 0

CC 0
type SD

SRC 0

FC  0

CISES 0

THDF 0

Charge 0

Multiplicity 1

Units A

H    0.000000    0.000000    0.000000
H    0.000000    0.000000    0.740000

"