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

In the QCP/examples/ folder one finds examples of input files.
Their content can be copied into INPUT.dat
The most straightforward example is h2_hf.txt which contains the input structure for an HF calculation in a minimal basis for H2
The other one is a similar calculation, extended with a CCSD calculation.

For a detailed description of calculation features, see QCP/documentation/documentation.txt

