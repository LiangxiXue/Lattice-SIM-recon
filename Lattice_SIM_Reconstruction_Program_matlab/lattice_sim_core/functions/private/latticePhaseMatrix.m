function W = latticePhaseMatrix()
%LATTICEPHASEMATRIX Return the five-frame Lattice-SIM demodulation matrix.

phasePairs = [
    0,       0
    0,       2*pi/3
    0,       4*pi/3
    2*pi/3,  0
    4*pi/3,  2*pi/3
];

W = zeros(5, 5);
for idx = 1:5
    phiS = phasePairs(idx, 1);
    phiT = phasePairs(idx, 2);
    W(idx, :) = [1, exp(1i * phiS), exp(-1i * phiS), ...
        exp(1i * phiT), exp(-1i * phiT)];
end
end
