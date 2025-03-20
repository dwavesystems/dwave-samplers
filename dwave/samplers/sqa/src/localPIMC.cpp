/*
Copyright 2024 D-Wave

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

/*
   Path Integral Monte Carlo code for analysis of a finite temperature transverse field Ising model,
   Defined by partition function Z = Trace[exp(-invTemp H)], and scaled Hamiltonian
   H = 1/2 sum_{i,j} J_{ij} \sigma^z_i \sigma^z_j + \sum_i [h_i \sigma^z_i - \Gamma\sigma^x_i]

   Methods exploit either single qubit Swendsen-Wang updates, or multi-qubit Swendsen-Wang updates,
   in the latter case specifically for regular and independent 1d ferromagnetic subsequences. These
   match the methods explored in A King et al. https://arxiv.org/abs/1911.03446

   Authors: Jack Raymond, Stephen Face
   Copyright: D-Wave Systems
   License: Apache 2
   Last modification: April 1 2022

   See also README.md localPIMC.h and main.cpp
*/

#include <cstdint>
#include <math.h>
#include <vector>
#include <stdexcept>
#include "localPIMC.h"


using namespace std;


// Returns the energy of a given state and problem
// @param state a int8_t array containing the spin state to compute the energy of
// @param h vector of h or field value on each variable
// @param coupler_starts an int vector containing the variables of one side of
//        each coupler in the problem
// @param coupler_ends an int vector containing the variables of the other side 
//        of each coupler in the problem
// @param coupler_weights a double vector containing the weights of the 
//        couplers in the same order as coupler_starts and coupler_ends
// @return A double corresponding to the energy for `state` on the problem
//        defined by h and the couplers passed in
double get_state_energy(
    std::int8_t* state,
    const vector<double>& h,
    const vector<int>& coupler_starts,
    const vector<int>& coupler_ends,
    const vector<double>& coupler_weights
) {
    double energy = 0.0;
    // sum the energy due to local fields on variables
    for (unsigned int var = 0; var < h.size(); var++) {
        energy += state[var] * h[var];
    }
    // sum the energy due to coupling weights
    for (unsigned int c = 0; c < coupler_starts.size(); c++) {
        energy += state[coupler_starts[c]] * coupler_weights[c] * state[coupler_ends[c]];
    }
    return energy;
}


/* Public functions */

localPIMC::localPIMC(int Lperiodic, double invTempOverJ, double GammaOverJ, int initialCondition, int qubitsPerUpdate0,
                     int qubitsPerChain0, unsigned int seed) {
    // Constructor parameterized to match paper experiments:
    qubitsPerUpdate = qubitsPerUpdate0;
    qubitsPerChain = qubitsPerChain0;
    assert(qubitsPerChain >= qubitsPerUpdate);
    numVar = Lperiodic * (Lperiodic / 2 + 3) * qubitsPerChain;
    constructCouplingMatrix(Lperiodic, invTempOverJ);
    initializeWorldLines(initialCondition, Lperiodic, qubitsPerChain);
    const double Jchain = -1.8;
    invTempH = std::vector<double>(numVar, 0);
    invTempJchain = Jchain * invTempOverJ;
    invTempGamma = GammaOverJ * invTempOverJ;
    initPRNG(seed);

    invTemp = 1;
}
localPIMC::localPIMC(double Gamma, double Jchain, int qubitsPerUpdate0, int qubitsPerChain0,
                     std::vector<std::vector<int> > adjMat0, std::vector<std::vector<double> > invTempJ0,
                     std::vector<double> invTempH0, std::vector<int> classicalInitialCondition, unsigned int seed) {
    // Generic constructor
    qubitsPerUpdate = qubitsPerUpdate0;
    qubitsPerChain = qubitsPerChain0;
    assert(qubitsPerChain >= qubitsPerUpdate);
    numVar = adjMat0.size();
    adjMat.swap(adjMat0);
    invTempJ.swap(invTempJ0);
    invTempH.swap(invTempH0);
    invTempJchain = Jchain;
    invTempGamma = Gamma;
    breaks.resize(numVar);
    firstSlice.swap(classicalInitialCondition);
    initPRNG(seed);
    
    invTemp = 1;
};
void localPIMC::run(int nSweeps) {
    if (qubitsPerUpdate == 1) {
        std::uniform_int_distribution<> randomQubitIndex(0, numVar - 1);
        for (int sweepI = 0; sweepI < nSweeps * numVar; sweepI++) {
            qubitUpdate(randomQubitIndex(prng));
        }
    } else {
        int numChains = numVar / qubitsPerChain;
        std::uniform_int_distribution<> randomChainIndex(0, numChains - 1);
        for (int sweepI = 0; sweepI < nSweeps * numChains; sweepI++) {
            chainUpdate(randomChainIndex(prng));
        }
    }
}
void localPIMC::run(const std::vector<double> & HdField,
                    const std::vector<double> & HpField,
                    const int nSweepsPerField,
		    std::int8_t *statistics,
		    const int evaluateEvery) {
    if (qubitsPerUpdate == 1) {
        std::uniform_int_distribution<> randomQubitIndex(0, numVar - 1);
        for(int schedI = 0; schedI < HdField.size(); schedI++){
            invTempGamma = HdField[schedI];
            invTemp = HpField[schedI];
            for (int sweepI = 0; sweepI < nSweepsPerField * numVar; sweepI++)
                qubitUpdate(randomQubitIndex(prng));
	    if(evaluateEvery !=0 && schedI%evaluateEvery==0){
	        readSlice(statistics);
	        statistics+=numVar;
	    }
        }
    } else {
        int numChains = numVar / qubitsPerChain;
        std::uniform_int_distribution<> randomChainIndex(0, numChains - 1);
        for(int schedI = 0; schedI < HdField.size(); schedI++){
            invTempGamma = HdField[schedI];
            invTemp = HpField[schedI];
            for (int sweepI = 0; sweepI < nSweepsPerField * numVar; sweepI++) {
                chainUpdate(randomChainIndex(prng));
            }
	    if(evaluateEvery!=0 && schedI%evaluateEvery==0){
	        readSlice(statistics);
	        statistics+=numVar;
	    }
        }
    }
}


void localPIMC::reinitClassical(std::int8_t *state) {
  for(unsigned int i=0; i<firstSlice.size(); i++) {
    firstSlice[i] = state[i];
    breaks[i].resize(0);
  }
}
int localPIMC::reinitQuantum(std::int8_t *state, int *num_breaks, int *breaks_buffer) {
  int num_breaks_total = 0;
  for(unsigned int i=0; i<firstSlice.size(); i++) {
    firstSlice[i] = state[i];
    breaks[i].resize(num_breaks[i]);
    for(unsigned int j=0; j<firstSlice.size(); j++){
      breaks[i][j] = *breaks_buffer++;
    }
    num_breaks_total += num_breaks[i];
  }
  return(num_breaks_total);
}
void localPIMC::readSlice(std::int8_t *state) {
  for(unsigned int i=0; i<firstSlice.size(); i++) {
    state[i] = firstSlice[i];
  }
}
int localPIMC::readBreaks(int *num_breaks, int *breaks_buffer, int capacity) {
  int num_breaks_total=0;
  for(unsigned int i=0; i<breaks.size(); i++) {
    num_breaks[i] = breaks[i].size();
    if(capacity >= num_breaks[i]) {
        for(unsigned int j=0; j<breaks[i].size(); j++) {
            *breaks_buffer++ = breaks[i][j];
	}
	capacity -= num_breaks[i];
    }
    num_breaks_total += breaks[i].size();
  }
  return(num_breaks_total);
}

/* Private functions */

void localPIMC::constructCouplingMatrix(int Lperiodic, double invTemp0) {
    // Construct sparse coupling matrix for cylindrical lattices used in paper, either triangular or square
    // octagonal.
    assert(Lperiodic % 6 == 0);
    int Lopen = (Lperiodic / 2 + 3);
    assert(qubitsPerChain == 1 || qubitsPerChain == 4);
    int disp_i[] = {0, 1, 1}, disp_j[] = {1, 0, 1};
    int chainConnectionFrom[] = {(qubitsPerChain > 1 ? qubitsPerChain - 2 : 0), qubitsPerChain - 1, qubitsPerChain - 1},
        chainConnectionTo[] = {0, (qubitsPerChain > 1 ? 1 : 0), 0};
    invTempJ.resize(numVar);
    adjMat.resize(numVar);
    for (int i0 = 0; i0 < Lperiodic; i0++) {
        for (int j0 = 0; j0 < Lopen; j0++) {
            for (int couplerOrientationI = 0; couplerOrientationI < 3; couplerOrientationI++) {
                int i1 = (i0 + disp_i[couplerOrientationI]) % Lperiodic;
                int j1 = j0 + disp_j[couplerOrientationI];
                if (j1 < Lopen) {
                    int linearIndex0 = (i0 * Lopen + j0) * qubitsPerChain + chainConnectionFrom[couplerOrientationI];
                    int linearIndex1 = (i1 * Lopen + j1) * qubitsPerChain + chainConnectionTo[couplerOrientationI];
                    adjMat[linearIndex1].push_back(linearIndex0);
                    adjMat[linearIndex0].push_back(linearIndex1);

                    if ((j0 % (Lopen - 1)) || (j1 % (Lopen - 1))) {
                        invTempJ[linearIndex1].push_back(invTemp0);
                        invTempJ[linearIndex0].push_back(invTemp0);
                    } else {
                        // Cylindrical boundary condition:
                        invTempJ[linearIndex1].push_back(invTemp0 / 2);
                        invTempJ[linearIndex0].push_back(invTemp0 / 2);
                    }
                }
            }
        }
    }
    if (qubitsPerUpdate == 1 && qubitsPerChain > 1) {
        // Enumerate also couplings within chain if qubitUpdate instead of chainUpdate
        for (int n = 0; n < numVar; n += qubitsPerChain) {
            for (int k = 0; k < qubitsPerChain - 1; k++) {
                adjMat[n + k].push_back(n + k + 1);
                adjMat[n + k + 1].push_back(n + k);
                invTempJ[n + k].push_back(invTempJchain);
                invTempJ[n + k + 1].push_back(invTempJchain);
            }
        }
    }
}

void localPIMC::initializeWorldLines(int initialCondition, int Lperiodic, int qubitsPerChain) {
    // Initial state of Markov Chain, clockwise (1), counterclockwise (-1) or unwound (0)
    // Periodicity required in vertical direction (cylindrical lattice), dimension is multiple of 6
    assert(Lperiodic % 6 == 0 && Lperiodic >= 6);
    assert(initialCondition <= 1 && -1 <= initialCondition);
    int Lopen = 3 * (Lperiodic / 6 + 1);
    numVar = Lperiodic * Lopen * qubitsPerChain;
    // See paper description, 6 blocks, each block uses a different pseudoSpin orientation:
    int alignedMask[] = {1, 1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, -1};  // Check paper!
    std::vector<int> tripartiteClassification = makeTripartiteClassification(Lperiodic);
    int blockSize = numVar / (6 * qubitsPerChain);
    int blockState = 0;
    int nChains = Lperiodic * Lopen;
    firstSlice.resize(numVar);
    breaks.resize(numVar);
    for (int n = 0; n < nChains; n++) {
        breaks[n].resize(0);
        assert(tripartiteClassification.size() == nChains);
        assert(blockState * 3 + tripartiteClassification[n] < 18);
        int thisSpin = alignedMask[blockState * 3 + tripartiteClassification[n]];
        for (int k = 0; k < qubitsPerChain; k++) {
            firstSlice[qubitsPerChain * n + k] = thisSpin;
        }
        if (n % blockSize == blockSize - 1) {
            // ordered (initialCondition=0) leaves pseudoSpin unchanged, otherwise rotated
            blockState = (blockState + initialCondition + 6) % 6;
        }
    }
}
void localPIMC::addHToEffectiveField(std::vector<double>& effectiveFieldPerDomain,
                                     const std::vector<int>& allInterfaces, double H) const {
    H /= numTrotterSlices;
    effectiveFieldPerDomain[0] += (numTrotterSlices + allInterfaces.front() - allInterfaces.back()) * H;
    for (unsigned int interfaceI = 1; interfaceI < allInterfaces.size(); interfaceI++) {
        effectiveFieldPerDomain[interfaceI] += (allInterfaces[interfaceI] - allInterfaces[interfaceI - 1]) * H;
    }
}

void localPIMC::addHToEffectiveField(std::vector<double>& effectiveFieldPerDomain,
                                     const std::vector<int>& componentLabels, int componentOffset,
                                     const std::vector<int>& allInterfaces, double H) const {
    H /= numTrotterSlices;
    effectiveFieldPerDomain[componentLabels[componentOffset]] +=
        (numTrotterSlices + allInterfaces.front() - allInterfaces.back()) * H;
    for (unsigned int interfaceI = 1; interfaceI < allInterfaces.size(); interfaceI++) {
        effectiveFieldPerDomain[componentLabels[componentOffset + interfaceI]] +=
            (allInterfaces[interfaceI] - allInterfaces[interfaceI - 1]) * H;
    }
}

void localPIMC::addJToEffectiveField(std::vector<double>& effectiveFieldPerDomain,
                                     const std::vector<int>& allInterfaces, int neighbor, double Js) const {
    Js /= numTrotterSlices;
    std::vector<int> allInterfacesPair(allInterfaces.size() + breaks[neighbor].size());
    std::merge(breaks[neighbor].begin(), breaks[neighbor].end(), allInterfaces.begin(), allInterfaces.end(),
               allInterfacesPair.begin());
    effectiveFieldPerDomain[0] += (numTrotterSlices + allInterfacesPair.front() - allInterfacesPair.back()) * Js;
    unsigned int interfaceI = 0;
    unsigned int interfaceAllI = 0;
    for (; interfaceAllI + 1 < allInterfacesPair.size(); interfaceAllI++) {
        if (allInterfacesPair[interfaceAllI] == allInterfaces[interfaceI]) {
            if (++interfaceI == allInterfaces.size()) break;
        } else {
            Js *= -1;
        }
        effectiveFieldPerDomain[interfaceI] +=
            (allInterfacesPair[interfaceAllI + 1] - allInterfacesPair[interfaceAllI]) * Js;
    }
    for (; interfaceAllI + 1 < allInterfacesPair.size(); interfaceAllI++) {
        effectiveFieldPerDomain[0] += (allInterfacesPair[interfaceAllI + 1] - allInterfacesPair[interfaceAllI]) * Js;
        Js *= -1;
    }
}
void localPIMC::addJToEffectiveField(std::vector<double>& effectiveFieldPerDomain,
                                     const std::vector<int>& componentLabels, int componentOffset,
                                     const std::vector<int>& allInterfaces, int neighbor, double Js) const {
    Js /= numTrotterSlices;
    std::vector<int> allInterfacesPair(allInterfaces.size() + breaks[neighbor].size());
    std::merge(breaks[neighbor].begin(), breaks[neighbor].end(), allInterfaces.begin(), allInterfaces.end(),
               allInterfacesPair.begin());
    effectiveFieldPerDomain[componentLabels[componentOffset]] +=
        (numTrotterSlices + allInterfacesPair.front() - allInterfacesPair.back()) * Js;
    unsigned int interfaceI = 0;
    unsigned int interfaceAllI = 0;
    for (; interfaceAllI + 1 < allInterfacesPair.size(); interfaceAllI++) {
        if (allInterfacesPair[interfaceAllI] == allInterfaces[interfaceI]) {
            if (++interfaceI == allInterfaces.size()) break;
        } else {
            Js *= -1;
        }
        effectiveFieldPerDomain[componentLabels[componentOffset + interfaceI]] +=
            (allInterfacesPair[interfaceAllI + 1] - allInterfacesPair[interfaceAllI]) * Js;
    }
    for (; interfaceAllI + 1 < allInterfacesPair.size(); interfaceAllI++) {
        effectiveFieldPerDomain[componentLabels[componentOffset]] +=
            (allInterfacesPair[interfaceAllI + 1] - allInterfacesPair[interfaceAllI]) * Js;
        Js *= -1;
    }
}

void localPIMC::qubitUpdate(const int sp) {
    // Update 1 qubit system H = \nu_i \sigma^z_i + invTempGamma \sigma^x_i, where \nu_i are determined from
    // neighboring qubit states
    std::vector<int> breakProposals = makeBreakProposals();
    std::vector<int> allInterfaces;
    if (breakProposals.size() + breaks[sp].size() > 1) {
        allInterfaces.resize(breakProposals.size() + breaks[sp].size());
        std::merge(breakProposals.begin(), breakProposals.end(), breaks[sp].begin(), breaks[sp].end(),
                   allInterfaces.begin());
    } else
        allInterfaces.push_back(numTrotterSlices);  // Book keeping, less branching

    int nDom = allInterfaces.size();
    std::vector<double> effectiveFieldPerDomain(nDom, 0);
    addHToEffectiveField(effectiveFieldPerDomain, allInterfaces, invTempH[sp]);
    for (unsigned int neighI = 0; neighI < adjMat[sp].size(); neighI++) {
        unsigned int neighbor = adjMat[sp][neighI];
        addJToEffectiveField(effectiveFieldPerDomain, allInterfaces, neighbor,
                             firstSlice[neighbor] * invTempJ[sp][neighI]);
    }
    // Sample spin on boundary spanning domain
    int sValue = GibbsSamplePM1(effectiveFieldPerDomain[0]);
    firstSlice[sp] = sValue;
    breaks[sp].resize(0);
    // Sample states
    for (unsigned int interfaceI = 1; interfaceI < allInterfaces.size(); interfaceI++) {
        // Sample spin on subsequent domains, if flipped record interface:
        if (sValue * GibbsSamplePM1(effectiveFieldPerDomain[interfaceI]) != 1) {
            sValue *= -1;
            breaks[sp].push_back(allInterfaces[interfaceI - 1]);
        }
    }
    if (sValue != firstSlice[sp]) breaks[sp].push_back(allInterfaces.back());
}
//
void localPIMC::chainUpdate(const int sp) {
    // Update multi-qubits system H = invTempJchain \sum_i \sigma^z_i \sigma^z_{i+1} + \sum_i [invTemp_i \sigma^x_i
    // + \nu_i \sigma^Z_i],
    // by Swendsen-Wang where \nu_i are determined from neighboring qubit states
    std::vector<std::vector<int> > allInterfacesEveryChain(qubitsPerChain);
    std::vector<int> nDomains(qubitsPerChain, 0);
    int nDomTotal = 0;
    for (int chainI = 0; chainI < qubitsPerChain; chainI++) {
        int qubitI = qubitsPerChain * sp + chainI;
        std::vector<int> breakProposals = makeBreakProposals();
        // Joint list of interfaces with
        if (breakProposals.size() + breaks[qubitI].size() > 1) {
            allInterfacesEveryChain[chainI].resize(breakProposals.size() + breaks[qubitI].size());
            std::merge(breakProposals.begin(), breakProposals.end(), breaks[qubitI].begin(), breaks[qubitI].end(),
                       allInterfacesEveryChain[chainI].begin());
        } else {
            allInterfacesEveryChain[chainI].push_back(numTrotterSlices);  // Book keeping, less branching later
        }
        nDomains[chainI] = nDomTotal;
        nDomTotal += allInterfacesEveryChain[chainI].size();
    }
    // Build connectivity of qubit level domains, by Swendsen-Wang rule
    std::vector<std::vector<int> > domainGraph(nDomTotal);
    for (int chainI = 0; chainI < qubitsPerChain - 1; chainI++) {
        makeDomainGraph(nDomains[chainI], nDomains[chainI + 1], sp, chainI, allInterfacesEveryChain, domainGraph);
    }
    // Find components
    int nComponents = 0;
    std::vector<int> componentLabels(nDomTotal, -1);
    for (unsigned int root = 0; root < componentLabels.size(); root++) {
        if (componentLabels[root] == -1) depthFirstDomainAssignment(domainGraph, componentLabels, nComponents++, root);
    }
    // Consolidate external fields
    std::vector<double> effectiveFieldAll(nComponents, 0);
    for (int chainI = 0; chainI < qubitsPerChain; chainI++) {
        int qubitI = qubitsPerChain * sp + chainI;
        addHToEffectiveField(effectiveFieldAll, componentLabels, nDomains[chainI], allInterfacesEveryChain[chainI],
                             invTempH[qubitI]);  // OPTIONAL FOR LATTICE CODE:
        for (unsigned int neighI = 0; neighI < adjMat[qubitI].size(); neighI++) {
            int neighbor = adjMat[qubitI][neighI];
            addJToEffectiveField(effectiveFieldAll, componentLabels, nDomains[chainI], allInterfacesEveryChain[chainI],
                                 neighbor, firstSlice[neighbor] * invTempJ[qubitI][neighI]);
        }
    }
    // sample domain spin according to effective fields
    std::vector<int> sValues(nComponents, 0);
    for (int componentLabel = 0; componentLabel < nComponents; componentLabel++) {
        sValues[componentLabel] = GibbsSamplePM1(effectiveFieldAll[componentLabel]);
    }
    // Map back assignment to each qubit and record interfaces
    for (int chainI = 0; chainI < qubitsPerChain; chainI++) {
        int qubitI = qubitsPerChain * sp + chainI;
        int dom = componentLabels[nDomains[chainI]];
        firstSlice[qubitI] = sValues[dom];
        breaks[qubitI].resize(0);
        int sValue = firstSlice[qubitI];
        for (unsigned int domainI = 1; domainI < allInterfacesEveryChain[chainI].size(); domainI++) {
            dom = componentLabels[nDomains[chainI] + domainI];
            if (sValue * sValues[dom] != 1) {
                sValue *= -1;
                breaks[qubitI].push_back(allInterfacesEveryChain[chainI][domainI - 1]);
            }
        }
        if (sValue != firstSlice[qubitI]) {
            breaks[qubitI].push_back(allInterfacesEveryChain[chainI].back());
        }
    }
}

void localPIMC::depthFirstDomainAssignment(const std::vector<std::vector<int> >& domainGraph,
                                           std::vector<int>& componentLabels, int componentLabel, int root) const {
    componentLabels[root] = componentLabel;
    for (unsigned int leaf = 0; leaf < domainGraph[root].size(); leaf++)
        if (componentLabels[domainGraph[root][leaf]] == -1)
            depthFirstDomainAssignment(domainGraph, componentLabels, componentLabel, domainGraph[root][leaf]);
}

void localPIMC::makeDomainGraph(int zeroChainIndex, int firstChainIndex, int sp, int chainI,
                                const std::vector<std::vector<int> >& allInterfacesEveryChain,
                                std::vector<std::vector<int> >& domainGraph) const {
    // If unfrustrated, attempt merge as function of interface size:
    std::uniform_real_distribution<double> probability(
        0.0, 1.0);  // Check if this is slow.. norm by max value would be good enough.
    int qubitI = sp * qubitsPerChain + chainI;
    int s0s1 = firstSlice[qubitI] * firstSlice[qubitI + 1];
    int nProposedInterfaces = allInterfacesEveryChain[chainI].size() + allInterfacesEveryChain[chainI + 1].size();
    std::vector<int> allInterfacesPair(nProposedInterfaces);
    std::merge(allInterfacesEveryChain[chainI].begin(), allInterfacesEveryChain[chainI].end(),
               allInterfacesEveryChain[chainI + 1].begin(), allInterfacesEveryChain[chainI + 1].end(),
               allInterfacesPair.begin());
    if (allInterfacesPair[0] == numTrotterSlices) {
        // No non-trivial domains, simple join:
        if (1 == s0s1 && pNotJoin(numTrotterSlices) < probability(prng)) {
            domainGraph[zeroChainIndex].push_back(firstChainIndex);
            domainGraph[firstChainIndex].push_back(zeroChainIndex);
        }
    } else {
        // At most one book keeping domain,remove:
        if (allInterfacesPair.back() == numTrotterSlices) {
            // Remove unnecessary book keeping values
            allInterfacesPair.pop_back();
        }
        // Special handling of boundary spanning domain (periodic in imaginary time):
        if (1 == s0s1 &&
            pNotJoin(numTrotterSlices - allInterfacesPair.back() + allInterfacesPair.front()) < probability(prng)) {
            domainGraph[zeroChainIndex].push_back(firstChainIndex);
            domainGraph[firstChainIndex].push_back(zeroChainIndex);
        }

        unsigned int posProposal0 = 0, posExisting0 = 0, posProposal1 = 0, posExisting1 = 0,
            chain0valid = (allInterfacesEveryChain[chainI][0] != numTrotterSlices);
        unsigned int interfaceI = 0;
        while (++interfaceI < allInterfacesPair.size()) {
            if (chain0valid && allInterfacesPair[interfaceI - 1] == allInterfacesEveryChain[chainI][posProposal0]) {
                if (allInterfacesEveryChain[chainI].size() == ++posProposal0) {
                    chain0valid = 0;
                    posProposal0 = 0;
                }
                if (breaks[qubitI].size() > posExisting0 &&
                    allInterfacesPair[interfaceI - 1] == breaks[qubitI][posExisting0]) {
                    s0s1 = s0s1 * (-1);
                    posExisting0++;  // Part of existing boundary
                }
            } else {
                if (allInterfacesEveryChain[chainI + 1].size() == ++posProposal1) {
                    posProposal1 = 0;
                }
                if (breaks[qubitI + 1].size() > posExisting1 &&
                    allInterfacesPair[interfaceI - 1] == breaks[qubitI + 1][posExisting1]) {
                    s0s1 = s0s1 * (-1);
                    posExisting1++;  // Part of existing boundary
                }
            }
            if (s0s1 == 1 &&
                pNotJoin(allInterfacesPair[interfaceI] - allInterfacesPair[interfaceI - 1]) < probability(prng)) {
                domainGraph[zeroChainIndex + posProposal0].push_back(firstChainIndex + posProposal1);
                domainGraph[firstChainIndex + posProposal1].push_back(zeroChainIndex + posProposal0);
            }
        }
    }
}

std::vector<int> localPIMC::makeBreakProposals() const {
    // Sampling of conditional distribution, exploiting continuous limit invTempGamma << numTrotterSlices:
    std::vector<int> breakProposals;
    std::uniform_real_distribution<double> probability(0.0, 1.0);
    double nInterfacesScale = numTrotterSlices / invTempGamma;
    double position = -nInterfacesScale * log(probability(prng));
    while (position < numTrotterSlices) {
        breakProposals.push_back(int(position));
        position = -nInterfacesScale * log(probability(prng)) + breakProposals.back() + 1;
    }
    return breakProposals;
}

std::vector<int> localPIMC::makeTripartiteClassification(int Lperiodic) const {
    int Lopen = 3 * (Lperiodic / 6 + 1);
    std::vector<int> tripartiteClassification(Lperiodic * Lopen);
    for (int i = 0; i < Lperiodic; i++) {
        for (int j = 0; j < Lopen; j++) {
            tripartiteClassification[i * Lopen + j] = (i + j) % 3;
        }
    }
    return tripartiteClassification;
}

int localPIMC::GibbsSamplePM1(const double effectiveField) const {
    // Sample from P(s) ~ exp( - effField *s) = (1 - s tanh(effField))/2, effField is the energy associated to state
    // +1 scaled to beta=1.
    std::uniform_real_distribution<double> magThresh(-1.0, 1.0);
    return tanh(invTemp*effectiveField) > magThresh(prng) ? -1 : 1;
}

double localPIMC::pNotJoin(int nOverlaps) const { return exp((2 * invTempJchain * nOverlaps * invTemp) / numTrotterSlices); }

void localPIMC::initPRNG(unsigned int seed) const {
    if (seed) {
        prng = std::mt19937(seed);
    } else {
        // seed = 0 signals use of random device
        std::random_device r;
        std::seed_seq seedSeq{r(), r(), r(), r(), r(), r(), r(), r()};
        prng = std::mt19937(seedSeq);
    }
}


// Perform simulated annealing on a general problem
// @param states a std::int8_t array of size num_samples * number of variables in the
//        problem. Will be overwritten by this function as samples are filled
//        in. The initial state of the samples are used to seed the simulated
//        annealing runs.
// @param energies a double array of size num_samples. Will be overwritten by
//        this function as energies are filled in.
// @param num_samples the number of samples to get.
// @param h vector of h or field value on each variable
// @param coupler_starts an int vector containing the variables of one side of
//        each coupler in the problem
// @param coupler_ends an int vector containing the variables of the other side 
//        of each coupler in the problem
// @param coupler_weights a double vector containing the weights of the couplers
//        in the same order as coupler_starts and coupler_ends
// @param sweeps_per_beta The number of sweeps to perform at each field value.
//        Total number of sweeps is `sweeps_per_beta` * length of
//        `Hp_field`.
// @param Hp_field A list of the longitudinal field values to run `sweeps_per_beta`
//        sweeps at.
// @param Hd_field A list of the transverse field values to run `sweeps_per_beta`
//        sweeps at.
// @param Gamma Transverse field per qubit. 
// @param chain_coupler_strength The coupling strength between qubits on the
//        same chain.
// @param qubits_per_chain Number of qubits per chain.
// @param qubits_per_update Number of qubits addressed per local (in space) update.
// @param interrupt_callback A function that is invoked between each run of simulated annealing
//        if the function returns True then it will stop running.
// @param interrupt_function A pointer to contents that are passed to interrupt_callback.
// @return the number of samples taken. If no interrupt occured, will equal num_samples.

int general_simulated_annealing(
    std::int8_t* states,
    double* energies,
    const bool project_inputs,
    const bool project_outputs,
    int *num_breaks,
    int *breaks_in,
    int *breaks_buffer_out,
    int breaks_buffer_size,
    const int num_samples,
    const vector<double> h,
    const vector<int> coupler_starts,
    const vector<int> coupler_ends,
    const vector<double> coupler_weights,
    const int sweeps_per_beta,
    const vector<double> HpField,
    const vector<double> HdField,
    const double Gamma, 
    const double chain_coupler_strength,
    const int qubits_per_chain,
    const int qubits_per_update,
    const unsigned int seed,
    std::int8_t* statistics,
    const int schedule_sample_interval,
    callback interrupt_callback,
    void * const interrupt_function
) {
    // TO DO 
    // assert len(states) == num_samples*num_vars*sizeof(std::int8_t)
    // assert len(coupler_starts) == len(coupler_ends) == len(coupler_weights)
    // assert max(coupler_starts + coupler_ends) < num_vars
    
    // the number of variables in the problem
    const int num_vars = h.size();
    if (!((coupler_starts.size() == coupler_ends.size()) &&
                (coupler_starts.size() == coupler_weights.size()))) {
        throw runtime_error("coupler vectors have mismatched lengths");
    }
    
    // neighbors is a vector of vectors, such that neighbors[i][j] is the jth
    // neighbor of variable i
    vector<vector<int>> neighbors(num_vars);
    // neighbour_couplings is another vector of vectors with the same structure
    // except neighbour_couplings[i][j] is the weight on the coupling between i
    // and its jth neighbor
    vector<vector<double>> neighbour_couplings(num_vars);

    // build the neighbors, and neighbour_couplings vectors by
    // iterating over the inputted coupler vectors
    for (unsigned int cplr = 0; cplr < coupler_starts.size(); cplr++) {
        int u = coupler_starts[cplr];
        int v = coupler_ends[cplr];

        if ((u < 0) || (v < 0) || (u >= num_vars) || (v >= num_vars)) {
            throw runtime_error("coupler indexes contain an invalid variable");
        }

        // add v to u's neighbors list and vice versa
        neighbors[u].push_back(v);
        neighbors[v].push_back(u);
        // add the weights
        neighbour_couplings[u].push_back(coupler_weights[cplr]);
        neighbour_couplings[v].push_back(coupler_weights[cplr]);
    }

    //Would make sense to build once here and then reuse. However,
    //PRNG might be issue if this loop is parallelized. For some
    //reason error is raised when passing pimc as object.

    std::vector<int> classical_initial_condition=std::vector<int>(num_vars); //copy-in for this format.
    localPIMC pimc(Gamma, chain_coupler_strength, qubits_per_update, qubits_per_chain, 
                   neighbors, neighbour_couplings,
                   h, classical_initial_condition, seed);
    
    // get the simulated annealing samples
    int sample_index = 0;
    int *p_breaks_in = breaks_in;
    int *p_breaks_buffer_out = breaks_buffer_out;
    int stats_per_anneal = schedule_sample_interval ?
      (1 + (HdField.size()-1)/schedule_sample_interval)*num_vars : 0;
    while (sample_index < num_samples) {
        // states is a giant spin array that will hold the resulting states for
        // all the samples, so we need to get the location inside that vector
        // where we will store the sample for this sample_index
        std::int8_t *p_states = states + sample_index*num_vars;
        int *p_num_breaks = num_breaks + sample_index*num_vars;
        if(project_inputs) {
          pimc.reinitClassical(p_states);
        } else {
          int num_breaks = pimc.reinitQuantum(p_states, p_num_breaks, p_breaks_in);
          p_breaks_in += num_breaks;
        }
        pimc.run(HdField, HpField, sweeps_per_beta,
		 statistics + sample_index*stats_per_anneal, schedule_sample_interval); //Add sample_every
        pimc.readSlice(p_states);
        // compute the energy of the sample and store it in `energies`
        energies[sample_index] = get_state_energy(p_states, h, coupler_starts, 
                                                coupler_ends, coupler_weights);
        
        if(!project_outputs) {
	  int num_breaks = pimc.readBreaks(p_num_breaks, p_breaks_buffer_out, breaks_buffer_size - (p_breaks_buffer_out-breaks_buffer_out));
	  p_breaks_buffer_out += num_breaks;
	}
        sample_index++;

        // if interrupt_function returns true, stop sampling
        if (interrupt_function && interrupt_callback(interrupt_function)) break;
    }

    // return the number of samples we actually took
    return sample_index;
}
