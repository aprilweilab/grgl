/*
// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton
*/

// Modified slightly to C++-ify it by Drew DeHaas:
// * use std::vector instead of malloc()
// * throw exceptions instead of exit(), which makes it friendlier for APIs.
// * change integer types to stdint/size_t

#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <vector>

double SNPHWE(int64_t obs_hets, int64_t obs_hom1, int64_t obs_hom2) {
    if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) {
        std::stringstream ssErr;
        ssErr << "HWE received a negative count: " << obs_hets << ", " << obs_hom1 << ", " << obs_hom2;
        throw std::runtime_error(ssErr.str());
    }

    int64_t obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
    int64_t obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

    int64_t rare_copies = 2 * obs_homr + obs_hets;
    int64_t genotypes = obs_hets + obs_homc + obs_homr;

    std::vector<double> het_probs(rare_copies + 1);

    for (size_t i = 0; i <= rare_copies; i++)
        het_probs[i] = 0.0;

    /* start at midpoint */
    int64_t mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

    /* check to ensure that midpoint and rare alleles have same parity */
    if ((rare_copies & 1) ^ (mid & 1))
        mid++;

    int64_t curr_hets = mid;
    int64_t curr_homr = (rare_copies - mid) / 2;
    int64_t curr_homc = genotypes - curr_hets - curr_homr;

    het_probs[mid] = 1.0;
    double sum = het_probs[mid];
    for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
        het_probs[curr_hets - 2] =
            het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
        sum += het_probs[curr_hets - 2];

        /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
        curr_homr++;
        curr_homc++;
    }

    curr_hets = mid;
    curr_homr = (rare_copies - mid) / 2;
    curr_homc = genotypes - curr_hets - curr_homr;
    for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
        het_probs[curr_hets + 2] =
            het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0));
        sum += het_probs[curr_hets + 2];

        /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
        curr_homr--;
        curr_homc--;
    }

    for (size_t i = 0; i <= rare_copies; i++)
        het_probs[i] /= sum;

    double p_hwe = 0.0;
    /*  p-value calculation for p_hwe  */
    for (size_t i = 0; i <= rare_copies; i++) {
        if (het_probs[i] > het_probs[obs_hets])
            continue;
        p_hwe += het_probs[i];
    }

    p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;
    return p_hwe;
}
