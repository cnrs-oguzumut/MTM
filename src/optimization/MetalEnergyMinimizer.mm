#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <simd/simd.h>

#include "../../include/optimization/MetalEnergyMinimizer.h"
#include "../../include/optimization/LatticeOptimizer.h"
#include "../../include/lattice_energy/Strain_Energy_LatticeCalculator.h"
#include "../../include/lattice_energy/EnergyFunctions.h"

#include <iostream>
#include <vector>
#include <cmath>

namespace metal_gpu {

struct GPUElement {
    int solver_idx[3];
    simd_float2 fixed_pos_or_trans[3];
    simd_float2 dNdX[3];
    float area;
    float padding; // 16-byte alignment
};

struct GPUParams {
    int n_elements;
    int n_vars;
    float beta;
    float K;
    float burgers;
    float zero_energy;
    float normalisation;
    float padding;
};

static const char* metal_shader_source = R"(
#include <metal_stdlib>
using namespace metal;

struct GPUElement {
    int solver_idx[3];
    float2 fixed_pos_or_trans[3];
    float2 dNdX[3];
    float area;
    float padding;
};

struct GPUParams {
    int n_elements;
    int n_vars;
    float beta;
    float K;
    float burgers;
    float zero_energy;
    float normalisation;
    float padding;
};

inline float calc_phi(float c11, float c22, float c12, float beta, float K, float burgers) {
    float detC = c11 * c22 - c12 * c12;
    float energy = -K * (log(detC / burgers) - detC / burgers) +
        beta * (
            pow(1.0f / pow(detC, 1.5f) * pow(c11 - 4.0f * c12 + c22, 3.0f) * (1.0f / 9.0f) -
                pow(c11 - c22, 2.0f) * (1.0f / pow(detC, 1.5f)) * (c11 - 4.0f * c12 + c22), 2.0f) * 9.46969696969697e-4f
            - pow((pow(c11 - c22, 2.0f) * 0.25f) / detC + (pow(c11 - 4.0f * c12 + c22, 2.0f) * (1.0f / 12.0f)) / detC, 3.0f) * (41.0f / 99.0f)
            + ((pow(c11 - c22, 2.0f) * 0.25f) / detC + (pow(c11 - 4.0f * c12 + c22, 2.0f) * (1.0f / 12.0f)) / detC) * (1.0f / (detC * detC)) * pow(c11 - c12 + c22, 4.0f) * (1.0f / 81.0f)
            - (1.0f / pow(detC, 1.5f) * pow(c11 - 4.0f * c12 + c22, 3.0f) * (1.0f / 9.0f) - pow(c11 - c22, 2.0f) * (1.0f / pow(detC, 1.5f)) * (c11 - 4.0f * c12 + c22)) *
              ((pow(c11 - c22, 2.0f) * 0.25f) / detC + (pow(c11 - 4.0f * c12 + c22, 2.0f) * (1.0f / 12.0f)) / detC) * (1.0f / sqrt(detC)) * (c11 - c12 + c22) * (7.0f / 198.0f)
        )
        + pow(1.0f / pow(detC, 1.5f) * pow(c11 - 4.0f * c12 + c22, 3.0f) * (1.0f / 9.0f) - pow(c11 - c22, 2.0f) * (1.0f / pow(detC, 1.5f)) * (c11 - 4.0f * c12 + c22), 2.0f) * (17.0f / 528.0f)
        + pow((pow(c11 - c22, 2.0f) * 0.25f) / detC + (pow(c11 - 4.0f * c12 + c22, 2.0f) * (1.0f / 12.0f)) / detC, 3.0f) * (4.0f / 11.0f)
        - (1.0f / pow(detC, 1.5f) * pow(c11 - 4.0f * c12 + c22, 3.0f) * (1.0f / 9.0f) - pow(c11 - c22, 2.0f) * (1.0f / pow(detC, 1.5f)) * (c11 - 4.0f * c12 + c22)) * (1.0f / pow(detC, 1.5f)) * pow(c11 - c12 + c22, 3.0f) * (1.0f / 27.0f)
        + (1.0f / pow(detC, 1.5f) * pow(c11 - 4.0f * c12 + c22, 3.0f) * (1.0f / 9.0f) - pow(c11 - c22, 2.0f) * (1.0f / pow(detC, 1.5f)) * (c11 - 4.0f * c12 + c22)) *
          ((pow(c11 - c22, 2.0f) * 0.25f) / detC + (pow(c11 - 4.0f * c12 + c22, 2.0f) * (1.0f / 12.0f)) / detC) * (1.0f / sqrt(detC)) * (c11 - c12 + c22) * (8.0f / 33.0f);
    return energy;
}

inline float3 calc_dphi(float c11, float c22, float c12, float beta, float K, float burgers) {
    float detC = c11 * c22 - c12 * c12;
    float detC_1_5 = pow(detC, 1.5f);
    float detC_2_5 = pow(detC, 2.5f);
    float detC_sq = detC * detC;
    float sqrt_detC = sqrt(detC);

    float termA = (pow(c11 - c22, 2.0f) * 0.25f) / detC + (pow(c11 - 4.0f * c12 + c22, 2.0f) * (1.0f / 12.0f)) / detC;
    float termB = 1.0f / detC_1_5 * pow(c11 - 4.0f * c12 + c22, 3.0f) * (1.0f / 9.0f) - pow(c11 - c22, 2.0f) * (1.0f / detC_1_5) * (c11 - 4.0f * c12 + c22);

    float dPhi00 = -beta*(pow(termA, 2.0f)*((c11*(1.0f/6.0f)-c12*(2.0f/3.0f)+c22*(1.0f/6.0f))/detC+(c11*0.5f-c22*0.5f)/detC-c22*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/12.0f)-c22*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.25f)*(41.0f/33.0f)
        + termB*(pow(c11-c22, 2.0f)/detC_1_5 - (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/3.0f) + c22*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/6.0f) + (1.0f/detC_1_5)*(2.0f*c11-2.0f*c22)*(c11-4.0f*c12+c22) - c22*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*1.5f)*(1.0f/528.0f)
        - termA*(1.0f/detC_sq)*pow(c11-c12+c22, 3.0f)*(4.0f/81.0f)
        - (1.0f/detC_sq)*pow(c11-c12+c22, 4.0f)*((c11*(1.0f/6.0f)-c12*(2.0f/3.0f)+c22*(1.0f/6.0f))/detC+(c11*0.5f-c22*0.5f)/detC-c22*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/12.0f)-c22*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.25f)*(1.0f/81.0f)
        + termB*termA*(1.0f/sqrt_detC)*(7.0f/198.0f)
        + c22*termA*(1.0f/pow(detC, 3.0f))*pow(c11-c12+c22, 4.0f)*(2.0f/81.0f)
        - termA*(1.0f/sqrt_detC)*(c11-c12+c22)*(pow(c11-c22, 2.0f)/detC_1_5 - (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/3.0f) + c22*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/6.0f) + (1.0f/detC_1_5)*(2.0f*c11-2.0f*c22)*(c11-4.0f*c12+c22) - c22*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*1.5f)*(7.0f/198.0f)
        + termB*(1.0f/sqrt_detC)*(c11-c12+c22)*((c11*(1.0f/6.0f)-c12*(2.0f/3.0f)+c22*(1.0f/6.0f))/detC+(c11*0.5f-c22*0.5f)/detC-c22*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/12.0f)-c22*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.25f)*(7.0f/198.0f)
        - c22*termB*termA*(1.0f/detC_1_5)*(c11-c12+c22)*(7.0f/396.0f))
        + K*(c22/burgers - c22/detC)
        + pow(termA, 2.0f)*((c11*(1.0f/6.0f)-c12*(2.0f/3.0f)+c22*(1.0f/6.0f))/detC+(c11*0.5f-c22*0.5f)/detC-c22*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/12.0f)-c22*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.25f)*(12.0f/11.0f)
        - termB*(pow(c11-c22, 2.0f)/detC_1_5 - (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/3.0f) + c22*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/6.0f) + (1.0f/detC_1_5)*(2.0f*c11-2.0f*c22)*(c11-4.0f*c12+c22) - c22*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*1.5f)*(17.0f/264.0f)
        + (1.0f/detC_1_5)*pow(c11-c12+c22, 3.0f)*(pow(c11-c22, 2.0f)/detC_1_5 - (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/3.0f) + c22*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/6.0f) + (1.0f/detC_1_5)*(2.0f*c11-2.0f*c22)*(c11-4.0f*c12+c22) - c22*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*1.5f)*(1.0f/27.0f)
        - termB*(1.0f/detC_1_5)*pow(c11-c12+c22, 2.0f)*(1.0f/9.0f)
        + termB*termA*(1.0f/sqrt_detC)*(8.0f/33.0f)
        + c22*termB*(1.0f/detC_2_5)*pow(c11-c12+c22, 3.0f)*(1.0f/18.0f)
        - termA*(1.0f/sqrt_detC)*(c11-c12+c22)*(pow(c11-c22, 2.0f)/detC_1_5 - (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/3.0f) + c22*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/6.0f) + (1.0f/detC_1_5)*(2.0f*c11-2.0f*c22)*(c11-4.0f*c12+c22) - c22*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*1.5f)*(8.0f/33.0f)
        + termB*(1.0f/sqrt_detC)*(c11-c12+c22)*((c11*(1.0f/6.0f)-c12*(2.0f/3.0f)+c22*(1.0f/6.0f))/detC+(c11*0.5f-c22*0.5f)/detC-c22*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/12.0f)-c22*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.25f)*(8.0f/33.0f)
        - c22*termB*termA*(1.0f/detC_1_5)*(c11-c12+c22)*(4.0f/33.0f);

    float dPhi11 = beta*(pow(termA, 2.0f)*((c11*(-1.0f/6.0f)+c12*(2.0f/3.0f)-c22*(1.0f/6.0f))/detC+(c11*0.5f-c22*0.5f)/detC+c11*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/12.0f)+c11*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.25f)*(41.0f/33.0f)
        + termB*(-pow(c11-c22, 2.0f)/detC_1_5 + (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/3.0f) - c11*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/6.0f) + (1.0f/detC_1_5)*(2.0f*c11-2.0f*c22)*(c11-4.0f*c12+c22) + c11*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*1.5f)*(1.0f/528.0f)
        + termA*(1.0f/detC_sq)*pow(c11-c12+c22, 3.0f)*(4.0f/81.0f)
        - (1.0f/detC_sq)*pow(c11-c12+c22, 4.0f)*((c11*(-1.0f/6.0f)+c12*(2.0f/3.0f)-c22*(1.0f/6.0f))/detC+(c11*0.5f-c22*0.5f)/detC+c11*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/12.0f)+c11*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.25f)*(1.0f/81.0f)
        - termB*termA*(1.0f/sqrt_detC)*(7.0f/198.0f)
        - c11*termA*(1.0f/pow(detC, 3.0f))*pow(c11-c12+c22, 4.0f)*(2.0f/81.0f)
        - termA*(1.0f/sqrt_detC)*(c11-c12+c22)*(-pow(c11-c22, 2.0f)/detC_1_5 + (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/3.0f) - c11*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/6.0f) + (1.0f/detC_1_5)*(2.0f*c11-2.0f*c22)*(c11-4.0f*c12+c22) + c11*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*1.5f)*(7.0f/198.0f)
        + termB*(1.0f/sqrt_detC)*(c11-c12+c22)*((c11*(-1.0f/6.0f)+c12*(2.0f/3.0f)-c22*(1.0f/6.0f))/detC+(c11*0.5f-c22*0.5f)/detC+c11*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/12.0f)+c11*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.25f)*(7.0f/198.0f)
        + c11*termB*termA*(1.0f/detC_1_5)*(c11-c12+c22)*(7.0f/396.0f))
        + K*(c11/burgers - c11/detC)
        - pow(termA, 2.0f)*((c11*(-1.0f/6.0f)+c12*(2.0f/3.0f)-c22*(1.0f/6.0f))/detC+(c11*0.5f-c22*0.5f)/detC+c11*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/12.0f)+c11*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.25f)*(12.0f/11.0f)
        + termB*(-pow(c11-c22, 2.0f)/detC_1_5 + (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/3.0f) - c11*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/6.0f) + (1.0f/detC_1_5)*(2.0f*c11-2.0f*c22)*(c11-4.0f*c12+c22) + c11*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*1.5f)*(17.0f/264.0f)
        - (1.0f/detC_1_5)*pow(c11-c12+c22, 3.0f)*(-pow(c11-c22, 2.0f)/detC_1_5 + (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/3.0f) - c11*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/6.0f) + (1.0f/detC_1_5)*(2.0f*c11-2.0f*c22)*(c11-4.0f*c12+c22) + c11*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*1.5f)*(1.0f/27.0f)
        - termB*(1.0f/detC_1_5)*pow(c11-c12+c22, 2.0f)*(1.0f/9.0f)
        + termB*termA*(1.0f/sqrt_detC)*(8.0f/33.0f)
        + c11*termB*(1.0f/detC_2_5)*pow(c11-c12+c22, 3.0f)*(1.0f/18.0f)
        + termA*(1.0f/sqrt_detC)*(c11-c12+c22)*(-pow(c11-c22, 2.0f)/detC_1_5 + (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/3.0f) - c11*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/6.0f) + (1.0f/detC_1_5)*(2.0f*c11-2.0f*c22)*(c11-4.0f*c12+c22) + c11*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*1.5f)*(8.0f/33.0f)
        - termB*(1.0f/sqrt_detC)*(c11-c12+c22)*((c11*(-1.0f/6.0f)+c12*(2.0f/3.0f)-c22*(1.0f/6.0f))/detC+(c11*0.5f-c22*0.5f)/detC+c11*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/12.0f)+c11*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.25f)*(8.0f/33.0f)
        - c11*termB*termA*(1.0f/detC_1_5)*(c11-c12+c22)*(4.0f/33.0f);

    float dPhi01 = pow(termA, 2.0f)*((c11*(-2.0f/3.0f)+c12*(8.0f/3.0f)-c22*(2.0f/3.0f))/detC + c12*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/6.0f) + c12*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.5f)*(12.0f/11.0f)
        - K*(2.0f*c12/burgers - 2.0f*c12/detC)
        + termB*(pow(c11-c22, 2.0f)/detC_1_5*4.0f - (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(4.0f/3.0f) + c12*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/3.0f) - c12*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*3.0f)*(17.0f/264.0f)
        - beta*(pow(termA, 2.0f)*((c11*(-2.0f/3.0f)+c12*(8.0f/3.0f)-c22*(2.0f/3.0f))/detC + c12*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/6.0f) + c12*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.5f)*(41.0f/33.0f)
            - termB*(pow(c11-c22, 2.0f)/detC_1_5*4.0f - (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(4.0f/3.0f) + c12*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/3.0f) - c12*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*3.0f)*(1.0f/528.0f)
            + termA*(1.0f/detC_sq)*pow(c11-c12+c22, 3.0f)*(4.0f/81.0f)
            - (1.0f/detC_sq)*((c11*(-2.0f/3.0f)+c12*(8.0f/3.0f)-c22*(2.0f/3.0f))/detC + c12*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/6.0f) + c12*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.5f)*pow(c11-c12+c22, 4.0f)*(1.0f/81.0f)
            - termB*termA*(1.0f/sqrt_detC)*(7.0f/198.0f)
            - c12*termA*(1.0f/pow(detC, 3.0f))*pow(c11-c12+c22, 4.0f)*(4.0f/81.0f)
            + termA*(1.0f/sqrt_detC)*(c11-c12+c22)*(pow(c11-c22, 2.0f)/detC_1_5*4.0f - (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(4.0f/3.0f) + c12*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/3.0f) - c12*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*3.0f)*(7.0f/198.0f)
            + termB*(1.0f/sqrt_detC)*((c11*(-2.0f/3.0f)+c12*(8.0f/3.0f)-c22*(2.0f/3.0f))/detC + c12*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/6.0f) + c12*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.5f)*(c11-c12+c22)*(7.0f/198.0f)
            + c12*termB*termA*(1.0f/detC_1_5)*(c11-c12+c22)*(7.0f/198.0f))
        - (1.0f/detC_1_5)*pow(c11-c12+c22, 3.0f)*(pow(c11-c22, 2.0f)/detC_1_5*4.0f - (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(4.0f/3.0f) + c12*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/3.0f) - c12*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*3.0f)*(1.0f/27.0f)
        + termB*(1.0f/detC_1_5)*pow(c11-c12+c22, 2.0f)*(1.0f/9.0f)
        - termB*termA*(1.0f/sqrt_detC)*(8.0f/33.0f)
        - c12*termB*(1.0f/detC_2_5)*pow(c11-c12+c22, 3.0f)*(1.0f/9.0f)
        + termA*(1.0f/sqrt_detC)*(c11-c12+c22)*(pow(c11-c22, 2.0f)/detC_1_5*4.0f - (1.0f/detC_1_5)*pow(c11-4.0f*c12+c22, 2.0f)*(4.0f/3.0f) + c12*(1.0f/detC_2_5)*pow(c11-4.0f*c12+c22, 3.0f)*(1.0f/3.0f) - c12*pow(c11-c22, 2.0f)*(1.0f/detC_2_5)*(c11-4.0f*c12+c22)*3.0f)*(8.0f/33.0f)
        + termB*(1.0f/sqrt_detC)*((c11*(-2.0f/3.0f)+c12*(8.0f/3.0f)-c22*(2.0f/3.0f))/detC + c12*(1.0f/detC_sq)*pow(c11-4.0f*c12+c22, 2.0f)*(1.0f/6.0f) + c12*pow(c11-c22, 2.0f)*(1.0f/detC_sq)*0.5f)*(c11-c12+c22)*(8.0f/33.0f)
        + c12*termB*termA*(1.0f/detC_1_5)*(c11-c12+c22)*(8.0f/33.0f);

    dPhi01 *= 0.5f;

    return float3(dPhi00, dPhi11, dPhi01);
}

kernel void compute_energy_and_forces(
    constant GPUElement* elements       [[buffer(0)]],
    constant float*      x_coords       [[buffer(1)]],
    device atomic_float* grad_forces    [[buffer(2)]],
    device atomic_float* total_energy   [[buffer(3)]],
    constant GPUParams&  params         [[buffer(4)]],
    uint eid                            [[thread_position_in_grid]],
    uint tid_in_tg                      [[thread_position_in_threadgroup]],
    uint simd_lane                      [[thread_index_in_simdgroup]],
    uint simd_group                     [[simdgroup_index_in_threadgroup]])
{
    threadgroup float simd_totals[8];
    if (tid_in_tg < 8) {
        simd_totals[tid_in_tg] = 0.0f;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float elem_e_weighted = 0.0f;

    if (eid < (uint)params.n_elements) {
        GPUElement elem = elements[eid];
        int n_vars = params.n_vars;

        // 1. Gather coordinates of the 3 vertices
        float2 pos[3];
        for (int v = 0; v < 3; v++) {
            int s_idx = elem.solver_idx[v];
            if (s_idx < 0) {
                pos[v] = elem.fixed_pos_or_trans[v];
            } else {
                pos[v] = float2(x_coords[s_idx], x_coords[n_vars + s_idx]) + elem.fixed_pos_or_trans[v];
            }
        }

        // 2. Deformation gradient F = pos * dNdX (2x3 * 3x2)
        float F00 = pos[0].x * elem.dNdX[0].x + pos[1].x * elem.dNdX[1].x + pos[2].x * elem.dNdX[2].x;
        float F01 = pos[0].x * elem.dNdX[0].y + pos[1].x * elem.dNdX[1].y + pos[2].x * elem.dNdX[2].y;
        float F10 = pos[0].y * elem.dNdX[0].x + pos[1].y * elem.dNdX[1].x + pos[2].y * elem.dNdX[2].x;
        float F11 = pos[0].y * elem.dNdX[0].y + pos[1].y * elem.dNdX[1].y + pos[2].y * elem.dNdX[2].y;

        // 3. Metric tensor C = F^T * F
        float c11 = F00 * F00 + F10 * F10;
        float c12 = F00 * F01 + F10 * F11;
        float c22 = F01 * F01 + F11 * F11;

        // 4. Lagrange reduction (with fast path)
        float m00 = 1.0f, m01 = 0.0f;
        float m10 = 0.0f, m11 = 1.0f;

        if (!(c12 >= 0.0f && c22 >= c11 && 2.0f * c12 <= c11)) {
            int iter = 0;
            while ((c12 < 0.0f || c22 < c11 || 2.0f * c12 > c11) && iter < 16) {
                if (c12 < 0.0f) {
                    c12 = -c12;
                    m01 = -m01;
                    m11 = -m11;
                }
                if (c22 < c11) {
                    float tmp_c = c11; c11 = c22; c22 = tmp_c;
                    float tmp_m0 = m00; m00 = m01; m01 = tmp_m0;
                    float tmp_m1 = m10; m10 = m11; m11 = tmp_m1;
                }
                if (2.0f * c12 > c11) {
                    float d12 = c12 - c11;
                    float d22 = c22 + c11 - 2.0f * c12;
                    c12 = d12;
                    c22 = d22;
                    m01 -= m00;
                    m11 -= m10;
                }
                iter++;
            }
        }

        // 5. Energy and derivative on reduced tensor
        float phi = calc_phi(c11, c22, c12, params.beta, params.K, params.burgers);
        float elem_e = (phi - params.zero_energy) / params.normalisation;
        elem_e_weighted = elem_e * elem.area;

        float3 dphi = calc_dphi(c11, c22, c12, params.beta, params.K, params.burgers);
        float dE00 = dphi.x / params.normalisation;
        float dE11 = dphi.y / params.normalisation;
        float dE01 = dphi.z / params.normalisation;

        // 6. Stress P = 2.0 * F * M * dE * M^T
        float MdE00 = m00 * dE00 + m01 * dE01;
        float MdE01 = m00 * dE01 + m01 * dE11;
        float MdE10 = m10 * dE00 + m11 * dE01;
        float MdE11 = m10 * dE01 + m11 * dE11;

        float A00 = MdE00 * m00 + MdE01 * m01;
        float A01 = MdE00 * m10 + MdE01 * m11;
        float A10 = MdE10 * m00 + MdE11 * m01;
        float A11 = MdE10 * m10 + MdE11 * m11;

        float P00 = 2.0f * (F00 * A00 + F01 * A10);
        float P01 = 2.0f * (F00 * A01 + F01 * A11);
        float P10 = 2.0f * (F10 * A00 + F11 * A10);
        float P11 = 2.0f * (F10 * A01 + F11 * A11);

        // 7. Compute nodal forces and atomically scatter
        float area = elem.area;
        for (int v = 0; v < 3; v++) {
            int s_idx = elem.solver_idx[v];
            if (s_idx >= 0) {
                float fx = (P00 * elem.dNdX[v].x + P01 * elem.dNdX[v].y) * area;
                float fy = (P10 * elem.dNdX[v].x + P11 * elem.dNdX[v].y) * area;

                atomic_fetch_add_explicit(&grad_forces[s_idx], fx, memory_order_relaxed);
                atomic_fetch_add_explicit(&grad_forces[n_vars + s_idx], fy, memory_order_relaxed);
            }
        }
    }

    // SIMD + Threadgroup Energy Reduction
    float sum_simd = simd_sum(elem_e_weighted);
    if (simd_lane == 0 && simd_group < 8) {
        simd_totals[simd_group] = sum_simd;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (tid_in_tg == 0) {
        float tg_total = 0.0f;
        for (uint s = 0; s < 8; ++s) {
            tg_total += simd_totals[s];
        }
        atomic_fetch_add_explicit(total_energy, tg_total, memory_order_relaxed);
    }
}
)";

class MetalContext {
public:
    id<MTLDevice> device = nil;
    id<MTLCommandQueue> queue = nil;
    id<MTLComputePipelineState> pso = nil;

    id<MTLBuffer> elem_buf = nil;
    id<MTLBuffer> x_buf = nil;
    id<MTLBuffer> grad_buf = nil;
    id<MTLBuffer> energy_buf = nil;
    id<MTLBuffer> params_buf = nil;

    size_t cached_element_count = 0;
    int cached_n_vars = 0;
    const void* cached_elements_ptr = nullptr;
    Eigen::Matrix2d cached_F_ext = Eigen::Matrix2d::Zero();
    bool initialized = false;

    static MetalContext& instance() {
        static MetalContext ctx;
        return ctx;
    }

    bool init() {
        if (initialized) return (device != nil && pso != nil);
        initialized = true;

        @autoreleasepool {
            device = MTLCreateSystemDefaultDevice();
            if (!device) {
                std::cerr << "[MetalGPU] No Metal default device found." << std::endl;
                return false;
            }

            NSError* error = nil;
            NSString* src = [NSString stringWithUTF8String:metal_shader_source];
            MTLCompileOptions* opts = [[MTLCompileOptions alloc] init];
            [opts setLanguageVersion:MTLLanguageVersion3_0];
            id<MTLLibrary> lib = [device newLibraryWithSource:src options:opts error:&error];
            if (!lib) {
                std::cerr << "[MetalGPU] Shader compilation failed: "
                          << [[error localizedDescription] UTF8String] << std::endl;
                return false;
            }

            id<MTLFunction> func = [lib newFunctionWithName:@"compute_energy_and_forces"];
            pso = [device newComputePipelineStateWithFunction:func error:&error];
            if (!pso) {
                std::cerr << "[MetalGPU] Pipeline state creation failed: "
                          << [[error localizedDescription] UTF8String] << std::endl;
                return false;
            }

            queue = [device newCommandQueue];
            params_buf = [device newBufferWithLength:sizeof(GPUParams) options:MTLResourceStorageModeShared];
            energy_buf = [device newBufferWithLength:sizeof(float) options:MTLResourceStorageModeShared];
        }
        return true;
    }

    void invalidate_cache() {
        cached_element_count = 0;
        cached_n_vars = 0;
        cached_elements_ptr = nullptr;
        cached_F_ext = Eigen::Matrix2d::Zero();
        elem_buf = nil;
        x_buf = nil;
        grad_buf = nil;
    }
};

bool is_available() {
    return MetalContext::instance().init();
}

void invalidate_mesh_cache() {
    MetalContext::instance().invalidate_cache();
}

void compute_energy_and_gradient(
    const alglib::real_1d_array &x,
    double &func,
    alglib::real_1d_array &grad,
    UserData *userData)
{
    // If not Strain_Energy_LatticeCalculator or no Metal, fall back to CPU reference
    if (dynamic_cast<Strain_Energy_LatticeCalculator*>(&userData->calculator) == nullptr) {
        minimize_energy_with_triangles(x, func, grad, userData);
        return;
    }

    MetalContext& ctx = MetalContext::instance();
    if (!ctx.init()) {
        minimize_energy_with_triangles(x, func, grad, userData);
        return;
    }

    @autoreleasepool {
        const size_t active_count = userData->active_elements.size();
        const int n_vars = static_cast<int>(userData->interior_mapping.size());
        const void* elements_ptr = userData->elements.data();
        const Eigen::Matrix2d& F_ext = userData->F_external;


        // Check if mesh buffers need repacking / reallocating
        bool need_repack = (ctx.cached_element_count != active_count ||
                            ctx.cached_n_vars != n_vars ||
                            ctx.cached_elements_ptr != elements_ptr ||
                            !ctx.cached_F_ext.isApprox(F_ext, 1e-12) ||
                            ctx.elem_buf == nil);

        if (need_repack) {
            size_t elem_bytes = active_count * sizeof(GPUElement);
            ctx.elem_buf = [ctx.device newBufferWithLength:elem_bytes options:MTLResourceStorageModeShared];

            GPUElement* gpu_elems = static_cast<GPUElement*>([ctx.elem_buf contents]);
            const auto& full_mapping = userData->full_mapping;
            const auto& points = userData->points;

            for (size_t i = 0; i < active_count; ++i) {
                size_t elem_idx = userData->active_elements[i];
                const ElementTriangle2D& elem = userData->elements[elem_idx];
                const auto& dNdX = elem.getDNdX();
                double area = elem.getReferenceArea();

                GPUElement gpu_el;
                gpu_el.area = static_cast<float>(area);
                gpu_el.padding = 0.0f;

                for (int v = 0; v < 3; ++v) {
                    int global_node_idx = elem.getNodeIndex(v);
                    const auto& [original_idx, solver_idx] = full_mapping[global_node_idx];
                    gpu_el.solver_idx[v] = solver_idx;

                    Eigen::Vector2d eff_trans = elem.getEffectiveTranslation(v);
                    if (solver_idx == -1) {
                        Eigen::Vector2d fixed_pos = points[original_idx].coord + eff_trans;
                        gpu_el.fixed_pos_or_trans[v] = simd_make_float2(static_cast<float>(fixed_pos.x()),
                                                                        static_cast<float>(fixed_pos.y()));
                    } else {
                        gpu_el.fixed_pos_or_trans[v] = simd_make_float2(static_cast<float>(eff_trans.x()),
                                                                        static_cast<float>(eff_trans.y()));
                    }

                    gpu_el.dNdX[v] = simd_make_float2(static_cast<float>(dNdX(v, 0)),
                                                      static_cast<float>(dNdX(v, 1)));
                }
                gpu_elems[i] = gpu_el;
            }

            ctx.cached_element_count = active_count;
            ctx.cached_n_vars = n_vars;
            ctx.cached_elements_ptr = elements_ptr;
            ctx.cached_F_ext = F_ext;
        }


        // Reallocate coordinates and gradient buffers if n_vars changed
        size_t grad_bytes = 2 * n_vars * sizeof(float);
        if (ctx.grad_buf == nil || ctx.grad_buf.length < grad_bytes) {
            ctx.x_buf = [ctx.device newBufferWithLength:grad_bytes options:MTLResourceStorageModeShared];
            ctx.grad_buf = [ctx.device newBufferWithLength:grad_bytes options:MTLResourceStorageModeShared];
        }

        // Copy x coordinates to shared buffer as floats
        float* x_ptr = static_cast<float*>([ctx.x_buf contents]);
        const int total_dofs = 2 * n_vars;
        for (int i = 0; i < total_dofs; ++i) {
            x_ptr[i] = static_cast<float>(x[i]);
        }

        // Zero output buffers
        std::memset([ctx.grad_buf contents], 0, grad_bytes);
        std::memset([ctx.energy_buf contents], 0, sizeof(float));

        // Update parameters
        GPUParams* p = static_cast<GPUParams*>([ctx.params_buf contents]);
        p->n_elements = static_cast<int>(active_count);
        p->n_vars = n_vars;
        p->beta = static_cast<float>(energy_functions::current_beta);
        p->K = static_cast<float>(energy_functions::current_K);
        p->burgers = static_cast<float>(energy_functions::current_burgers);
        p->zero_energy = static_cast<float>(userData->zero_energy);
        p->normalisation = static_cast<float>(userData->calculator.getUnitCellArea());
        p->padding = 0.0f;


        // Encode and dispatch Metal compute command
        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

        [enc setComputePipelineState:ctx.pso];
        [enc setBuffer:ctx.elem_buf offset:0 atIndex:0];
        [enc setBuffer:ctx.x_buf offset:0 atIndex:1];
        [enc setBuffer:ctx.grad_buf offset:0 atIndex:2];
        [enc setBuffer:ctx.energy_buf offset:0 atIndex:3];
        [enc setBuffer:ctx.params_buf offset:0 atIndex:4];

        MTLSize gridSize = MTLSizeMake(active_count, 1, 1);
        NSUInteger tgSize = std::min((NSUInteger)256, (NSUInteger)ctx.pso.maxTotalThreadsPerThreadgroup);
        MTLSize threadgroupSize = MTLSizeMake(tgSize, 1, 1);
        [enc dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
        [enc endEncoding];


        [cmd commit];
        while ([cmd status] < MTLCommandBufferStatusCompleted) {
            // spin wait
        }


        // Retrieve results into func and grad
        func = static_cast<double>(*static_cast<float*>([ctx.energy_buf contents]));

        const float* g_ptr = static_cast<float*>([ctx.grad_buf contents]);
        for (int i = 0; i < total_dofs; ++i) {
            grad[i] = static_cast<double>(g_ptr[i]);
        }

    }
}

} // namespace metal_gpu
