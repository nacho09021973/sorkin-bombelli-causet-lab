#include <cuda_runtime.h>

#include <cmath>

extern "C" __global__ void compute_energy_kernel(
    const unsigned char* z,
    const double* x,
    const double* r,
    double rave,
    double* out,
    int n,
    int dim
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n * n;
    if (idx >= total) {
        return;
    }

    int i = idx / n;
    int j = idx % n;
    if (i >= j) {
        out[idx] = 0.0;
        return;
    }

    double rij = r[i] - r[j];
    double xij_sq = 0.0;
    for (int k = 0; k < dim; ++k) {
        double diff = x[i * dim + k] - x[j * dim + k];
        xij_sq += diff * diff;
    }

    double s2 = -(rij * rij) + xij_sq;
    double xij = sqrt(fmax(xij_sq, 0.0));
    unsigned char zij = z[i * n + j];
    double value = 0.0;
    const double roottwo = 1.4142135623730950488;

    if (zij) {
        if (s2 > 0.0) {
            value = (xij + rij) / (roottwo * rave);
        } else if (rij > 0.0) {
            value = sqrt(s2 + 2.0 * (rij * rij)) / rave;
        } else {
            value = 0.0;
        }
    } else {
        if (s2 > 0.0) {
            value = 0.0;
        } else {
            value = (fabs(rij) - xij) / (roottwo * rave);
        }
    }

    out[idx] = value;
}
