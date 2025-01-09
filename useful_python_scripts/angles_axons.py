import numpy as np
from scipy.special import hyp0f1
from scipy.integrate import quad
from scipy.special import erfi
import pandas as pd
from scipy.interpolate import interp1d

def c2toKappa(c2, tol=1e-3, kappa_bounds=(0, 64)):
    if c2 == 1.0:
        kappa = np.inf
        return kappa
    elif c2 < 1/3:
        kappa = 0
        return kappa
    kappas = np.arange(kappa_bounds[0], kappa_bounds[1], tol)
    Fs = np.sqrt(np.pi) / 2 * np.exp(-kappas) * erfi(np.sqrt(kappas))
    c2s = 1 / (2 * np.sqrt(kappas) * Fs) - 1 / (2 * kappas)
    c2s[-1] = 1
    idx_c2 = np.nanargmin(np.abs(c2 - c2s))
    kappa = kappas[idx_c2]
    if kappa < 0:
        kappa = 0
    return kappa


def KappaToC2(kappa):
    # from kappa, calculate c2, the mean cos ^ 2 of the angle between axons and
    # main bundle orientation(an easier metric of orientations dispersion, c2
    # varies between 1 / 3(isotropic) and 1(perfectly parallel axons)
    Fs = np.sqrt(np.pi) / 2 * np.exp(-kappa) * erfi(np.sqrt(kappa))
    c2 = 1 / (2 * np.sqrt(kappa) * Fs) - 1 / (2 * kappa)
    return c2

# Sample from a density using rejection sampling
def sample_from_density(p, kappa, x_range, num_samples, M=None):

    xmin, xmax = x_range
    if M is None:
        x = np.linspace(xmin, xmax, 1000)
        M = np.max(p(x, kappa))

    samples = []
    while len(samples) < num_samples:
        x_candidate = np.random.uniform(xmin, xmax)
        y = np.random.uniform(0, M)
        if y <= p(x_candidate, kappa):
            samples.append(x_candidate)

    return np.array(samples)

# Define the normalized PDF
def p(x, kappa):
    normalization_constant = 1 / hyp0f1(0.5, 1.5, np.array([kappa]))
    return normalization_constant * np.exp(kappa * (np.cos(x)**2))

# Compute the normalized PDF
def normalized_pdf(x, pdf, norm_const):
    return pdf(x) / norm_const

# Verify the normalization
def normalization_constant(kappa):
    integral, _ = quad(p, -np.pi/2, np.pi/2, args=(kappa,))
    return integral

# Draw angles to match a specific c2
def draw_angles(c2, num_samples=1):
    if not (0 <= c2 <= 1):
        raise ValueError("c2 must be between 0 and 1.")
    kappa = c2toKappa(c2)
    if kappa < 0:
        kappa = 0
    angles = sample_from_density(p, kappa, [-np.pi / 2, np.pi / 2], num_samples)
    return angles

def angle_from_cdf(kappa, df, num_samples=1000):
    """
    Computes angles corresponding to random values for the Watson distribution with `kappa`.
    
    Args:
        kappa (float): The Watson distribution's concentration parameter.
        df (pd.DataFrame): The CDF table with angles as columns (degrees) and kappas as index.
        num_samples (int): Number of angles to sample.

    Returns:
        list: List of sampled angles in degrees.
    """
    max_kappa = df.index.max()
    if (kappa > max_kappa) :
        return np.zeros(num_samples)
    min_kappa = df.index.min()
    if (kappa < min_kappa) :
        return np.ones(num_samples)*np.nan
    
    # Interpolate for the kappa row
    cdf_interpolator = interp1d(
        df.index, df.values, axis=0, kind="quadratic", fill_value="extrapolate"
    )
    interpolated_row = cdf_interpolator(kappa)

    # Convert column names to numeric angles
    angles = df.columns.astype(float)

    # Interpolate the angle for random x values
    angle_interpolator = interp1d(
        interpolated_row,
        angles,
        kind="quadratic",
        bounds_error=False,
        fill_value=(angles[0], angles[-1]),
    )

    # Generate random values and map them to angles
    random_values = np.random.uniform(0, 1, num_samples)
    sampled_angles = angle_interpolator(random_values)

    return sampled_angles.tolist()

# Main testing
if __name__ == "__main__":


    # Sample DataFrame
    data = {
        5: [0.025, 0.055, 0.110, 0.210, 0.380, 0.620],
        10: [0.095, 0.200, 0.370, 0.610, 0.850, 0.980],
        15: [0.20, 0.39, 0.65, 0.88, 0.99, 1.00],
        30: [0.56, 0.84, 0.98, 1.00, 1.00, 1.00],
        45: [0.80, 0.97, 1.00, 1.00, 1.00, 1.00],
        60: [0.91, 0.99, 1.00, 1.00, 1.00, 1.00],
        75: [0.97, 1.00, 1.00, 1.00, 1.00, 1.00]
    }
    df = pd.DataFrame(data, index=[4, 8, 16, 32, 64, 128])

    print(df)

    # Example usage
    c2 = 0.85
    kappa = c2toKappa(c2)
    print(kappa)
    angles = draw_angles(c2, num_samples=1)
    cos_2 = []
    for angle in angles:
        cos_2.append(np.cos(np.radians(angle))**2)
        
    print(np.mean(cos_2))
