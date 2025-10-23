import numpy as np
import math

print("=" * 80)
print("TRAC BOOM DESIGN CALCULATOR")
print("=" * 80)
print()

# ==============================================================================
# USER DEFINED PARAMETERS
# ==============================================================================
# NOTE: All calculations performed in SI units (meters, Pascals)

print("INPUT PARAMETERS:")
print("-" * 80)

# Mandrel Specifications
R_mandrel = 0.025  # [m] Mandrel radius - adjust based on available stock
print(f"Mandrel radius:              {R_mandrel * 1000:.2f} mm")
print(f"Mandrel diameter:            {R_mandrel * 2000:.2f} mm")

# Cross-Section Geometry
w = 0.005  # [m] Web width (bonded region between flanges)
fh = 20 * 1e-3  # [m] Target flattened height of deployed boom
print(f"Web width:                   {w * 1000:.2f} mm")
print(f"Target flattened height:     {fh * 1000:.2f} mm")

# Laminate Properties
t = 0.2 * 1e-3  # [m] Total laminate thickness per flange
print(f"Laminate thickness:          {t * 1e6:.1f} μm ({t * 1000:.3f} mm)")

# Boom Length
L = 2.1  # [m] Total boom length
print(f"Boom length:                 {L:.2f} m")

# Material Properties (Carbon Fiber/Epoxy)
E = 128e9  # [Pa] Longitudinal elastic modulus (typical for carbon/epoxy)
G12 = 7.5e9  # [Pa] In-plane shear modulus
tau_y = 70e6  # [Pa] Matrix yield shear strength
phi = math.radians(5)  # [rad] Fiber misalignment angle (5° typical)
print(f"Elastic modulus (E):         {E / 1e9:.1f} GPa")
print(f"Shear modulus (G12):         {G12 / 1e9:.1f} GPa")
print(f"Matrix yield strength:       {tau_y / 1e6:.0f} MPa")
print(f"Fiber misalignment angle:    {math.degrees(phi):.1f}°")

# Manufacturing Constant
k_spring_in = 0.835  # Spring-in factor (validate with test coupons!)
print(f"Spring-in factor:            {k_spring_in:.3f}")

# Design Strain Limits
eps_yield = 0.01  # [unitless] Maximum allowable strain (1.0% - moderate)
print(f"Design strain limit:         {eps_yield * 100:.1f}%")

print()

# ==============================================================================
# SECTION 1: CURED GEOMETRY CALCULATIONS
# ==============================================================================
print("\n" + "=" * 80)
print("SECTION 1: CURED BOOM GEOMETRY")
print("=" * 80)

# Actual flange radius after spring-in
R_cured = R_mandrel * k_spring_in
print(f"Actual flange radius (r):    {R_cured * 1000:.2f} mm")

# Flange opening angle (in radians)
# Derived from: fh = 2*R_cured*sin(theta/2) + w
theta = 2 * math.asin((fh - w) / (2 * R_cured))
theta_deg = math.degrees(theta)
print(f"Flange opening angle (θ):    {theta_deg:.1f}° ({theta:.3f} rad)")

# Arc length of each flange
s = R_cured * theta
print(f"Flange arc length:           {s * 1000:.2f} mm")

# Total perimeter of cross-section
perimeter = 2 * s + w
print(f"Cross-section perimeter:     {perimeter * 1000:.2f} mm")

# Cross-sectional area (approximate)
A_cross = 2 * t * (s + w / 2)
print(f"Cross-sectional area:        {A_cross * 1e6:.3f} mm²")

# Fabric width needed for layup (with 50% margin on each side)
width = (s + w) * 1.5 * 2
print(f"Fabric width needed:         {width * 1000:.1f} mm (with margins)")

print()

# ==============================================================================
# SECTION 2: STRAIN ANALYSIS
# ==============================================================================
print("=" * 80)
print("SECTION 2: STRAIN ANALYSIS")
print("=" * 80)

# Transverse strain during flange flattening
eps_22 = t / (2 * R_cured)
print(f"Transverse strain (ε₂₂):     {eps_22:.6f} ({eps_22 * 100:.3f}%)")
if eps_22 <= eps_yield:
    print(f"  ✓ PASS - Below yield limit of {eps_yield * 100:.1f}%")
else:
    print(f"  ✗ FAIL - Exceeds yield limit of {eps_yield * 100:.1f}%")
    print(f"  → Reduce thickness or increase flange radius")

# Minimum hub radius from strain limit
R_hub_min = t / (2 * eps_yield)
print(f"\nMinimum hub radius (strain): {R_hub_min * 1000:.2f} mm")
print(f"Minimum hub diameter:        {R_hub_min * 2000:.2f} mm")

# Longitudinal strain at minimum hub radius
eps_11 = t / (2 * R_hub_min)
print(f"Longitudinal strain (ε₁₁):   {eps_11:.6f} ({eps_11 * 100:.3f}%)")

print()

# ==============================================================================
# SECTION 3: STRESS ANALYSIS
# ==============================================================================
print("=" * 80)
print("SECTION 3: STRESS ANALYSIS")
print("=" * 80)

# Stresses in fully coiled region (simple beam theory)
long_stress_calc = -E * t / R_hub_min  # Negative = compressive
trans_stress_calc = E * t / (2 * R_cured)

print("Fully Coiled Region (Analytical):")
print(f"  Longitudinal stress (σ_xx): {long_stress_calc / 1e6:.1f} MPa (compressive)")
print(f"  Transverse stress (σ_yy):   {trans_stress_calc / 1e6:.1f} MPa")

# Stresses in transition region (with stress concentration)
stress_concentration_factor = 2.8  # Typical from FEA studies
long_stress_max = abs(long_stress_calc) * stress_concentration_factor
trans_stress_max = abs(trans_stress_calc) * stress_concentration_factor

print(f"\nTransition Region (with {stress_concentration_factor:.1f}× concentration):")
print(f"  Peak longitudinal stress:   {long_stress_max / 1e6:.1f} MPa (compressive)")
print(f"  Peak transverse stress:     {trans_stress_max / 1e6:.1f} MPa")

print()

# ==============================================================================
# SECTION 4: FAILURE ANALYSIS
# ==============================================================================
print("=" * 80)
print("SECTION 4: MATERIAL FAILURE CRITERION")
print("=" * 80)

# Fiber microbuckling failure criterion
gamma_y = tau_y / G12  # Yield shear strain
compressive_failure_stress = G12 / (1 + phi / gamma_y)

print(f"Yield shear strain (γ_y):    {gamma_y:.6f}")
print(f"Compressive failure (σ_c):   {compressive_failure_stress / 1e6:.1f} MPa")
print()

# Check if transition region stress exceeds failure criterion
safety_factor_transition = compressive_failure_stress / long_stress_max
print(f"Safety Factor Analysis:")
print(f"  Transition stress:          {long_stress_max / 1e6:.1f} MPa")
print(f"  Material strength:          {compressive_failure_stress / 1e6:.1f} MPa")
print(f"  Safety factor:              {safety_factor_transition:.2f}")

if safety_factor_transition >= 1.2:
    print(f"  ✓ SAFE - Adequate safety margin")
elif safety_factor_transition >= 1.0:
    print(f"  ⚠ MARGINAL - Consider design improvements")
else:
    print(f"  ✗ FAILURE LIKELY - Redesign required!")
    print(f"  → Increase hub radius or reduce thickness")

print()

# ==============================================================================
# SECTION 5: MINIMUM HUB RADIUS (COMPREHENSIVE)
# ==============================================================================
print("=" * 80)
print("SECTION 5: MINIMUM HUB RADIUS DETERMINATION")
print("=" * 80)

K_stress = 2.0  # Conservative stress concentration factor for design

# Constraint 1: Strain-based minimum radius
R_strain = t / (2 * eps_yield)
print(f"Constraint 1 - Strain Limit:")
print(f"  R_min (strain):             {R_strain * 1000:.2f} mm")

# Constraint 2: Strength-based minimum radius
R_stress = K_stress * E * t / compressive_failure_stress
print(f"\nConstraint 2 - Material Strength (with K={K_stress:.1f}):")
print(f"  R_min (strength):           {R_stress * 1000:.2f} mm")

# Constraint 3: Conservative empirical formula
R_conservative = 450 * t  # Conservative rule of thumb
print(f"\nConstraint 3 - Empirical Formula (450×t):")
print(f"  R_min (empirical):          {R_conservative * 1000:.2f} mm")

# Final design minimum (take maximum of all constraints)
R_min_hub = max(R_strain, R_stress, R_conservative)
D_min_hub = 2 * R_min_hub

print(f"\n{'=' * 40}")
print(f"DESIGN MINIMUM HUB RADIUS:   {R_min_hub * 1000:.2f} mm")
print(f"DESIGN MINIMUM HUB DIAMETER: {D_min_hub * 1000:.2f} mm")
print(f"{'=' * 40}")

# Check against user's assumed hub radius
if R_hub_min < R_min_hub:
    print(f"\n⚠ WARNING: Initially calculated R_hub_min ({R_hub_min * 1000:.2f} mm)")
    print(f"           is less than comprehensive R_min ({R_min_hub * 1000:.2f} mm)")
    print(f"  → Use {R_min_hub * 1000:.2f} mm as minimum hub radius")

print()

# ==============================================================================
# SECTION 6: FLANGE FLATTENABILITY CHECK
# ==============================================================================
print("=" * 80)
print("SECTION 6: FLANGE FLATTENABILITY")
print("=" * 80)

# Minimum flange radius to allow flattening
r_min_flange = t / (2 * eps_yield)
print(f"Minimum flange radius (r_min): {r_min_flange * 1000:.2f} mm")
print(f"Actual flange radius (r):      {R_cured * 1000:.2f} mm")

if R_cured >= r_min_flange:
    margin = (R_cured / r_min_flange - 1) * 100
    print(f"  ✓ FLATTENABLE - {margin:.1f}% margin")
else:
    print(f"  ✗ NOT FLATTENABLE - Increase flange radius or reduce thickness")

print()

# ==============================================================================
# SECTION 7: PACKAGING ANALYSIS
# ==============================================================================
print("=" * 80)
print("SECTION 7: PACKAGING ANALYSIS")
print("=" * 80)

# Number of wraps around hub (approximate)
# Assumes boom flattens to height fh and coils with no gaps
n_wraps = L / (np.pi * D_min_hub) if D_min_hub > 0 else 0
print(f"Number of wraps (approx):    {n_wraps:.1f} wraps")

# Packaged outer diameter (approximate)
# Each wrap adds approximately 2×(flattened height)
D_packaged = D_min_hub + 2 * n_wraps * fh
print(f"Packaged outer diameter:     {D_packaged * 1000:.1f} mm")

# Packaging volume (cylindrical approximation)
V_packaged = np.pi * (D_packaged / 2)**2 * fh  # Volume of flat coil
print(f"Packaging volume (approx):   {V_packaged * 1e6:.1f} cm³")

# Packaging efficiency metric
V_deployed = A_cross * L  # Volume of material deployed
packing_ratio = V_deployed / V_packaged * 100 if V_packaged > 0 else 0
print(f"Packing ratio:               {packing_ratio:.2f}%")

print()

# ==============================================================================
# SECTION 8: MANDREL SPECIFICATIONS FOR MANUFACTURING
# ==============================================================================
print("=" * 80)
print("SECTION 8: MANDREL MANUFACTURING SPECIFICATIONS")
print("=" * 80)

print("For Two-Step Process (Cylindrical Mandrels):")
print(f"  • Material:       Aluminum 6061-T6 or Carbon Foam")
print(f"  • Radius:         {R_mandrel * 1000:.2f} mm")
print(f"  • Diameter:       {R_mandrel * 2000:.2f} mm")
print(f"  • Length:         {(L + 0.075) * 1000:.0f} mm (boom + 75mm handling)")
print(f"  • Surface finish: Ra < 0.8 μm")
print(f"  • Concentricity:  ± 0.05 mm")

print(f"\nFor Co-Curing Process (U-Shaped Mandrels):")
print(f"  • Inner radius:   {R_mandrel * 1000:.2f} mm")
print(f"  • Opening angle:  {theta_deg / 2:.1f}° (half of boom angle)")
print(f"  • Base width:     {w * 1000 / 2:.2f} mm (half web width)")
print(f"  • Require:        Two matching left/right mandrels")

print()

# ==============================================================================
# SECTION 9: DESIGN RECOMMENDATIONS
# ==============================================================================
print("=" * 80)
print("SECTION 9: DESIGN RECOMMENDATIONS & NEXT STEPS")
print("=" * 80)

print("\n1. FABRICATION PRIORITIES:")
print(f"   • Fabricate test coupons (300-500 mm length)")
print(f"   • Validate spring-in factor (currently {k_spring_in:.3f})")
print(f"   • Test coiling at R = {R_min_hub * 1000 * 1.3:.1f} mm (1.3× minimum)")
print(f"   • Gradually reduce radius in 5mm increments")

print("\n2. QUALITY CONTROL:")
print(f"   • Measure cured flange radius (target: {R_cured * 1000:.2f} mm ± 0.3 mm)")
print(f"   • Verify opening angle (target: {theta_deg:.1f}° ± 2°)")
print(f"   • Check web width (target: {w * 1000:.2f} mm ± 0.5 mm)")
print(f"   • Inspect for fiber misalignment (keep φ < 5°)")

print("\n3. STRESS MITIGATION:")
if safety_factor_transition < 1.5:
    print(f"   ⚠ Consider these improvements:")
    print(f"   • Reduce laminate thickness (try t = {t * 0.75 * 1e6:.0f} μm)")
    print(f"   • Variable curvature cross-section (24% stress reduction)")
    print(f"   • Use [0/90/0] instead of [0/90]s (21% stress reduction)")
    print(f"   • Increase hub radius to {R_min_hub * 1.2 * 1000:.1f} mm")
else:
    print(f"   ✓ Design appears robust with SF = {safety_factor_transition:.2f}")

print("\n4. MATERIAL SELECTION:")
print(f"   • Carbon/Epoxy prepreg: {t * 1e6:.0f} μm thick")
print(f"   • Adhesive film for web: Match CTE to prepreg")
print(f"   • Release agent: Teflon spray or Frekote")

print()
print("=" * 80)
print("END OF TRAC BOOM DESIGN CALCULATION")
print("=" * 80)
