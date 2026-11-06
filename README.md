# 01_engine_mount_mechanical_design

**Lightweight Engine Mount Strut | Mechanical Design, Structural Analysis, Fabrication and Test Plan**

Design and analysis of a tubular strut for a representative reusable launch vehicle engine mount, combining analytical mechanics, Python based parametric optimization, joint sizing, and a planned FEA to test correlation workflow.

**Key Results**

- 36% lower shank mass than the lightest solid round bar meeting the same analytical criteria.
- 44 mm OD × 2.5 mm wall 6061 T6 aluminum tube selected through a two variable design sweep.
- Low slenderness ratio accounted for by using the Johnson parabolic formula rather than Euler equations.
- +0.30 yield margin and +0.40 ultimate margin at factored design loads, with a calculated buckling margin of +0.12 against a 10% requirement.
- 70 to 83 kN predicted failure load across the modeled load eccentricity range of 0.5 to 2 mm.
- Controlled sequential failure design: end fittings sized to keep joint failure from preceding the predicted column failure, enabling a meaningful future test of analytical and FEA predictions.
- Conducted ANSYS finite-element and buckling analyses, validated by hand calculations, to predict yield and ultimate margins of 0.3 and 0.4

**Engineering Approach**

- Derived a representative strut load from engine thrust, strut count, angle, and geometric scaling.
- Used Python to sweep tube diameter and wall thickness against strength, column buckling, and wall slenderness constraints.
- Applied the Johnson column curve and secant formula to account for inelastic buckling behavior and eccentric loading.
- Performed Finite-element analysis in ANSYS Mechanical, tested mesh convergence, cross-checked results with hand calculations
- Defined a planned instrumented compression test, prediction correlation, and post test redesign workflow.

**Project Status**

Analysis, simulation, and test planning complete. Fabrication and physical testing will start once author starts Mechanical Engineering Masters Program at USC in Jan 2027.

**Files**
- Report/ Full technical report, assumptions, calculations, and test plan.
- Code/ Python sizing, analysis, and plotting scripts.
- Figures/ Load path, geometry, sizing map, and predicted response.
- ANSYS FEA/ Results from Finite-element Methods using ANSYS.
- CAD/ Design drawing using Solidworks.
