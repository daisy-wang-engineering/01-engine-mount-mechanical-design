# 01_engine_mount_mechanical_design

Lightweight Engine Mount Strut | Mechanical Design, Structural Analysis, Fabrication and Test Plan

Design and analysis of a tubular strut for a representative reusable launch vehicle engine mount, combining analytical mechanics, Python based parametric optimization, joint sizing, and a planned FEA to test correlation workflow.

Key Results
36% lower shank mass than the lightest solid round bar meeting the same analytical criteria.
44 mm OD × 2.5 mm wall 6061 T6 aluminum tube selected through a two variable design sweep.
low slenderness ratio accounted for by using the Johnson parabolic formula rather than Euler equations. 
+0.30 yield margin and +0.40 ultimate margin at factored design loads, with a calculated buckling margin of +0.12 against a 10% requirement.
70 to 83 kN predicted failure load across the modeled load eccentricity range of 0.5 to 2 mm.
Failure-driven design approach: end fittings sized to keep joint failure from preceding the predicted column failure, enabling a meaningful future test of analytical and FEA predictions.

Engineering Approach
Derived a representative strut load from engine thrust, strut count, angle, and geometric scaling.
Used Python to sweep tube diameter and wall thickness against strength, column buckling, and wall slenderness constraints.
Applied the Johnson column curve and secant formula to account for inelastic buckling behavior and eccentric loading.
Sized pin and lug features using published structural design methods.
Defined a planned FEA, instrumented compression test, prediction correlation, and post test redesign workflow.

Project Status
Analysis and test planning complete. FEA, fabrication, and physical testing will start once author starts Mechanical Engineering Masters Program at USC in Jan 2027.

Files
report/ Full technical report, assumptions, calculations, and test plan.
code/ Python sizing, analysis, and plotting scripts.
figures/ Load path, geometry, sizing map, and predicted response.
