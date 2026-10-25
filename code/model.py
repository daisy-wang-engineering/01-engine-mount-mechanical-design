"""
Common model: 6061-T6 tubular strut with integral clevis-tongue ends, loaded in axial compression through pins.
Units: N, mm, MPa, kg.  Imported by the sweep, end-fitting and prediction scripts.
"""
import numpy as np

# ---------- material 6061-T6 ----------
E, G, rho = 68.9e3, 26.0e3, 2.70e-6            # MPa, MPa, kg/mm^3
Fty_d, Ftu_d = 240.0, 290.0                    # design (minimum) yield / ultimate, MPa
Fty_t, Ftu_t = 276.0, 310.0                    # typical values (used for test prediction only)
Fsu, Fbru, Fbry = 207.0, 607.0, 386.0          # shear; bearing ultimate / yield (e/D = 2.0)
# ---------- pin steel (4130 values from Air Force Stress Manual worked example) ----------
ksi = 6.89476
Ftu_P, Fty_P, Fsu_P = 125*ksi, 103*ksi, 82*ksi
kbP = 1.56                                      # plastic bending coefficient, ductile pin
# ---------- factors and criteria ----------
FS_y, FS_u, MS_buck_min = 1.25, 1.40, 0.10
K = 1.0                                         # pinned-pinned effective length factor
L = 400.0                                       # pin-centre to pin-centre length, mm
e_ecc = 1.0                                     # load eccentricity used in analysis, mm
# ---------- load definition ----------
T_full, D_eng = 845e3, 920.0                    # Merlin 1D sea-level thrust (N) and diameter (mm)
N_str, alpha = 4, np.radians(20.0)
lam = L/D_eng
P_full = T_full/(N_str*np.cos(alpha))
P_L = lam**2*P_full                             # limit load on the test strut
P_y, P_u = FS_y*P_L, FS_u*P_L

# ---------- round tube section ----------
def tube(D, t):
    d = D - 2*t
    A = np.pi/4*(D**2 - d**2); I = np.pi/64*(D**4 - d**4); r = np.sqrt(I/A)
    return A, I, r, D/2

# ---------- column models ----------
def Cc(Fty): return np.pi*np.sqrt(2*E/Fty)
def P_crit(A, r, Fty=Fty_d):
    S = K*L/r
    sig = np.where(S >= Cc(Fty), np.pi**2*E/S**2, Fty - Fty**2/(4*np.pi**2*E)*S**2)
    return sig*A, S
def phi(P, A, r): return K*L/(2*r)*np.sqrt(P/(A*E))
def sigma_secant(P, A, r, c, e=e_ecc):
    f = phi(P, A, r)
    return np.where(f < np.pi/2, P/A*(1 + e*c/r**2/np.cos(np.minimum(f, 1.55))), np.inf)
def delta_secant(P, A, r, e=e_ecc):
    f = phi(P, A, r)
    return np.where(f < np.pi/2, e*(1/np.cos(np.minimum(f, 1.55)) - 1), np.inf)
def sigma_local(D, t):                          # classical axial shell buckling stress, 0.605 E t / R (R = mean radius)
    return 0.605*E*t/((D - t)/2)
def first_yield(A, r, c, e, Fty=Fty_t):
    Pc, _ = P_crit(A, r, Fty); Ps = np.linspace(1, 0.9999*Pc, 60000)
    sg = sigma_secant(Ps, A, r, c, e)
    return Ps[np.argmax(sg >= Fty)] if np.any(sg >= Fty) else Ps[-1]     # no yield before the Johnson load -> ceiling

# ---------- lug and pin (Air Force Stress Manual, double shear, axial) ----------
def pin_capacity(Dp, t_fork, t_tongue, g):
    arm = t_fork/2 + t_tongue/4 + g                        # eq. 9-13
    return dict(arm=arm,
                ult_bend = 0.1963*kbP*Dp**3*Ftu_P/arm,       # eq. 9-15
                elastic  = 0.1963*Dp**3*Fty_P/arm,
                shear    = 1.571*Dp**2*Fsu_P)                # eq. 9-12
def tongue_capacity(Dp, e, w, t, Kn=0.74, K_b=1.33, Ftu=Ftu_d, Fty=Fty_d):
    Fb = min(K_b*Ftu, 1.304*K_b*Fty, Fbru, 1.304*Fbry) if e/Dp >= 1.5 else np.nan
    return dict(bearing=Fb*Dp*t, net=Kn*Ftu*(w - Dp)*t)
