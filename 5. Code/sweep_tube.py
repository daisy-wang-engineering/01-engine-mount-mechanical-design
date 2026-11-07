"""Two-variable sizing sweep (outside diameter D, wall t) for the round tubular strut shank."""
import numpy as np, json
from model import *

Ds = np.arange(28.0, 60.01, 0.5); ts = np.arange(1.5, 12.01, 0.25)
Dg, tg = np.meshgrid(Ds, ts, indexing="ij")

def evaluate(D, t, e=e_ecc, D_over_t_max=25.0, t_min=2.0):
    A, I, r, c = tube(D, t)
    u_y = sigma_secant(P_y, A, r, c, e)/Fty_d                  # yield at 1.25 x limit
    u_u = sigma_secant(P_u, A, r, c, e)/Ftu_d                  # ultimate at 1.4 x limit
    Pc, S = P_crit(A, r)
    u_b = (1 + MS_buck_min)*P_u/Pc                             # buckling margin >= 10 %
    u_l = 3*Fty_d/sigma_local(D, t)                            # classical local buckling stress >= 3 x yield
    ok = (u_y <= 1) & (u_u <= 1) & (u_b <= 1) & (u_l <= 1) & (D/t <= D_over_t_max) & (t >= t_min)
    return ok, A, dict(yield_=u_y, ult=u_u, buck=u_b, local=u_l)

def optimum(e=e_ecc, D_max=60.0, **kw):
    ok, A, u = evaluate(Dg, tg, e, **kw); ok &= (Dg <= D_max)
    m = np.where(ok, A, np.inf); i = np.unravel_index(np.argmin(m), m.shape)
    ut = {k: float(v[i]) for k, v in u.items()}
    return Dg[i], tg[i], A[i], rho*A[i]*L, ut, max(ut, key=ut.get)

if __name__ == "__main__":
    print(f"lambda={lam:.4f}  P_full={P_full/1e3:.1f} kN  P_L={P_L/1e3:.2f}  P_y={P_y/1e3:.2f}  P_u={P_u/1e3:.2f} kN")
    print("\n e(mm) Dmax(mm) | OD(mm) wall(mm)  A(mm2) mass(kg) | util yield ult buck local | governing   [D/t<=25, t>=2.0]")
    rows = {}
    for e in (0.5, 1.0, 2.0, 5.0):
        for Dm in (40, 50, 60):
            D, t, A, m, ut, gov = optimum(e, Dm)
            print(f" {e:4.1f}   {Dm:4d}    | {D:5.1f}  {t:5.2f}   {A:6.1f}  {m:6.3f}  | {ut['yield_']:.2f} {ut['ult']:.2f} {ut['buck']:.2f} {ut['local']:.2f} | {gov}")
            rows[f"e{e}_D{Dm}"] = dict(D=float(D), t=float(t), A=float(A), mass=float(m), gov=gov)
    json.dump(rows, open("sweep_tube_results.json", "w"), indent=1)
    print("\nSensitivity to the wall-slenderness bound (e = 1 mm, OD <= 50 mm):")
    for dt in (15, 20, 25, 30):
        D, t, A, m, ut, gov = optimum(1.0, 50, D_over_t_max=dt)
        print(f" D/t<={dt}: OD={D:.1f} t={t:.2f} A={A:.0f} mm2 mass={m:.3f} kg  util(y,u,b)=({ut['yield_']:.2f},{ut['ult']:.2f},{ut['buck']:.2f})  governing: {gov}")
    for Dsol in np.arange(15, 40, 0.5):                       # lightest solid round bar, for comparison
        A = np.pi/4*Dsol**2; r = Dsol/4; c = Dsol/2; Pc, _ = P_crit(A, r)
        if sigma_secant(P_y, A, r, c)/Fty_d <= 1 and sigma_secant(P_u, A, r, c)/Ftu_d <= 1 and (1 + MS_buck_min)*P_u/Pc <= 1:
            print(f"lightest solid round bar: D={Dsol} mm, A={A:.0f} mm2, mass={rho*A*L:.3f} kg"); break
