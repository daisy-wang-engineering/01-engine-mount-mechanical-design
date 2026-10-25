"""Baseline selection, end-fitting sizing, mass/stiffness, and test predictions for the tubular strut."""
import numpy as np, json
from model import *
from sweep_tube import optimum, evaluate, Dg, tg      # sweep helpers

print("\n=== SENSITIVITY TO WALL-SLENDERNESS BOUND (e=1 mm, OD<=50) ===")
for dt in (15, 20, 25, 30):
    D, t, A, m, ut, gov = optimum(1.0, 50, D_over_t_max=dt)
    print(f" D/t<={dt}: OD={D:.1f} t={t:.2f} A={A:.0f} mass={m:.3f} kg util(y,u,b)=({ut['yield_']:.2f},{ut['ult']:.2f},{ut['buck']:.2f}) gov={gov}")

# ---------------- baseline shank (D/t <= 20, OD <= 50) ----------------
D, t, A_opt, m_opt, ut, gov = optimum(1.0, 50, D_over_t_max=20)
D, t = 44.0, 2.5                                    # sweep optimum OD 44.0; wall 2.25 rounded up to 2.5 mm for boring tolerance
A, I, r, c = tube(D, t); d = D - 2*t
print(f"\nBASELINE SHANK: OD={D} wall={t} ID={d}  A={A:.1f} mm2  I={I:.0f} mm4  r={r:.3f} mm  S={K*L/r:.2f}  Cc(design)={Cc(Fty_d):.2f}")
Pc_d, S = P_crit(A, r); Pc_t, _ = P_crit(A, r, Fty_t)
print(f" Johnson Pcr design={Pc_d/1e3:.1f} kN (MS vs 1.1Pu {Pc_d/(1.1*P_u)-1:+.2f}); typical ceiling={Pc_t/1e3:.1f} kN; Euler ideal={np.pi**2*E*I/(K*L)**2/1e3:.1f} kN")
for nm, P, F in (("limit", P_L, None), ("1.25 x limit", P_y, Fty_d), ("1.4 x limit", P_u, Ftu_d)):
    f = float(phi(P, A, r)); sg = float(sigma_secant(P, A, r, c)); dl = float(delta_secant(P, A, r))
    print(f" {nm:13s} P={P/1e3:6.2f} kN  P/A={P/A:6.1f}  phi={f:.4f} sec={1/np.cos(f):.4f} delta={dl:.3f} mm  sigma_max={sg:.1f} MPa" + (f"  MS={F/sg-1:+.2f}" if F else ""))
print(f" local buckling: classical sigma_cl={sigma_local(D,t):.0f} MPa = {sigma_local(D,t)/Fty_d:.1f} x design yield; D/t={D/t:.0f}")
print(f" shank mass per mm of tube = {rho*A:.5f} kg/mm")

# ---------------- end fittings ----------------
P_end, P_pin_el = 1.25*Pc_t, 1.10*Pc_t
print(f"\nEND REQUIREMENTS: strut ceiling {Pc_t/1e3:.1f} kN -> ends >= {P_end/1e3:.1f} kN, pin elastic >= {P_pin_el/1e3:.1f} kN")
t_fork, g = 10.0, 0.5
best = None
for Dh in (50.0, 56.0, 60.0, 63.5):                     # head bar stock diameter (tongue is milled from it)
    for Dp in (16, 18, 20, 22, 25, 28, 30):
        for tl in np.arange(10, 0.8*Dh+0.1, 1.0):
            w = np.sqrt(Dh**2 - tl**2); e = float(np.ceil(1.6*Dp/5)*5)
            pn = pin_capacity(Dp, t_fork, tl, g); tc = tongue_capacity(Dp, e, w, tl)
            ok = (pn['elastic'] >= P_pin_el) and (pn['ult_bend'] >= P_end) and (pn['shear'] >= P_end) and (tc['bearing'] >= P_end) and (tc['net'] >= P_end)
            if ok:
                vol = tl*w*(e+20)
                if best is None or vol < best[0]: best = (vol, Dh, Dp, tl, w, e, pn, tc)
vol, Dh, Dp, tl, w, e_h, pn, tc = best
print(f" minimum-volume tongue: head stock D={Dh}, pin D={Dp} mm, tongue thickness={tl:.0f}, chord width={w:.1f}, e={e_h:.0f} (e/D={e_h/Dp:.2f}, a={e_h-Dp/2:.1f})")
# round to a clean drawing: choose from the search result
Dp, tl, e_h = float(Dp), float(tl), float(e_h); w = np.sqrt(Dh**2 - tl**2)
pn = pin_capacity(Dp, t_fork, tl, g); tcd = tongue_capacity(Dp, e_h, w, tl); tct = tongue_capacity(Dp, e_h, w, tl, Ftu=Ftu_t, Fty=Fty_t)
print(f" pin: arm={pn['arm']:.2f} mm  ult bending={pn['ult_bend']/1e3:.0f}  elastic={pn['elastic']/1e3:.0f}  shear={pn['shear']/1e3:.0f} kN")
print(f" tongue (design min): bearing={tcd['bearing']/1e3:.1f}  net={tcd['net']/1e3:.1f} kN ; typical: bearing={tct['bearing']/1e3:.1f} net={tct['net']/1e3:.1f} kN")
for nm, cap, req in (("tongue bearing", tcd['bearing'], P_end), ("tongue net section", tcd['net'], P_end), ("pin shear", pn['shear'], P_end), ("pin ult bending", pn['ult_bend'], P_end), ("pin elastic", pn['elastic'], P_pin_el)):
    print(f"   {nm:19s} MS vs requirement = {cap/req-1:+.2f}   | vs ultimate design load {cap/P_u-1:+.2f}")
print(f" tear-out sanity (2 planes, Fsu={Fsu}): {2*Fsu*tl*(e_h-Dp/2*np.cos(np.radians(40)))/1e3:.0f} kN")

# ---------------- geometry, mass, stiffness ----------------
L_tot = L + 2*e_h; L_tongue = e_h + 20.0; L_taper = 30.0; L_tube = L_tot - 2*(L_tongue + L_taper)
R = Dh/2; zt = tl/2
A_tongue = 2*(zt*np.sqrt(R**2 - zt**2) + R**2*np.arcsin(zt/R)); A_head = np.pi*Dh**2/4; A_tube_end = np.pi*D**2/4
vol_tube = A*L_tube
vol_taper = 2*(np.pi*L_taper/12)*(Dh**2 + Dh*D + D**2)             # frustum, solid
vol_tongue = 2*A_tongue*L_tongue - 2*np.pi/4*Dp**2*tl
mass = rho*(vol_tube + vol_taper + vol_tongue)
print(f"\nGEOMETRY: total length {L_tot:.0f}, tongue zone {L_tongue:.0f}, taper {L_taper:.0f}, tube {L_tube:.0f} mm; head stock {Dh} mm; tube OD {D} wall {t}")
print(f"MASS = {mass:.3f} kg  (tube {rho*vol_tube:.3f}, tapers {rho*vol_taper:.3f}, tongues {rho*vol_tongue:.3f})")
A_net_h = tl*(w - Dp)
comp = L_tube/(E*A) + 2*L_taper/(E*(A_head+A_tube_end)/2) + 2*20/(E*A_net_h) + 2*(L_tongue-20)/(E*A_tongue)
print(f"AXIAL STIFFNESS: EA/L(tube, L=400)={E*A/L/1e3:.1f} kN/mm ; segmented-bar model={1/comp/1e3:.1f} kN/mm")

# ---------------- predictions (typical properties) ----------------
print("\nPREDICTION (typical): ceiling (Johnson) = %.1f kN" % (Pc_t/1e3))
pred = {}
for ee in (0.5, 1.0, 2.0):
    fy = first_yield(A, r, c, ee); pred[ee] = fy
    print(f"  e={ee}: first yield = {fy/1e3:.1f} kN  ({fy/P_L:.2f} x limit)   delta at limit = {float(delta_secant(P_L,A,r,ee)):.3f} mm")
print(f"  planning load for machine = 1.5 x ceiling = {1.5*Pc_t/1e3:.0f} kN = {1.5*Pc_t/4448.2:.0f} kip ; fork bearing at that load = {1.5*Pc_t/(Dp*2*t_fork):.0f} MPa")
print(f"  peak/mass at e=1: {pred[1.0]/1e3/mass:.0f} kN/kg")
json.dump(dict(D=D,t=t,Dh=Dh,Dp=Dp,tl=tl,e_h=e_h,mass=float(mass),ceiling=float(Pc_t),first_yield={str(k):float(v) for k,v in pred.items()}), open("design.json","w"), indent=1)
