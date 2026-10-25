"""Figures: load-path / free-body diagrams, geometry drawing, sizing map, predicted load-deflection curves."""
import numpy as np, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, Polygon
from model import *
import sweep_tube as sw
dsg = json.load(open("design.json")); D, t, Dh, Dp, tl, e_h = dsg["D"], dsg["t"], dsg["Dh"], dsg["Dp"], dsg["tl"], dsg["e_h"]
plt.rcParams.update({"font.size": 9})
def arr(ax, p, q, c="red", lw=2, **k): ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=14, color=c, lw=lw, **k))

# ---- Fig 1: free-body diagrams ----
fig, ax = plt.subplots(1, 3, figsize=(13, 4.6))
a = ax[0]; a.set_xlim(0, 10); a.set_ylim(0, 8); a.axis("off"); a.set_title("(a) Vehicle-level load path")
a.add_patch(Rectangle((1, 6.2), 8, 0.6, fc="#dddddd", ec="k")); a.text(5, 6.9, "thrust frame / tank bulkhead", ha="center")
a.add_patch(Polygon([[4.3, 2.0], [5.7, 2.0], [6.1, 0.6], [3.9, 0.6]], fc="#bbbbbb", ec="k")); a.text(5, 0.2, "engine", ha="center")
for x0 in (1.8, 3.4, 6.6, 8.2):
    a.plot([x0, 4.6 if x0 < 5 else 5.4], [6.2, 2.0], "k-", lw=3)
arr(a, (5, 1.2), (5, 5.2), "tab:blue", 3); a.text(5.2, 4.0, "T", color="tab:blue", fontsize=12)
a.text(0.2, 1.2, "N = 4 struts\nat angle alpha = 20 deg\nP = T / (N cos alpha)", fontsize=9)
b = ax[1]; b.set_xlim(0, 10); b.set_ylim(0, 8); b.axis("off"); b.set_title("(b) Strut: two-force member, eccentric load")
b.plot([5, 5], [1, 7], "k-", lw=7)
xs = np.linspace(1, 7, 50); b.plot(5 + 0.5*np.sin(np.pi*(xs-1)/6), xs, "b:", lw=2)
arr(b, (5.8, 7.6), (5.8, 7.05), "red"); arr(b, (5.8, 0.4), (5.8, 0.95), "red")
b.text(6.2, 7.5, "P", color="red", fontsize=12); b.text(6.2, 0.2, "P", color="red", fontsize=12)
b.annotate("", xy=(5, 7.3), xytext=(5.8, 7.3), arrowprops=dict(arrowstyle="<->")); b.text(5.3, 7.75, "e", fontsize=11)
b.annotate("", xy=(5, 4), xytext=(5.5, 4), arrowprops=dict(arrowstyle="<->", color="b")); b.text(5.6, 4, "delta", color="b")
b.annotate("", xy=(3.2, 1), xytext=(3.2, 7), arrowprops=dict(arrowstyle="<->")); b.text(0.2, 4, "L = 400 mm\n(pin to pin)", fontsize=9)
b.text(8.0, 5.5, "sigma_max = P/A [1 + (e c/r^2)\n sec(KL/2r sqrt(P/AE))]", fontsize=8, ha="center")
c = ax[2]; c.set_xlim(0, 12); c.set_ylim(0, 8); c.axis("off"); c.set_title("(c) Pin joint: double shear")
c.add_patch(Rectangle((1, 3.2), 4.5, 1.6, fc="#cfe2f3", ec="k")); c.text(2.8, 4.0, "strut tongue\n(thickness t_l)", ha="center", va="center")
for y0 in (4.9, 2.3): c.add_patch(Rectangle((4, y0), 5, 0.9, fc="#eeeeee", ec="k"))
c.text(7.8, 5.35, "fork plate (t1)", ha="center"); c.text(7.8, 2.75, "fork plate (t1)", ha="center")
c.add_patch(Rectangle((4.6, 2.3), 0.6, 3.5, fc="#999999", ec="k")); c.text(4.9, 1.8, "pin", ha="center")
arr(c, (0.2, 4.0), (1.0, 4.0), "red"); c.text(0.1, 4.4, "P", color="red", fontsize=12)
arr(c, (9.8, 5.35), (9.1, 5.35), "red"); arr(c, (9.8, 2.75), (9.1, 2.75), "red"); c.text(9.9, 5.6, "P/2", color="red"); c.text(9.9, 3.0, "P/2", color="red")
c.text(6.0, 0.7, "pin bending arm = t1/2 + t_l/4 + g", ha="center")
plt.tight_layout(); plt.savefig("fig1_free_body_diagrams.png", dpi=170); plt.close()

# ---- Fig 2: geometry ----
L_tot = L + 2*e_h; Lt = e_h + 20.0; Lk = 30.0
fig, ax = plt.subplots(figsize=(13, 3.6)); ax.set_aspect("equal"); ax.axis("off")
xs_ = [0, Lt, Lt+Lk, L_tot-Lt-Lk, L_tot-Lt, L_tot]
ax.add_patch(Polygon([[xs_[0], -tl/2], [xs_[1], -tl/2], [xs_[1], -Dh/2], [xs_[2], -D/2], [xs_[3], -D/2], [xs_[4], -Dh/2], [xs_[4], -tl/2], [xs_[5], -tl/2], [xs_[5], tl/2], [xs_[4], tl/2], [xs_[4], Dh/2], [xs_[3], D/2], [xs_[2], D/2], [xs_[1], Dh/2], [xs_[1], tl/2], [xs_[0], tl/2]], fc="#cfe2f3", ec="k"))
ax.add_patch(Rectangle((xs_[2]+0.0, -(D/2-t)), xs_[3]-xs_[2], D-2*t, fc="white", ec="k", ls="--"))
for xc in (e_h, L_tot - e_h):
    ax.add_patch(Circle((xc, 0), 2.0, fc="k"))
    ax.plot([xc, xc], [-Dh/2-4, Dh/2+4], "k-.", lw=0.6)
ax.annotate("", xy=(e_h, -Dh/2-8), xytext=(L_tot-e_h, -Dh/2-8), arrowprops=dict(arrowstyle="<->")); ax.text(L_tot/2, -Dh/2-14, f"pin centre distance L = {L:.0f} mm", ha="center")
ax.annotate("", xy=(0, Dh/2+8), xytext=(L_tot, Dh/2+8), arrowprops=dict(arrowstyle="<->")); ax.text(L_tot/2, Dh/2+12, f"overall length {L_tot:.0f} mm", ha="center")
ax.text(L_tot/2, 0, f"tube OD {D:.0f} x wall {t:.1f} (bored)", ha="center", va="center", fontsize=9)
ax.text(Lt/2+5, Dh/2+0, f"tongue {tl:.0f} mm thick\npin hole {Dp:.0f} H7, e = {e_h:.0f} mm", ha="center", va="top", fontsize=8)
ax.set_xlim(-20, L_tot+20); ax.set_ylim(-Dh/2-22, Dh/2+22); ax.set_title("Strut geometry, edge view of the tongues (schematic; shoulder radii not drawn)")
plt.tight_layout(); plt.savefig("fig2_geometry.png", dpi=170); plt.close()

# ---- Fig 3: sizing map and constraint utilisation ----
ok, A_, u = sw.evaluate(sw.Dg, sw.tg, 1.0, D_over_t_max=20)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
ax[0].contourf(sw.Dg, sw.tg, ok, levels=[-.5, .5, 1.5], colors=["#f4cccc", "#d9ead3"])
cs = ax[0].contour(sw.Dg, sw.tg, rho*A_*L, levels=[0.25, 0.3, 0.35, 0.4, 0.5, 0.6], colors="gray", linewidths=.7); ax[0].clabel(cs, fmt="%.2f kg")
ax[0].contour(sw.Dg, sw.tg, u["buck"], levels=[1], colors="b", linestyles="--"); ax[0].contour(sw.Dg, sw.tg, np.maximum(u["yield_"], u["ult"]), levels=[1], colors="r")
ax[0].plot(D, t, "ko"); ax[0].plot(44, 2.25, "k^")
ax[0].set_xlabel("tube outside diameter (mm)"); ax[0].set_ylabel("wall thickness (mm)"); ax[0].set_ylim(1.5, 8); ax[0].set_xlim(28, 60)
ax[0].set_title("Feasible region (green, e = 1 mm, D/t <= 20)\nblue dashed: buckling limit, red: strength limit")
Dv = np.arange(30, 60.1, 0.5); tv = np.full_like(Dv, 2.5)
okv, Av, uv = sw.evaluate(Dv, tv, 1.0, D_over_t_max=100)
ax[1].plot(Dv, uv["yield_"], label="yield (1.25 x limit)"); ax[1].plot(Dv, uv["ult"], label="ultimate (1.4 x limit)"); ax[1].plot(Dv, uv["buck"], label="buckling (MS >= 10%)")
ax[1].axhline(1, color="k", lw=.8); ax[1].set_ylim(0, 2); ax[1].set_xlabel("tube OD at 2.5 mm wall (mm)"); ax[1].set_ylabel("demand / allowable"); ax[1].legend(); ax[1].set_title("Active criterion along the 2.5 mm wall line")
plt.tight_layout(); plt.savefig("fig3_sizing_map.png", dpi=170); plt.close()

# ---- Fig 4: predicted load vs lateral offset ----
A, I, r, c = tube(D, t); fig, ax = plt.subplots(figsize=(6.4, 4.4))
Pc_t, _ = P_crit(A, r, Fty_t)
for e in (0.5, 1.0, 2.0):
    Ps = np.linspace(0, 0.995*Pc_t, 800); ax.plot(delta_secant(Ps, A, r, e) + e, Ps/1e3, label=f"e = {e} mm")
ax.axhline(P_L/1e3, color="g", ls=":", label="limit load"); ax.axhline(P_u/1e3, color="orange", ls="--", label="1.4 x limit")
ax.axhline(Pc_t/1e3, color="k", lw=.7, label="Johnson ceiling (typical)")
ax.set_xlim(0, 8); ax.set_xlabel("mid-length lateral offset (mm)"); ax.set_ylabel("axial load (kN)"); ax.set_title("Predicted load vs lateral offset (secant model)"); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig("fig4_prediction.png", dpi=170); plt.close()
print("figures written")
