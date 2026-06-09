#!/usr/bin/env python3
"""Build a PDF progress note on the GP regression of rotational rate coefficients."""
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, PageBreak)

styles = getSampleStyleSheet()
H1 = styles["Heading1"]; H2 = styles["Heading2"]; N = styles["BodyText"]
N.spaceAfter = 6
small = ParagraphStyle("small", parent=N, fontSize=8, leading=10)
mono = ParagraphStyle("mono", parent=N, fontName="Courier", fontSize=8, leading=10)

story = []


def P(t, s=N): story.append(Paragraph(t, s))
def SP(h=8): story.append(Spacer(1, h))


def tbl(data, widths=None):
    t = Table(data, colWidths=widths, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#33455e")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8.5),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#eef2f7")]),
        ("ALIGN", (1, 0), (-1, -1), "CENTER"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
    ]))
    story.append(t)


# ---------------------------------------------------------------- title
P("Predicting rotational rate coefficients k(T) from model potentials", H1)
P("Gaussian-process regression &mdash; progress note", styles["Heading3"])
P("Atom + linear-rigid-rotor (MOLSCAT) dataset. Out-of-sample validation, "
  "feature-set ablation, and example k(T) predictions with uncertainty.", small)
SP(6)

# ---------------------------------------------------------------- goal/data
P("1. Goal and dataset", H2)
P("<b>Goal.</b> Predict state-to-state de-excitation rate coefficients "
  "k(i&rarr;f, T) for an atom colliding with a linear rigid rotor directly from "
  "the interaction potential, bypassing the close-coupling scattering "
  "calculation &mdash; with a calibrated uncertainty on every prediction.")
P("<b>Dataset.</b> 1000 Sobol-sampled model potentials of Born&ndash;Mayer + "
  "dispersion form, V(R,&theta;) = V0(R) + V1(R)P<sub>1</sub> + V2(R)P<sub>2</sub>, with "
  "V0 = A e<super>-aR</super> - C6 R<super>-6</super>, an odd wall anisotropy &eta; and "
  "tail anisotropy a<sub>1</sub> in V1, and an even anisotropy a<sub>2</sub> in "
  "V2 = a<sub>2</sub>V0. Rates were computed with MOLSCAT, thermally averaged, "
  "and fitted (64&thinsp;838 rows over He and H<sub>2</sub> partners).")
P("<b>Target representation.</b> Each k(T) curve is stored as "
  "log<sub>10</sub>k(T) = &Sigma;<sub>n=0..5</sub> a<sub>n</sub> (T<super>-1/6</super>)<super>n</super>. "
  "Over 5&ndash;50 K the monomials in T<super>-1/6</super> are strongly collinear "
  "(condition number ~10<super>7</super>, coefficient correlation ~1.0), so the raw "
  "a<sub>n</sub> are poor regression targets; a decorrelated Chebyshev basis on "
  "the same interval (condition ~4.5) is the conditioned alternative. Here we "
  "regress log<sub>10</sub>k at fixed temperatures, which sidesteps the issue.")

# ---------------------------------------------------------------- features
P("2. Feature sets", H2)
P("<b>Fingerprint (raw).</b> arcsinh(V(R,&theta;)/200) sampled on a fixed grid "
  "R &isin; {2,3,3.5,4,5,6,8,10} &Aring; &times; &theta; &isin; {0,90,180}&deg; "
  "(24 numbers) + log B + &mu;. No physics distilled.")
P("<b>Geometry (derived).</b> ~8 physical descriptors computed from the "
  "potential: well depths &epsilon;<sub>0,90,180</sub>; the odd/even depth "
  "contrasts d&epsilon;<sub>odd</sub>=&epsilon;<sub>0</sub>-&epsilon;<sub>180</sub>, "
  "d&epsilon;<sub>even</sub>; minimum radius R<sub>m</sub>; wall radius "
  "R<sub>wall</sub>; log B; &mu;; transition gap &Delta;E.")
P("<b>DWBA+ overlap.</b> Geometry plus distorted-wave-Born-motivated coupling "
  "overlaps: &int; V1 over the V0 well (depth- and plain-weighted), its "
  "V0-normalised form, and &eta;, a<sub>1</sub>. (V2 = a<sub>2</sub>V0 carries "
  "only the scalar a<sub>2</sub>, so it adds no independent radial feature.)")

# ---------------------------------------------------------------- results A
P("3. Feature-set ablation (5-fold CV R<super>2</super>)", H2)
P("Target: log<sub>10</sub>k(He, 1&rarr;0, 10 K), N = 891 potentials. GBT = "
  "gradient-boosted trees; GP = Gaussian process (Mat&eacute;rn, ARD).", small)
tbl([["feature set", "# feat", "GBT R2", "GP R2"],
     ["raw fingerprint", "26", "0.283", "0.559"],
     ["geometry (derived)", "10", "0.896", "0.903"],
     ["+ DWBA overlap", "15", "0.905", "0.901"]],
    widths=[2.4 * inch, 0.9 * inch, 1.0 * inch, 1.0 * inch])
SP(4)
P("The raw fingerprint reproduces the prior ~0.5 baseline; physics-derived "
  "geometry features take the same target and CV to ~0.90. The DWBA overlap "
  "adds only +0.002&ndash;0.010, smaller than the fold-to-fold scatter "
  "(&plusmn;0.04) &mdash; not a statistically real gain.", small)
SP(6)
P("Robustness across the &Delta;j=1 ladder and the H<sub>2</sub> partner "
  "(GBT R<super>2</super>):", N)
tbl([["target", "N", "fingerprint", "geometry", "+DWBA", "&Delta;DWBA"],
     ["He 1->0", "891", "0.283", "0.896", "0.905", "+0.010"],
     ["He 2->1", "899", "0.398", "0.904", "0.913", "+0.010"],
     ["He 3->2", "891", "0.435", "0.884", "0.886", "+0.002"],
     ["H2 1->0", "884", "0.384", "0.932", "0.935", "+0.003"]],
    widths=[1.2 * inch, 0.6 * inch, 1.1 * inch, 1.0 * inch, 0.8 * inch, 0.9 * inch])
SP(4)
P("<b>Conclusion.</b> The derived geometry set is sufficient and preferred: "
  "equal accuracy to DWBA+ within noise, ~8 interpretable features already in "
  "the library (no extra computation), lower overfitting risk. Simpler is "
  "better.", N)
SP(10)

# ---------------------------------------------------------------- validation
P("4. Out-of-sample validation", H2)
P("Best model (GP, geometry + DWBA), out-of-fold predictions with the GP "
  "predictive &sigma;.", small)
tbl([["metric", "5 K", "10 K"],
     ["out-of-fold R2", "0.893", "0.903"],
     ["RMSE (dex)", "0.168", "0.163"],
     ["95% within (dex)", "0.35", "0.29"],
     ["coverage +-1 sigma (ideal 0.68)", "0.69", "0.71"],
     ["coverage +-2 sigma (ideal 0.95)", "0.91", "0.91"]],
    widths=[3.0 * inch, 1.0 * inch, 1.0 * inch])
SP(4)
P("Typical error ~0.16 dex (a factor ~1.45 in k); 95% of predictions within a "
  "factor ~2. <b>The GP uncertainty is well calibrated</b> (coverage matches "
  "the ideal), so low-&sigma; predictions are trustworthy and high-&sigma; ones "
  "are flagged. Mild heavy tails (&plusmn;2&sigma; = 0.91 &lt; 0.95) are the "
  "near-threshold cases; 5 K is marginally harder than 10 K, as expected.", small)
SP(6)
P("Validation panels &mdash; predicted vs actual, residuals, calibration, and "
  "|residual| vs &sigma; (10 K):", small)
story.append(Image("validation_He10.png", width=6.3 * inch, height=5.25 * inch))

story.append(PageBreak())

# ---------------------------------------------------------------- examples
P("5. Example k(T) predictions with uncertainty", H2)
P("Ten potentials spanning large/small V0 (depth &epsilon;<sub>90</sub>), V1 "
  "(odd anisotropy d&epsilon;<sub>odd</sub>) and V2 (even d&epsilon;<sub>even</sub>), "
  "for transitions 1&rarr;0, 3&rarr;2 (&Delta;j=1) and 2&rarr;0 (&Delta;j=2), "
  "He partner. All ten were held out of training, so each curve is "
  "out-of-sample. Black = dataset fit; blue dashed = GP mean; band = "
  "&plusmn;1&sigma;.", small)
SP(4)
story.append(Image("example_kT_curves.png", width=4.25 * inch, height=8.5 * inch))

story.append(PageBreak())

# ---------------------------------------------------------------- summary
P("6. Summary and status", H2)
for b in [
    "Predicting k(T) directly from the potential works: out-of-fold "
    "R<super>2</super> ~ 0.90 across &Delta;j=1 transitions and both partners.",
    "Representation dominates: derived geometry descriptors beat the raw "
    "V(R,&theta;) fingerprint by a wide margin (0.56 &rarr; 0.90).",
    "DWBA coupling-overlap features confirm the physics but add no "
    "statistically significant accuracy; the ~8-feature geometry model is the "
    "recommended, simplest sufficient choice.",
    "The Gaussian process supplies calibrated per-prediction uncertainty "
    "(coverage ~0.70 at 1&sigma;), turning point estimates into trustworthy, "
    "self-flagging predictions.",
    "Known limits: results are interpolation within the Sobol parameter box; "
    "a heavy tail of near-threshold cases is the main residual error; "
    "validated so far at single temperatures for He/H<sub>2</sub>.",
]:
    P("&bull;&nbsp; " + b, N)
SP(6)
P("<b>Next steps.</b> full-curve validation via the decorrelated Chebyshev "
  "coefficients; out-of-distribution test (hold out parameter-space edges, "
  "confirm &sigma; inflates); minimal-feature trim; inspect near-threshold "
  "outliers.", small)
SP(8)
P("<b>Artifacts.</b> regress_kT.py (ablation), analysis.py (fingerprint vs "
  "derived, ladder), make_valplot.py (validation_He{5,10}.png), "
  "make_examples.py (example_kT_curves.png). Data: potentials_library.csv, "
  "kT_coeffs.csv; evaluator potentials.py.", small)

doc = SimpleDocTemplate("GP_progress_note.pdf", pagesize=A4,
                        topMargin=0.7 * inch, bottomMargin=0.7 * inch,
                        leftMargin=0.8 * inch, rightMargin=0.8 * inch,
                        title="GP regression of rate coefficients - progress note")
doc.build(story)
print("wrote GP_progress_note.pdf")
