
import ROOT

eos_pref = "project/Cooling/temp_1.4_"

"""filenames = ["data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_1.4Ms_reduct_aoccdk_heavy_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_1.4Ms_reduct_aoccdk_light_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_1.4Ms_reduct_aoccdk_rotoch_heavy_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_1.4Ms_reduct_aoccdk_rotoch_light_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_2Ms_reduct_pureccdk_heavy_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_2Ms_reduct_pureccdk_light_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_2Ms_reduct_pureccdk_rotoch_heavy_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_2Ms_reduct_pureccdk_rotoch_light_atm.root"]
"""

filenames = [eos_pref + "eq.root",
                eos_pref + "noneq.root",
                eos_pref + "e_rh.root",
                eos_pref + "u_rh.root",
                eos_pref + "eu_rh.root"]

line_styles = [1, 9, 2, 2, 2]
line_colors = [ROOT.kBlack, ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen]

# extract "Graph" from every single rootfile
infs = []
graphs = []
for filename in filenames:
    f = ROOT.TFile(filename)
    infs.append(f)
    g = f.Get("Graph")
    graphs.append(g)

# create a canvas
c = ROOT.TCanvas("c", "c")
ROOT.gPad.SetLogy()
ROOT.gPad.SetLogx()
ROOT.gPad.SetTicks()
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetLeftMargin(0.11)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetBottomMargin(0.1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

# dummy histogram
h = ROOT.TH1D()
h.Draw('AXIS')
h.GetXaxis().SetTitle("t[yr]")
h.GetYaxis().SetTitle("T^{#infty}_{s}[K]")
h.GetYaxis().SetLabelFont(43)
h.GetYaxis().SetLabelSize(22)
h.GetYaxis().SetTitleFont(43)
h.GetYaxis().SetTitleSize(26)
h.GetYaxis().SetTitleOffset(0.5)
h.GetXaxis().SetLabelFont(43)
h.GetXaxis().SetLabelSize(22)
h.GetXaxis().SetTitleFont(43)
h.GetXaxis().SetTitleSize(26)
h.GetXaxis().SetTitleOffset(0.9)
h.GetYaxis().SetRangeUser(7e1, 8e6)
h.GetXaxis().SetLimits(1e-12, 1e10)

# draw graphs
for i, g in enumerate(graphs):
    g.SetLineColor(line_colors[i])
    g.SetLineStyle(line_styles[i])
    g.SetLineWidth(3)
    g.Draw("L")

# legend
legend = ROOT.TLegend(0.15, 0.2, 0.4, 0.6)
legend.SetBorderSize(0)
legend.SetTextFont(43)
legend.SetTextSize(27)
legend.SetFillStyle(0)
legend.SetMargin(0.35)
legend.SetHeader("DD2+RDF HS evolution, 1.4 M_{#odot}")
legend.AddEntry(graphs[0], "eq. cool.", "l")
legend.AddEntry(graphs[1], "noneq. cool.", "l")
legend.AddEntry(graphs[2], "eq. cool. + e imbalance", "l")
legend.AddEntry(graphs[3], "eq. cool. + u imbalance", "l")
legend.AddEntry(graphs[4], "eq. cool. + eu imbalance", "l")


legend.Draw()

# save canvas
c.SaveAs("pictures/superimpose.pdf")