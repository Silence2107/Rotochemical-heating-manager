
import ROOT

eos_pref = "data/rootfiles/cooling_wwo_rotoch_heat/dd2_rdf_"

"""filenames = ["data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_1.4Ms_reduct_aoccdk_heavy_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_1.4Ms_reduct_aoccdk_light_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_1.4Ms_reduct_aoccdk_rotoch_heavy_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_1.4Ms_reduct_aoccdk_rotoch_light_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_2Ms_reduct_pureccdk_heavy_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_2Ms_reduct_pureccdk_light_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_2Ms_reduct_pureccdk_rotoch_heavy_atm.root",
             "data/rootfiles/cooling_wwo_rotoch_heat/apr4_eos_2Ms_reduct_pureccdk_rotoch_light_atm.root"]
"""

filenames = [eos_pref + "1.4Ms_reduct_aoccdk_heavy_atm.root",
                eos_pref + "1.4Ms_reduct_aoccdk_light_atm.root",
                eos_pref + "1.4Ms_reduct_aoccdk_rotoch_heavy_atm.root",
                eos_pref + "1.4Ms_reduct_aoccdk_rotoch_light_atm.root",
                eos_pref + "2Ms_reduct_aoccdk_heavy_atm.root",
                eos_pref + "2Ms_reduct_aoccdk_light_atm.root",
                eos_pref + "2Ms_reduct_aoccdk_rotoch_heavy_atm.root",
             eos_pref + "2Ms_reduct_aoccdk_rotoch_light_atm.root"]

heating_line_styles = [1, 1, 9, 9, 1, 1, 9, 9]
mass_line_colors = [ROOT.kRed+1, ROOT.kRed-7, ROOT.kRed+1,
                    ROOT.kRed-7, ROOT.kBlue+1, ROOT.kBlue-7, ROOT.kBlue+1, ROOT.kBlue-7]
atm_mass_line_widths = [3, 2, 3, 2, 3, 2, 3, 2]

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
    g.SetLineColor(mass_line_colors[i])
    g.SetLineStyle(heating_line_styles[i])
    g.SetLineWidth(atm_mass_line_widths[i])
    g.Draw("L")

# legend for masses
legend = ROOT.TLegend(0.15, 0.15, 0.4, 0.35)
legend.SetBorderSize(0)
legend.SetTextFont(43)
legend.SetTextSize(27)
legend.SetFillStyle(0)
legend.SetMargin(0.35)
legend.AddEntry(graphs[0], "1.4 M_{#odot}", "l")
legend.AddEntry(graphs[4], "2.0 M_{#odot}", "l")

# legend for heating switch
legend2 = ROOT.TLegend(0.15, 0.4, 0.4, 0.6)
legend2.SetBorderSize(0)
legend2.SetTextFont(43)
legend2.SetTextSize(27)
legend2.SetFillStyle(0)
legend2.SetMargin(0.35)
legend2.SetHeader("DD2+RDF HS evolution")
legend2.AddEntry(graphs[0], "w/o RH", "l")
legend2.AddEntry(graphs[2], "w/ RH", "l")

# legend for atm mass
legend3 = ROOT.TLegend(0.45, 0.15, 0.7, 0.35)
legend3.SetBorderSize(0)
legend3.SetTextFont(43)
legend3.SetTextSize(27)
legend3.SetFillStyle(0)
legend3.SetMargin(0.35)
legend3.AddEntry(graphs[0], "Heavy atm.", "l")
legend3.AddEntry(graphs[1], "Light atm.", "l")

legend.Draw()
legend2.Draw()
legend3.Draw()

# save canvas
c.SaveAs("pictures/superimpose.pdf")