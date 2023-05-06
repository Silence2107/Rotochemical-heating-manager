
import ROOT

filenames = ["data/rootfiles/cooling/apr4_eos_2Ms_no_reduct.root",
             "data/rootfiles/cooling/apr4_eos_2Ms_reduct.root",
             "data/rootfiles/cooling/hybrid_probe_eos_2Ms_no_reduct.root",
             "data/rootfiles/cooling/hybrid_probe_eos_2Ms_reduct.root",
             "data/rootfiles/cooling/ist_eos_2Ms_no_reduct.root",
             "data/rootfiles/cooling/ist_eos_2Ms_reduct.root"]
            
eos_line_styles = [1, 1, 5, 5, 9, 9]
reduction_line_colors = [4, 2, 4, 2, 4, 2]

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
h.GetYaxis().SetRangeUser(7e2, 7e6)
h.GetXaxis().SetLimits(1e-12, 1e7)

# draw graphs
for i, g in enumerate(graphs):
    g.SetLineColor(reduction_line_colors[i])
    g.SetLineStyle(eos_line_styles[i])
    g.SetLineWidth(2)
    g.Draw("L")

# legend for EoS
legend = ROOT.TLegend(0.15, 0.1, 0.4, 0.3)
legend.SetBorderSize(0)
legend.SetTextFont(43)
legend.SetTextSize(27)
legend.SetFillStyle(0)
legend.SetMargin(0.35)
legend.AddEntry(graphs[0], "EOS: APR4", "l")
legend.AddEntry(graphs[2], "EOS: Hybrid with thin crust", "l")
legend.AddEntry(graphs[4], "EOS: IST + NZD-NV", "l")

# legend for reduction
legend2 = ROOT.TLegend(0.15, 0.4, 0.4, 0.6)
legend2.SetBorderSize(0)
legend2.SetTextFont(43)
legend2.SetTextSize(27)
legend2.SetFillStyle(0)
legend2.SetMargin(0.35)
legend2.AddEntry(graphs[0], "No reduction", "l")
legend2.AddEntry(graphs[1], "Max. reduction", "l")

legend.Draw()
legend2.Draw()

# save canvas
c.SaveAs("pictures/superimpose.pdf")