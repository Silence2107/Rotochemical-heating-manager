
import ROOT
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='Process m_r_diagram data into M-R preview.')
parser.add_argument('--file', type=str, help='rootfile to process')
parser.add_argument('--output', type=str, help='pdf output file')
args = parser.parse_args()

# read data
f = ROOT.TFile(args.file)
graph = f.Get("m_r_diagram")
c = ROOT.TCanvas("c", "c")
#ROOT.gPad.SetLogy()
#ROOT.gPad.SetLogx()
ROOT.gPad.SetTicks()
ROOT.gPad.SetTopMargin(0.05)
ROOT.gPad.SetLeftMargin(0.11)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetBottomMargin(0.1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

graph.Draw("AC")
graph.GetXaxis().SetTitle("R [km]")
graph.GetYaxis().SetTitle("M [Ms]")
graph.GetYaxis().SetLabelFont(43)
graph.GetYaxis().SetLabelSize(22)
graph.GetYaxis().SetTitleFont(43)
graph.GetYaxis().SetTitleSize(26)
graph.GetYaxis().SetTitleOffset(1.2)
graph.GetXaxis().SetLabelFont(43)
graph.GetXaxis().SetLabelSize(22)
graph.GetXaxis().SetTitleFont(43)
graph.GetXaxis().SetTitleSize(26)
graph.GetXaxis().SetTitleOffset(0.8)
graph.GetXaxis().SetLimits(0, 20)
graph.GetYaxis().SetRangeUser(0, 3)

c.SaveAs(args.output)