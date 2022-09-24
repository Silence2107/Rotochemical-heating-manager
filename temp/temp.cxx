

#include "../include/auxiliaries.h"

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>

int main()
{
    // probe cubic interpolation with 10 points
    std::vector<double> x;
    std::vector<double> y;
    for (double i = 0; i < 10; i += 1.0)
    {
        x.push_back(i);
        y.push_back(sin(i));
    }
    //double val = auxiliaries::interpolate(x, y, auxiliaries::InterpolationMode::kCubic, 3.00);
    //draw cubic interpolation
    std::vector<double> x2;
    std::vector<double> y2;
    for (double i = x.front(); i < x.back(); i += 0.01)
    {
        x2.push_back(i);
        y2.push_back(auxiliaries::interpolate(x, y, auxiliaries::InterpolationMode::kCubic, i));
    }
    // draw TMultiGraph
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TMultiGraph *mg = new TMultiGraph();
    TGraph *g1 = new TGraph(x.size(), x.data(), y.data());
    TGraph *g2 = new TGraph(x2.size(), x2.data(), y2.data());
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kRed);
    g2->SetLineColor(kBlue);
    mg->Add(g1);
    mg->Add(g2);
    mg->Draw("APL");
    c1->SaveAs("temp.pdf");
    return 0;
}