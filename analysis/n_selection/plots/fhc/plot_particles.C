#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"

void applyDuneStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    gStyle->SetTitleAlign(22);
    gStyle->SetTitleX(.5);
    gStyle->SetTitleY(.95);
    // gStyle->SetTitleAlign(33);   // Right-aligned, top
    // gStyle->SetTitleX(0.9);      // Move title to right horizontally
    // gStyle->SetTitleY(0.95);     // Optional: vertical position
    gStyle->SetTitleBorderSize(0);
    gStyle->SetFillColor(10);
    gStyle->SetCanvasColor(10);
    gStyle->SetPadColor(10);
    gStyle->SetFrameFillColor(10);
    gStyle->SetStatColor(10);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetTitleSize(.055, "xyz");
    gStyle->SetTitleSize(0.043, ""); 
    gStyle->SetTitleOffset(0.92, "xy");
    gStyle->SetTitleOffset(0.7, "z");
    gStyle->SetLabelSize(.04, "xyz");
    gStyle->SetLabelOffset(.005, "xyz");
    gStyle->SetHistLineWidth(2);
    gStyle->SetFuncColor(kRed);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetMarkerStyle(20);
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetLegendFont(42);
    gROOT->ForceStyle();
}

void drawWIPLabel() {
    TLatex* label = new TLatex(0.03, 0.08, "#font[62]{DUNE-NDLAr 2x2} Work In Progress - Simulation");
    label->SetNDC();
    label->SetTextFont(42);
    label->SetTextSize(0.03);
    label->Draw("SAME");
}
void drawVersionLabel() {
    TLatex* label = new TLatex(0.03, 0.05, "MiniRun 6.3 FHC 10E19 POT ");
    label->SetNDC();
    label->SetTextFont(42);
    label->SetTextSize(0.03);
    label->Draw("SAME");
}

void plot_particles(const char* filename) {
    applyDuneStyle();  // Apply DUNE styling

    TFile* f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    TTree* tree = (TTree*)f->Get("protons");
    if (!tree) {
        std::cerr << "Error: TTree 'protons' not found in file!" << std::endl;
        return;
    }

    const char* cut = "coincidence == 1 && ninduced == 1";
    Long64_t nPassingProtons = tree->Draw("rlen", cut, "goff");
    std::cout << "Number of protons passing the cut (coincidence == 1 && ninduced == 1): "
            << nPassingProtons << std::endl;

    // === Canvas 1 ===
    TCanvas* c1 = new TCanvas("c1", "Reco and True Proton Lengths", 800, 600);
    c1->SetGrid();
    c1->SetLogy();

    tree->Draw("rlen>>hReco(50, 0, 50)", cut, "HIST");
    tree->Draw("tlen>>hTrue(50, 0, 50)", cut, "HIST SAME");
    tree->Draw("rlen>>hFullReco(50, 0, 50)", "", "HIST SAME");

    TH1* hReco = (TH1*)gDirectory->Get("hReco");
    TH1* hTrue = (TH1*)gDirectory->Get("hTrue");
    TH1* hFullReco = (TH1*)gDirectory->Get("hFullReco");

    hReco->SetLineColor(kBlue);
    hTrue->SetLineColor(kRed);
    hFullReco->SetLineColor(kGreen);
    hReco->SetTitle("Neutron Induced Protons;Length [cm];Counts");

    TLegend* leg = new TLegend(0.70, 0.68, 0.90, 0.83);
    leg->AddEntry(hReco, "Reco Length (rlen)", "l");
    leg->AddEntry(hTrue, "True Length (tlen)", "l");
    leg->AddEntry(hFullReco, "Full Reco Length (all)", "l");
    leg->Draw();

    drawWIPLabel();  // Add watermark
    drawVersionLabel();  // Add version label
    c1->SaveAs("proton_lengths_comparison.png");

    // === Canvas 2 ===
    TCanvas* c2 = new TCanvas("c2", "Reco Length vs Distance", 800, 600);
    c2->SetGrid();
    tree->Draw("rdist:rlen>>h2D(50, 0, 50, 50, 0, 150)", cut, "COLZ");
    TH2* h2D = (TH2*)gDirectory->Get("h2D");
    if (h2D) h2D->SetTitle("Neutron Induced Protons;Reco Length [cm];Reco Distance [cm]");
    drawWIPLabel();
    drawVersionLabel();
    c2->SaveAs("proton_rlen_rdist_2D.png");

    // === Canvas 3 ===
    TCanvas* c3 = new TCanvas("c3", "Reco Length vs True Length", 800, 600);
    c3->SetGrid();
    tree->Draw("tlen:rlen>>h2D_rlen_tlen(50, 0, 50, 50, 0, 50)", cut, "COLZ");
    TH2* h2D_rlen_tlen = (TH2*)gDirectory->Get("h2D_rlen_tlen");
    if (h2D_rlen_tlen) h2D_rlen_tlen->SetTitle("Neutron Induced Protons;Reco Length [cm];True Length [cm]");
    drawWIPLabel();
    drawVersionLabel();
    c3->SaveAs("proton_rlen_tlen_2D.png");

    // === Canvas 4 ===
    TCanvas* c4 = new TCanvas("c4", "Reco Distance vs True Distance", 800, 600);
    c4->SetGrid();
    tree->Draw("tdist:rdist>>h2D_rdis_tdis(50, 0, 50, 50, 0, 50)", cut, "COLZ");
    TH2* h2D_rdis_tdis = (TH2*)gDirectory->Get("h2D_rdis_tdis");
    if (h2D_rdis_tdis) h2D_rdis_tdis->SetTitle("Neutron Induced Protons;Reco Distance [cm];True Distance [cm]");
    drawWIPLabel(); 
    drawVersionLabel();
    c4->SaveAs("proton_rdist_tdis_2D.png");

    f->Close();
}
