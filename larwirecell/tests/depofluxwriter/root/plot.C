// opinionated plotting

#include <TCanvas.h>
#include <TFile.h>
#include <TLine.h>
#include <TStyle.h>

struct MultiPdf {
  TCanvas canvas;
  std::string filename;

  MultiPdf(const std::string& filename) : canvas("c", "canvas", 500, 500), filename(filename)
  {
    canvas.Print((filename + "[").c_str(), "pdf");
  }
  ~MultiPdf() { close(); }
  void operator()()
  {
    canvas.Print(filename.c_str(), "pdf");
    canvas.Clear();
  }
  void close()
  {
    if (filename.empty()) { return; }
    canvas.Print((filename + "]").c_str(), "pdf");
    filename = "";
  }
};

struct Plotter {
  TFile* fp;
  MultiPdf pdf;
  Plotter(const std::string& histfile = "single_gen_dunefd_hist.root",
          const std::string& pdffile = "single_gen_dunefd_hist.pdf")
    : fp(TFile::Open(histfile.c_str(), "readonly")), pdf(pdffile)
  {}

  TH2D* get(const std::string& name)
  {
    auto h = (TH2D*)fp->Get(name.c_str());
    if (!h) { throw std::runtime_error("no hist: " + name); }
    return (TH2D*)h;
  }
};

void plot(const std::string& histfile = "single_gen_dunefd_hist.root",
          const std::string& pdffile = "single_gen_dunefd_hist.pdf")
{
  Plotter p(histfile, pdffile);

  p.pdf.canvas.SetLeftMargin(0.15);
  p.pdf.canvas.SetRightMargin(0.15);

  const double min_tick = 2350, max_tick = 2550;

  std::vector<std::string> sig_names = {"gauss", "wiener"};
  std::vector<std::vector<int>> special_channels = {{2602, 3186}, {3919, 3898}, {4286, 4288}};
  gStyle->SetOptStat(0);
  const int apa = 1;
  const char* uvw = "UVW";
  for (size_t iplane = 0; iplane < 3; ++iplane) {
    auto sc = p.get(Form("h_simpleSC_%c_%d", uvw[iplane], apa));
    for (const auto& sig_name : sig_names) {
      auto sig = p.get(Form("h_%s%d_%c", sig_name.c_str(), apa, uvw[iplane]));
      sig->SetTitle(Form("%s signals + IDE, %c-plane", sig_name.c_str(), uvw[iplane]));
      sig->SetXTitle("ticks");
      sig->SetYTitle("channels");
      sig->Draw("colz");
      sig->SetStats(false);
      sig->GetXaxis()->SetRangeUser(min_tick, max_tick);
      sc->Draw("colz same");
      TLine line;
      for (int ch : special_channels[iplane]) {
        line.DrawLine(min_tick, ch, max_tick, ch);
      }
      p.pdf();

      for (int ch : special_channels[iplane]) {
        TLegend leg(0.6, 0.8, 0.8, 0.9);

        {
          int sc_bin = sc->GetYaxis()->FindBin(ch);
          TH1D* sc_row = sc->ProjectionX("", sc_bin, sc_bin);
          sc_row->SetName(Form("simchannel_%d", ch));
          sc_row->SetTitle(
            Form("%s signals + IDE, %c-plane channel %d", sig_name.c_str(), uvw[iplane], ch));
          sc_row->GetXaxis()->SetRangeUser(min_tick, max_tick);
          sc_row->SetLineColor(kBlack);
          sc_row->SetXTitle("ticks");
          sc_row->SetYTitle("electrons");
          sc_row->Draw("hist");
          leg.AddEntry(sc_row, "IDE");
        }

        {
          int sig_bin = sig->GetYaxis()->FindBin(ch);
          TH1D* sig_row = sig->ProjectionX("", sig_bin, sig_bin);
          sig_row->SetName(Form("sigchannel_%d", ch));
          sig_row->Draw("hist same");
          sig_row->SetLineColor(kRed);
          leg.AddEntry(sig_row, "SigProc");
        }
        leg.Draw();

        p.pdf();
      }
    } // signals
    { // adc
      auto adc = p.get(Form("h_orig%d_%c", apa, uvw[iplane]));
      // median-subtract
      adc->Draw("colz");
      adc->SetStats(false);
      adc->GetXaxis()->SetRangeUser(min_tick, max_tick);
      TLine line;
      for (int ch : special_channels[iplane]) {
        line.DrawLine(min_tick, ch, max_tick, ch);
      }
      sc->Draw("same box");
      p.pdf();

      for (int ch : special_channels[iplane]) {
        TLegend leg(0.6, 0.8, 0.8, 0.9);
        {
          int sc_bin = sc->GetYaxis()->FindBin(ch);
          TH1D* sc_row = sc->ProjectionX("", sc_bin, sc_bin);
          sc_row->SetName(Form("simchannel_%d", ch));
          sc_row->SetTitle(Form("ADC + IDE, %c-plane channel %d", uvw[iplane], ch));
          sc_row->GetXaxis()->SetRangeUser(min_tick, max_tick);
          sc_row->SetXTitle("tick");
          sc_row->SetYTitle("electrons/tick");
          sc_row->SetLineColor(kBlack);
          sc_row->Draw("hist");
          leg.AddEntry(sc_row, "IDE");
          p.pdf.canvas.Update();
        }

        {
          int adc_bin = adc->GetYaxis()->FindBin(ch);
          TH1D* adc_row = adc->ProjectionX("", adc_bin, adc_bin);
          adc_row->SetName(Form("adcchannel_%d", ch));

          //scale to the pad coordinates
          float rightmax = 1.1 * adc_row->GetMaximum();
          float scale = gPad->GetUymax() / rightmax;
          adc_row->SetLineColor(kRed);
          adc_row->SetYTitle("ADC");
          adc_row->Scale(scale);
          adc_row->Draw("hist same");
          leg.AddEntry(adc_row, "ADC");

          //draw an axis on the right side
          TGaxis* axis = new TGaxis(gPad->GetUxmax(),
                                    gPad->GetUymin(),
                                    gPad->GetUxmax(),
                                    gPad->GetUymax(),
                                    0,
                                    rightmax,
                                    510,
                                    "+L");
          axis->SetLineColor(kRed);
          axis->SetLabelColor(kRed);
          axis->SetTitleColor(kRed);
          axis->SetTitle("ADC");
          axis->CenterLabels(true);
          axis->Draw();
        }
        leg.Draw();
        p.pdf();
      }
    }
  }
  p.pdf.close();
}

/*
h_daq1_U ticks=6000 chans=800
h_daq1_V ticks=6000 chans=800
h_daq1_W ticks=6000 chans=960

h_gauss1_U ticks=6000 chans=800
h_gauss1_V ticks=6000 chans=800
h_gauss1_W ticks=6000 chans=960

h_wiener1_U ticks=6000 chans=800
h_wiener1_V ticks=6000 chans=800
h_wiener1_W ticks=6000 chans=960

h_simpleSC_U_1 ticks=6000 chans=800
h_simpleSC_V_1 ticks=6000 chans=800
h_simpleSC_W_1 ticks=6000 chans=960
*/
