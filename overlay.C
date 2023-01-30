
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <math.h>

using std::cout;
using std::endl;

/* ===================================

This is the regular overlay script example for beginners.
Author: Sourabh Dube
Usage : At the root prompt
    [].x overlay.C

This file is modified to work in Root6.

The primary function in this file is overlay().

   =================================== */



// Declare all the other functions needed. See the function implementations for more
// comments about what the functions do.
void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill);
void decorate(TLegend *g, float textSize, TString legendheader);
float get_nevents(TH1F *hst, float bin_lo, float bin_hi);
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi);

// This is the primary function. It opens two histogram files, extracts two histograms from
// them, and plots them on top of each other (i.e. overlays them)
void overlay()
{

  //First declare the file names (just strings)
  TString file1 = "Outputs/hst_TTZ.root";
  TString file2 = "Outputs/hst_WZ.root";
  TString file3 = "Outputs/hst_ZZ.root";

  //OVL 2
  TString file4 = "Outputs/hst_TTZ.root";
  TString file5 = "Outputs/hst_WZ.root";
  TString file6 = "Outputs/hst_ZZ.root";
 
  //OVL 3
  TString file7 = "Outputs/hst_TTZ.root"; 
  TString file8 = "Outputs/hst_WZ.root";
  TString file9 = "Outputs/hst_ZZ.root";

  
  cout<<"test1"<<endl; //Declare other constants, strings that you might need here.

  //Declare the name of the plot that you want to overlay
  //(you can open the histogram file to see other names)
  TString plotname1 = "3_lp_evt_MET_Pt";
  TString plotname2 = "3_lp_evt_MET_Pt"; //(here we are picking the same plot from other file)
  TString plotname3 = "3_lp_evt_MET_Pt";

  //OVL 2
  TString plotname4 = "LT_3lp_evt";
  TString plotname5 = "LT_3lp_evt"; //(here we are picking the same plot from other file)
  TString plotname6 = "LT_3lp_evt";

  //OVL 3
  TString plotname7 = "HT_3lp_evt";
  TString plotname8 = "HT_3lp_evt"; //(here we are picking the same plot from other file)
  TString plotname9 = "HT_3lp_evt";
  
  
  // TString plotname4 = "DiJet_mass";
  //Also give fancy name for the axis titles
  // Axis titles for OVL 1
  TString x1title = "MET"; // Can use latex-type commands
                                               // in strings. For example
                                               // "#mu^{1} p_{T}"
  TString y1title = "Events"; // Or "Events"

  //Axis titles for OVL 2
  TString X2title = "LT"; // Can use latex-type commands
                                               // in strings. For example
                                               // "#mu^{1} p_{T}"
  TString Y2title = "Events"; // Or "Events"


  //Axis titles for OVL 3
  TString x3title = "HT";
  TString y3title = "Events";


    //Now let us open the files
  TFile *file_1 = new TFile(file1);
  TFile *file_2 = new TFile(file2);
  TFile *file_3 = new TFile(file3);

  // OVL 2
  TFile *file_4 = new TFile(file4);
  TFile *file_5 = new TFile(file5);
  TFile *file_6 = new TFile(file6);

  //OVL 3
  TFile *file_7 = new TFile(file7);
  TFile *file_8 = new TFile(file8);
  TFile *file_9 = new TFile(file9);

  

  //Now open the respective histograms from the files
  TH1F *h1 = (TH1F*)file_1->Get(plotname1);
  TH1F *h2 = (TH1F*)file_2->Get(plotname2);
  TH1F *h3 = (TH1F*)file_3->Get(plotname3);

  //OVL 2
  TH1F *h4 = (TH1F*)file_4->Get(plotname4);
  TH1F *h5 = (TH1F*)file_5->Get(plotname5);
  TH1F *h6 = (TH1F*)file_6->Get(plotname6);

  //OVL 3
  TH1F *h7 = (TH1F*)file_7->Get(plotname7);
  TH1F *h8 = (TH1F*)file_8->Get(plotname8);
  TH1F *h9 = (TH1F*)file_9->Get(plotname9);
 
  // TH1F *h4 = (TH1F*)file_4->Get(plotname4);
  cout<<"test2"<<endl;



  //Decorate the histograms using function decorate 
  // See function definition below for syntax
  // See https://root.cern.ch/root/html/TColor.html for color names.
  decorate(h1,x1title,y1title,"",kBlue,2,kRed-7,29,0);
  decorate(h2,x1title,y1title,"",kRed,2,kGreen+2,34,0);
  decorate(h3,x1title,y1title,"",kGreen+2,2,kGreen+2,34,0);

  //OVL 2
  decorate(h4,X2title,Y2title,"",kBlue,2,kRed-7,29,0);
  decorate(h5,X2title,Y2title,"",kRed,2,kGreen+2,34,0);
  decorate(h6,X2title,Y2title,"",kGreen+2,2,kGreen+2,34,0);

  //OVL 3
  decorate(h7,x3title,y3title,"",kBlue,2,kRed-7,29,0);
  decorate(h8,x3title,y3title,"",kRed,2,kGreen+2,34,0);
  decorate(h9,x3title,y3title,"",kGreen+2,2,kGreen+2,34,0);
  // decorate(h4,xtitle,ytitle,"",kYellow+4,2,kGreen+2,34,0);




  //Now let us set the last bin as the overflow bin
   int nbins = h2->GetNbinsX();
   h1->SetBinContent(nbins, h1->GetBinContent(nbins+1) + h1->GetBinContent(nbins) );
   h2->SetBinContent(nbins, h2->GetBinContent(nbins+1) + h2->GetBinContent(nbins) );
  // h3->SetBinContent(nbins, h3->GetBinContent(nbins+1) + h3->GetBinContent(nbins) );
  // h1->SetBinContent(0, h1->GetBinContent(1) + h1->GetBinContent(0) );
  //h2->SetBinContent(0, h2->GetBinContent(1) + h2->GetBinContent(0) );
  //h3->SetBinContent(0, h3->GetBinContent(1) + h3->GetBinContent(0) );
 
  


  // Now rebin the histograms if needed
  // Group nrebins bins together
  int nrebins = 1;
  h1->Rebin(nrebins);
  h2->Rebin(nrebins);
  h3->Rebin(nrebins);
  



  //Now we scale the histograms.
  /* Scaling can be done either by some outside measure (such as size of samples)
     or by normalizing the histograms such that they have the same integral.
     Conventionally, when comparing shapes, we normalize to 1, so let us 
     try that here.
  */
  // h1->Scale(1.0/h1->Integral());
  // h2->Scale(1.0/h2->Integral());
  // h3->Scale(1.0/h3->Integral());
  //  h4->Scale(1.0/h4->Integral());
  



  //Now let us declare a canvas 
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  cout<<"test3"<<endl;



  // Now we define a legend
  TLegend *lg1 = new TLegend(0.55,0.50,0.85,0.75,NULL,"NDC");
  decorate(lg1,0.05,""); // Decorated the legend using function below.
  lg1->AddEntry(h1,"MET for min 3lepton events in TTZ & ZZ sample","lf"); // Added the two entries for the two histograms
  lg1->AddEntry(h2,"MET for min 3lepton events in WZ sample","lp");           // we shall be drawing.
  // lg1->AddEntry(h3,"HT ZZ sample","lp");           // we shall be drawing.
  // lg1->AddEntry(h4,"Di-Jet Mass30","lp");           // we shall be drawing.



  //Now we Add two Histograms
  h1->Add(h3);
  //We set the stat box for the first one to be invisible first
  h1->SetStats(0);



  // Now to draw the histograms and the legend
  // int nrebins = 5;
  h1->Draw("hist");            
  h2->Draw("hist same");
  // h4->Draw("hist same");
  lg1->Draw();



  //OVL 2
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  cout<<"test3"<<endl;



  // Now we define a legend
  TLegend *lg2 = new TLegend(0.55,0.50,0.85,0.75,NULL,"NDC");
  decorate(lg2,0.05,""); // Decorated the legend using function below.
  lg2->AddEntry(h4,"LT for min 3lepton events in TTZ & ZZ sample","lf"); // Added the two entries for the two histograms
  lg2->AddEntry(h5,"LT for min 3lepton events in WZ sample","lp");       // we shall be drawing.
  h4->Add(h6);


  //We set the stat box for the first one to be invisible first
  h4->SetStats(0);




  // Now to draw the histograms and the legend
  // int nrebins = 5;
  h5->Draw("hist");            
  h4->Draw("hist same");
 
  lg2->Draw();



    //Canvas for OVL 3 
  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  cout<<"test3"<<endl;


  // Now we define a legend
  TLegend *lg3 = new TLegend(0.55,0.50,0.85,0.75,NULL,"NDC");
  decorate(lg1,0.05,""); // Decorated the legend using function below.
  lg3->AddEntry(h7,"HT for min 3lepton events in TTZ & ZZ sample","lf"); // Added the two entries for the two histograms
  lg3->AddEntry(h8,"HT for min 3lepton events in WZ sample","lp");           // we shall be drawing.
  // lg1->AddEntry(h3,"HT ZZ sample","lp");           // we shall be drawing.
  // lg1->AddEntry(h4,"Di-Jet Mass30","lp");           // we shall be drawing.
  h7->Add(h9);


  //We set the stat box for the first one to be invisible first
  h7->SetStats(0);



  // Now to draw the histograms and the legend
  // int nrebins = 5;
  h8->Draw("hist");            
  h7->Draw("hist same");
  // h4->Draw("hist same");
  lg3->Draw();

  /* To draw with stat boxes for each histogram
     -- Dont use SetStats(0)
     -- Then draw them, first one with option Draw(), next ones with option Draw("sames")
     -- The s at the end is for stats box
  */

  
  //following code Prints ratio in the root 


  //for MET overlay
  cout<< "For MET Histogram"<<endl;
  float N_wz_1 = h2->Integral(); // get_nevents(h2,33,91);// No upper bound is necessary
  float N_others_1 = h1->Integral(); //get_nevents(h1,33,91);
  cout<<"N_wz = "<<N_wz_1<<"   "<< "N_others = "<<N_others_1<<"   "<< "The ratio R_1 is "<< N_wz_1 / N_others_1 << endl;


  //for LT overlay
  cout<< "For LT Histogram"<<endl;
  float N_wz_2 =h5->Integral();// get_nevents(h5,110,160);
  float N_others_2 =h4->Integral();// get_nevents(h4,110,160);
  cout<<"N_wz = "<<N_wz_2<<"   "<< "N_others = "<<N_others_2<<"   "<< "The ratio R_2 is "<< N_wz_2 / N_others_2 << endl;
 


  //for HT overlay

 cout<< "For HT Histogram"<<endl;
 float N_wz_3 =h8->Integral();// get_nevents(h8,30,100);
 float N_others_3 = h7->Integral();//et_nevents(h7,30,100);
  cout<<"N_wz = "<<N_wz_3<<"   "<< "N_others = "<<N_others_3<<"   "<< "The ratio R_3 is "<< N_wz_3 / N_others_3 << endl;
 
  
}
            
/*
This function decorates the histogram, by setting the titles of the X and Y axes,
by setting the line and marker colors, the marker style, and by deciding whether
to fill in the histogram with a color or not.
 */
void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) {

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);

  h->SetLineColor(linecolor); 
  h->SetLineWidth(linewidth);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);

  if(tofill==1) h->SetFillColor(markercolor);
  h->SetMarkerSize(1.0);
  h->SetTitle(title);
}

/*
  This function decorates the legend, by setting the textsize, and filling in the
  legend with a pure white color.
*/
void decorate(TLegend *g, float textSize, TString legendheader)
{
  g->SetTextSize(textSize);
  g->SetFillStyle(1);
  g->SetFillColor(10);
  g->SetBorderSize(1);
  g->SetLineColor(10);
  //Usually legends should not have headers
  //but if you want one, uncomment the next line.
  g->SetHeader(legendheader);
}


// Here are a couple of other utility functions

// For a given histogram hst, return the number of entries between bin_lo and bin_hi
float get_nevents(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents += hst->GetBinContent(i);
  
  return nevents;         
}
// Partner function for above, returning the error for the above nevents
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents_err = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents_err += pow(hst->GetBinError(i),2);
  nevents_err = sqrt(nevents_err);
  
  return nevents_err;
}
