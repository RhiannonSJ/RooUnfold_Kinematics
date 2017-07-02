#include "roo_unfold.h"

using namespace std;

int mb_comparison(){


    TFile f_mb("/hepstore/rjones/Exercises/RooUnfold_Kinematics/working_dir/err_test/all-G16_02a.root");
    if(f_mb.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "=========================== MB file open ===========================" << endl;
    }
 
    TFile f_sbnd("/hepstore/rjones/Exercises/RooUnfold_Kinematics/working_dir/err_test/ddxsec_plot.root");
    if(f_sbnd.IsZombie()){
        std::cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "========================== SBND file open ==========================" << endl;
    }
 
    TH2D *h_mb   = (TH2D*) f_mb.Get("data:miniboone_nuccqe_2010:d2XS_dCosmu_dPmu:th2d");

    TH2D *h_sbnd = (TH2D*) f_sbnd.Get("h_ddxsec");

    gStyle->SetPalette(55);
    gStyle->SetNumberContours(250);

    TCanvas *c = new TCanvas();

    h_mb->SetStats(kFALSE);
    h_mb->Draw("colz");

    c->SaveAs("working_dir/err_test/mb_2D.png");

    TCanvas *c1 = new TCanvas();

    h_sbnd->SetStats(kFALSE);
    h_sbnd->Draw("colz");

    c1->SaveAs("working_dir/err_test/sbnd_2D.png");

    
    cout << " MB : " << h_mb->Integral("width");
    
    cout << ", SBND : " << h_sbnd->Integral("width");

    cout << ", SBND / MB : " << (h_sbnd->Integral("width") / h_mb->Integral("width")) << endl;

    int nx_mb   = h_mb->GetNbinsX();
    int nx_sbnd = h_sbnd->GetNbinsX();
    int ny_mb   = h_mb->GetNbinsY();
    int ny_sbnd = h_sbnd->GetNbinsY();

    if ( nx_mb != nx_sbnd || ny_mb != ny_sbnd ){
        
        cout << " nx mb : " << nx_mb << ", nx sbnd : " << nx_sbnd;
        cout << " | ny mb : " <<  ny_mb << ", ny sbnd : " << ny_sbnd << endl;

        exit(1);
    
    }

    // Get the slices
    // Y
    TCanvas *c_slice = new TCanvas();

    for ( int i = 1; i <= ny_mb; ++i ){
    
        double low_edge_y;
        double up_edge_y;

        low_edge_y = h_mb->GetYaxis()->GetBinLowEdge(i);
        up_edge_y = h_mb->GetYaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream conv;
        conv.clear();

        string title;
        title.clear();

        char file_name[1024];

        conv << setprecision(4) << "working_dir/err_test/mb_slices/Tmu_slice_" << low_edge_y << ";" << up_edge_y << ".png";
        title = conv.str();
        
        strcpy( file_name, title.c_str() );

        // For the title of the histogram
        stringstream hist;
        hist.clear();

        char hist_name[1024];
        string temp1;

        hist << setprecision(4) << "T_{#mu} slice: " << low_edge_y << "_" << up_edge_y;
        temp1 = hist.str();

        strcpy( hist_name, temp1.c_str() );

        TLegend * l_y = new TLegend(0.12, 0.68, 0.32, 0.88);
    
        TH1D *h_mb_slice    = h_mb->ProjectionX("h_mb_slice", i, i ); 
        TH1D *h_sbnd_slice  = h_sbnd->ProjectionX("h_sbnd_slice", i, i ); 
        
        l_y->AddEntry(h_mb_slice, "MiniBooNE","l");
        l_y->AddEntry(h_sbnd_slice, "SBND","l");

        
        for ( int j = 1; j <= nx_mb; ++j ){
            
            h_mb_slice->SetBinError(j, h_mb->GetBinError(j,i));
            h_sbnd_slice->SetBinError(j, h_sbnd->GetBinError(j,i));
        
        } 

        double max = 1.1 * TMath::Max(h_sbnd_slice->GetBinContent(h_sbnd_slice->GetMaximumBin()), h_mb_slice->GetBinContent(h_mb_slice->GetMaximumBin()));

        h_sbnd_slice->Draw("e1x0");
        h_mb_slice->Draw("e1x0same");
        h_sbnd_slice->SetTitle(hist_name);
        h_sbnd_slice->GetYaxis()->SetRangeUser(0,max);
        h_sbnd_slice->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        h_sbnd_slice->GetYaxis()->SetTitle("CC0#pi, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");
        h_sbnd_slice->SetLineColor(kBlue);
        h_mb_slice->SetLineColor(kRed);
        h_mb_slice->SetMarkerStyle(7);
        h_sbnd_slice->SetMarkerStyle(7);

        h_sbnd_slice->SetTitleOffset(1.2, "Y");
        h_sbnd_slice->SetStats(kFALSE);

        l_y->Draw();

        c_slice->SaveAs(file_name);

        delete h_mb_slice;
        delete h_sbnd_slice;
        delete l_y;
    } 

    // Get the slices
    // X
    TCanvas *c_slicex = new TCanvas();

    for ( int i = 1; i <= nx_mb; ++i ){
    
        double low_edge_x;
        double up_edge_x;

        low_edge_x = h_mb->GetXaxis()->GetBinLowEdge(i);
        up_edge_x = h_mb->GetXaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream convx;
        convx.clear();

        string titlex;
        titlex.clear();

        char file_namex[1024];

        convx << setprecision(4) << "working_dir/err_test/mb_slices/cosmu_slice_" << low_edge_x << ";" << up_edge_x << ".png";
        titlex = convx.str();
        
        strcpy( file_namex, titlex.c_str() );

        // For the title of the histogram
        stringstream histx;
        histx.clear();

        char hist_namex[1024];
        string temp1x;

        histx << setprecision(4) << "cos#theta_{#mu} slice: " << low_edge_x << "_" << up_edge_x;
        temp1x = histx.str();

        strcpy( hist_namex, temp1x.c_str() );

        TLegend * l_x = new TLegend(0.68, 0.68, 0.88, 0.88);
    
        TH1D *h_mb_slicex    = h_mb->ProjectionY("h_mb_slice", i, i ); 
        TH1D *h_sbnd_slicex  = h_sbnd->ProjectionY("h_sbnd_slice", i, i ); 
        
        l_x->AddEntry(h_mb_slicex, "MiniBooNE","l");
        l_x->AddEntry(h_sbnd_slicex, "SBND","l");

        
        for ( int j = 1; j <= ny_mb; ++j ){
            
            h_mb_slicex->SetBinError(j, h_mb->GetBinError(i,j));
            h_sbnd_slicex->SetBinError(j, h_sbnd->GetBinError(i,j));
        
        } 

        double maxx = 1.1 * TMath::Max(h_sbnd_slicex->GetBinContent(h_sbnd_slicex->GetMaximumBin()), h_mb_slicex->GetBinContent(h_mb_slicex->GetMaximumBin()));

        h_sbnd_slicex->Draw("e1x0");
        h_mb_slicex->Draw("e1x0same");
        h_sbnd_slicex->SetTitle(hist_namex);
        h_sbnd_slicex->GetYaxis()->SetRangeUser(0,maxx);
        h_sbnd_slicex->GetXaxis()->SetTitle("T_{#mu}");   
        h_sbnd_slicex->GetYaxis()->SetTitle("CC0#pi, d^{2}#sigma / dcos#theta_{#mu}dT_{#mu} [ 10^{-38} cm^{2} / GeV / n ]");
        h_sbnd_slicex->SetLineColor(kBlue);
        h_mb_slicex->SetLineColor(kRed);
        h_mb_slicex->SetMarkerStyle(7);
        h_sbnd_slicex->SetMarkerStyle(7);

        h_sbnd_slicex->SetTitleOffset(1.2, "Y");
        h_sbnd_slicex->SetStats(kFALSE);

        l_x->Draw();

        c_slicex->SaveAs(file_namex);

        delete h_mb_slicex;
        delete h_sbnd_slicex;
        delete l_x;
    } 

    delete c_slice;
    delete c_slicex;
    delete h_mb;
    delete h_sbnd;
    delete c;
    delete c1;

    return 0;
}
