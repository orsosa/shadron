{
  TTree *t=(TTree*) _file0->Get("outdata");

  t->SetAlias("xb","Q2/0.93827/2/Nu");
  TCut  KinCut = "(10.6-Nu)>2&&fE[0]>1&&fE[1]>1&&fE[0]/Nu>0.1&&fE[1]/Nu>0.1&&(fE[0]+fE[1])/Nu<0.95&&xFo>0&&theta0>10&&theta1>10&&Nu/10.6<0.8&&vzec<15&&1.05<sqrt(Mx2)";
  
  TCut DIS = "Q2>1&&W>2";
  TH1F *hphi_h0= new TH1F("hphi_h0","#phi_{R} helic 0",10,-200,200);
  hphi_h0->SetLineWidth(2);
  hphi_h0->GetXaxis()->SetTitle("#phi_{R}");
  hphi_h0->GetYaxis()->SetTitle("dN/d#phi_{R}");
  

  TH1F *hphi_h1= (TH1F*)hphi_h0->Clone("hphi_h1");
  hphi_h1->SetTitle("#phi_{R} helic 1");

  TH1F *hALU= (TH1F*) hphi_h0->Clone("hALU"); // diff/sum
  TH1F *hALU_s= (TH1F*) hphi_h0->Clone("hALU_s"); // sum
  TH1F *hALU_d= (TH1F*) hphi_h0->Clone("hALU_d"); // difference

  hALU->SetTitle("(N^{+} - N^{-})/(N^{+} + N^{-})/PB");
  hALU_s->SetTitle("(N^{+} + N^{-})");
  hALU_d->SetTitle("(N^{+} - N^{-})");

  hphi_h1->Sumw2(kFALSE);
  hphi_h0->Sumw2(kFALSE);

  ////////// bin 
  Float_t bin_low=0.65,bin_high=0.8;
  TCut binCut = Form("%f<M&&M<%f",bin_low,bin_high);
  /////////////////////////////////

  
  // t->Draw("phiR>>hphi_h0",DIS&&KinCut&&binCut&&"helic==0&&cos_theta_P0cm<0");
  // t->Draw("phiR>>hphi_h1",DIS&&KinCut&&binCut&&"helic==1&&cos_theta_P0cm<0");

  t->Draw("phiR>>hphi_h0",DIS&&KinCut&&binCut&&"helic==0");
  t->Draw("phiR>>hphi_h1",DIS&&KinCut&&binCut&&"helic==1");

  hphi_h1->Sumw2();
  hphi_h0->Sumw2();
  Float_t PB=0.8;

  hALU_s->Add(hphi_h1,hphi_h0,1,1);
  hALU_d->Add(hphi_h1,hphi_h0,1,-1);
  hALU->Divide(hALU_d,hALU_s,1,PB);

  TF1 *ff = new TF1("ff","[0]*sin( x*TMath::DegToRad() )",-180,180);

  gStyle->SetOptFit(111);
  ff->SetParName(0,"A_{LU}^{sin(#phi_{R})}");
  hALU->Draw("");
  hALU->Fit(ff);
  
  
  
}
