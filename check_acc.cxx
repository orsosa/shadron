{

  TNtuple *t = (TNtuple*) _file0->Get("ntuple_data");
  TString bname = "Zh";
  TString title_var="#phi_{H}";
  TString plot_var = "phiH";
  
  TH1F *hphi= new TH1F("hphi",title_var,400,-5,365);
  t->SetAlias("xb","Q2/0.93827/2/Nu");
  t->SetAlias("th_e","asin(sqrt(Q2/(10.6-Nu)/10.6))*TMath::RadToDeg()");
  t->SetAlias("phiH","(PhiPQ>0)*PhiPQ + (PhiPQ<0)*(PhiPQ+360)");
  TCut DCPip = "25<dcy_rot_0&&dcy_rot_0<135&&trajdczr0>0";
  TCut DCPim= "32<dcy_rot_0&&dcy_rot_0<145&&trajdczr0>0";
  
  TCut PidCut = "pid==211";
  TCut  KinCut = "Pe>2&&(P*P+0.14*0.14)>1.5&&Xf>0&&Nu/10.6<0.8&&Ze<15&&sqrt(W2p)>1.1&&10<th_e&&th_e<33" + PidCut;

  TF1 *ff = new TF1("ff","[0]/2 + [1]*sin(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad())",-5,360);
  ff->SetParameters(1,0.5,0.1);
  TCut DIS = "Q2>1&&W>2";

  TString be_str = "0.14 0.275 0.41 0.545 0.68 0.815 0.95";// bin edges z
  Int_t Nbe =  be_str.CountChar(' ') + 1;
  
  Double_t *be = new Double_t[Nbe];
  TString tok;
  Int_t k=0;
  Ssiz_t from = 0;
  while(be_str.Tokenize(tok,from," "))
  {
    be[k++] = atof(tok);
    cout<<tok<<endl;
  }

  TH1F *hC = new TH1F("hC_"+bname,"C over " + bname,Nbe-1,be);
  TH1F *hC2 = new TH1F("hC_"+bname,"C^{2} over " + bname,Nbe-1,be);
  TH1F *hB = new TH1F("hB_"+bname,"B over " + bname,Nbe-1,be);
  TH1F *hBC_2 = new TH1F("hBC_2_"+bname,"BC/2 over " + bname,Nbe-1,be);
  TH1F *hacc = new TH1F("hacc_"+bname,"acc over " + bname,Nbe-1,be);
  
  TCut bcut;
  for (int k=0;k<Nbe-1;k++)
  {
    ////////// bin /////////////
    Float_t bin_low=be[k], bin_high = be[k+1];
    TCut binCut = Form("%f<" + bname + "&&" + bname + "<%f",bin_low,bin_high);
    ///////////////////////////
    t->Draw(plot_var + ">>hphi",DIS&&KinCut&&binCut);
    hphi->Fit(ff,"S");
    
    Float_t C = ff->GetParameter(1)/ff->GetParameter(0);
    Float_t B = (ff->GetParameter(0) - ff->GetParameter(2))/ff->GetParameter(0);
    std::cout<<C<<"\t"<<B<<"\t"<<B*C/2.<<std::endl;
    
    hC->SetBinContent(k+1,C);
    hC2->SetBinContent(k+1,C*C);
    hB->SetBinContent(k+1,B);
    hBC_2->SetBinContent(k+1,B*C/2.);

      
  }

  hC->SetMarkerStyle(kFullDotLarge);
  hB->SetMarkerStyle(kFullDotLarge);
  hBC_2->SetMarkerStyle(kFullDotLarge);
  hacc->SetMarkerStyle(kFullDotLarge);

  //  hC->Draw("p");
  //hBC_2->Draw("psame");
  //hB->Draw("psame");
  hacc->Add(hB,hC2,1,-2.);
  hacc->Draw("p")
  

}
