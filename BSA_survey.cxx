////////////// Header can be modified //////////////////
///////////////// Testata può essere modificata /////////
TString bname = "Zh";
TString dirname= bname + "_bins";
//TString be_str = "0.2 0.48 0.64 0.8 2";// bin edges
TString be_str = "0.14 0.275 0.41 0.545 0.68 0.815 0.95";// bin edges z
//TString be_str = "0.15 0.35 0.55 1.0";// bin edges z
TCanvas c("c","c",960,640);

//TString title_var="#phi_{H} - #phi_{R}";
//TString plot_var = "phiH_phiR";
//TString title_var="#phi_{R}";
//TString plot_var = "phiR";
TString title_var="#phi_{H}";
TString plot_var = "phiH";

Float_t alu_hl = 0.1,alu_ll = -0.1;
TString HELIC_NEG = "-1";
///////////////////////////////////

void savepic (TString picname);
Int_t BSA_survey(TString fname="",TString tname="ntuple_data", Float_t field=-1){

  /////////// preamble: log file, bins, etc. ///////////
  /////////// preambolo: .../////////////
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
  
  time_t timestmp;
  time(&timestmp);
  
  struct tm *lt = localtime(&timestmp);
  TString date =Form("%02d%02d%02d_%02d%02d%02d",lt->tm_year%100,lt->tm_mon+1,lt->tm_mday,lt->tm_hour,lt->tm_min,lt->tm_sec);

  dirname += Form("_%02d%02d%02d",lt->tm_year%100,lt->tm_mon+1,lt->tm_mday);

  system("mkdir " + dirname + " 2>/dev/null");
  
  std::ofstream log(dirname + "/" + bname + "_" + date + ".log");

  log<<"BSA Analysis using bin on " + bname + "\n";
  log<<"edges: " + be_str <<std::endl;
  

  TH1F *hALU_bin = new TH1F("hALU_"+bname,"ALU over " + bname,Nbe-1,be);
  hALU_bin->GetXaxis()->SetTitle(bname);
  hALU_bin->GetYaxis()->SetTitle("A_{LU}^{sin(" + title_var + ")}");
  hALU_bin->GetYaxis()->SetRangeUser(alu_ll,alu_hl);
  ////////////////////////////////////////

  
  //  TFile fin(fname,"read"); 
  TChain *t = new TChain();
  t->Add(fname + "/" + tname);
  std::cout<<"Ntrees to be processed: "<<t->GetNtrees()<<std::endl;

  t->SetAlias("xb","Q2/0.93827/2/Nu");
  t->SetAlias("th_e","asin(sqrt(Q2/(10.6-Nu)/10.6))*TMath::RadToDeg()");
  t->SetAlias("phiH","(PhiPQ>0)*PhiPQ + (PhiPQ<0)*(PhiPQ+360)");
  TCut DCPip = "25<dcy_rot_0&&dcy_rot_0<135&&trajdczr0>0";
  TCut DCPim= "32<dcy_rot_0&&dcy_rot_0<145&&trajdczr0>0";
  
  //    TCut PidCut = "pid==-211" + DCPim;
  //  TCut PidCut = "pid==211" + DCPip;
  TCut PidCut = "pid==211";
  TCut thCut = (field<0)?"9.5<th_e&&th_e<33":"6.5<th_e&&th_e<28";
  TCut  KinCut = "Pe>2&&(P*P+0.14*0.14)>1.5&&Xf>0&&Nu/10.6<0.8&&Ze<15&&sqrt(W2p)>1.1" + thCut + PidCut;

  //  TCut  KinCut = "sqrt(Mx2)>1.05&&10<th_e&&th_e<30&&xb<0.1";
  
  //&&theta0>10&&theta1>10
  log<<"KinCut: " + TString(KinCut.GetTitle()) << std::endl;
  
  TCut DIS = "Q2>1&&W>2";
  TH1F *hphi_h0= new TH1F("hphi_h0",title_var + " helic -1",10,-5,365);
  hphi_h0->SetLineWidth(2);
  hphi_h0->GetXaxis()->SetTitle(title_var);
  hphi_h0->GetYaxis()->SetTitle("dN/d#phi");

  TH1F *hphi_h1= (TH1F*)hphi_h0->Clone("hphi_h1");
  hphi_h1->SetTitle(title_var + " helic 1");
  hphi_h0->SetLineColor(kRed);
  TLegend *leg = new TLegend(0.3,0.3,0.2,0.4);
  leg->AddEntry( hphi_h1,"helic 1 ");
  leg->AddEntry( hphi_h0,"helic -1");
    
  TH1F *hALU= (TH1F*) hphi_h0->Clone("hALU"); // diff/sum
  TH1F *hALU_s= (TH1F*) hphi_h0->Clone("hALU_s"); // sum
  TH1F *hALU_d= (TH1F*) hphi_h0->Clone("hALU_d"); // difference

  hALU->SetTitle("(N^{+} - N^{-})/(N^{+} + N^{-})/PB");
  hALU_s->SetTitle("(N^{+} + N^{-})");
  hALU_d->SetTitle("(N^{+} - N^{-})");

  hALU->GetYaxis()->SetTitle("A_{LU}");

  hphi_h1->Sumw2(kFALSE);
  hphi_h0->Sumw2(kFALSE);

  TH2F *h2sinp = new TH2F("h2sinp","sin+  vs " + bname,Nbe-1,be,300,-1,1);
  TH2F *h2sinm = new TH2F("h2sinm","-sin- vs " + bname,Nbe-1,be,300,-1,1);


  for (int k = 0;k<Nbe-1;k++)
  {
    ////////// bin 
    Float_t bin_low=be[k], bin_high = be[k+1];
    TCut binCut = Form("%f<" + bname + "&&" + bname + "<%f",bin_low,bin_high);
    /////////////////////////////////
  
    t->Draw(plot_var + ">>hphi_h0",DIS&&KinCut&&binCut&&("helic=="+HELIC_NEG));
    t->Draw(plot_var + ">>hphi_h1",DIS&&KinCut&&binCut&&"helic==1");


    hphi_h1->Sumw2();
    hphi_h0->Sumw2();
    Float_t PB=0.84;
    
    hALU_s->Add(hphi_h1,hphi_h0,1,1);
    hALU_d->Add(hphi_h1,hphi_h0,1,-1);
    hALU->Divide(hALU_d,hALU_s,1,PB);
    
    TF1 *ff = new TF1("ff","[0]*sin( x*TMath::DegToRad() )",-5,360);
    //TF1 *ff = new TF1("ff","[0]*sin( x*TMath::DegToRad() ) + [1]*sin( 2*x*TMath::DegToRad() )",0,360);

    
    gStyle->SetOptFit(111);
    ff->SetParName(0,"A_{LU}^{sin(" + title_var + ")}");
    //    ff->SetParName(1,"A_{LU}^{sin(2(" + title_var + "))}");
    hALU->Draw("");
    hALU->Fit(ff);
    savepic(Form("hALU_bin%d",k));

    hALU_bin->SetBinContent(k+1,ff->GetParameter(0));
    hALU_bin->SetBinError(k+1,ff->GetParError(0));

    hphi_h1->Draw();
    hphi_h0->Draw("same");
    leg->Draw();
    savepic(Form(plot_var + "_bin%d",k));
    
    log<<"bin"<<k<<" --> ["<<bin_low <<","<<bin_high<<"]\n";
  }
  gStyle->SetOptStat(0);
  hALU_bin->SetMarkerStyle(kFullDotLarge);
  hALU_bin->Draw();
  savepic("hALU_" + bname);

  Int_t Ne=-1;
  Ne = t->Draw("sin((" + plot_var + ")*TMath::DegToRad()):"+bname+">>h2sinp",DIS&&KinCut&&"helic==1");
  std::cout<<"Ne sinp "<< Ne <<std::endl;
  Ne = t->Draw("-sin((" + plot_var + ")*TMath::DegToRad()):"+bname+">>h2sinm",DIS&&KinCut&&("helic==" + HELIC_NEG));
  std::cout<<"Ne sinm "<< Ne <<std::endl;
  h2sinp->Draw();
  h2sinm->Draw();

  TH1F *m_sinp = (TH1F*)hALU_bin->Clone("m_sinp");
  TH1F *m_sinm = (TH1F*)hALU_bin->Clone("m_sinm");
  TH1F *sum_sin = (TH1F*)hALU_bin->Clone("sum_sin");
  
  m_sinp->SetTitle("<sin>^{+}");
  m_sinm->SetTitle("-<sin>^{-}");
  sum_sin->SetTitle("<sin>^{+} - <sin>^{-}");
    
  m_sinp = (TH1F *)h2sinp->ProfileX();
  m_sinm = (TH1F *)h2sinm->ProfileX();
  
  sum_sin->Add(m_sinp,m_sinm);
  sum_sin->GetYaxis()->SetRangeUser(alu_ll,alu_hl);
  sum_sin->Draw();
  savepic("sum_sin_" + bname);

  Double_t fit_m,fit_s,sum_m,sum_s;
  fit_m = hALU_bin->IntegralAndError(1,hALU_bin->GetNbinsX(),fit_s)/hALU_bin->GetNbinsX();
  sum_m = sum_sin->IntegralAndError(1,sum_sin->GetNbinsX(),sum_s)/sum_sin->GetNbinsX();

  fit_s = fit_s/sqrt(hALU_bin->GetNbinsX()-1);
  sum_s = sum_s/sqrt(sum_sin->GetNbinsX()-1);

  log<<"\n\n#type\tmean\trms"<<std::endl;
  log<<"fit\t"<< fit_m <<"\t"<< fit_s <<std::endl;
  log<<"sum_sin\t"<< sum_m <<"\t"<< sum_s <<std::endl;

  TFile outf(dirname + "/hALU_" + bname + ".root","recreate");
  hALU_bin->Write("",TObject::kOverwrite);
  h2sinp->Write("",TObject::kOverwrite);
  h2sinm->Write("",TObject::kOverwrite);
  m_sinp->Write("",TObject::kOverwrite);
  m_sinm->Write("",TObject::kOverwrite);
  sum_sin->Write("",TObject::kOverwrite);

  
  log.close();  
  return 0;
}

void savepic (TString picname)
{
  c.SaveAs(dirname + "/" + picname + ".gif");
  c.SaveAs(dirname + "/" + picname + ".C");
  c.SaveAs(dirname + "/" + picname + ".pdf");

}
