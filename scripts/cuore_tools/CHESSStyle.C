{

  //  gSystem->Load("/Users/snoplus/Work/SNO/qsno2rat/libqsno2rat.so");

  // Adopted from BaBar collaboration
  TStyle *CheSSStyle= new TStyle("CHESS","CHESS approved plots style");

  // use plain black on white colors
  CheSSStyle->SetFrameBorderMode(0);
  CheSSStyle->SetCanvasBorderMode(0);
  CheSSStyle->SetPadBorderMode(0);
  CheSSStyle->SetPadColor(0);
  CheSSStyle->SetCanvasColor(0);
  CheSSStyle->SetStatColor(0);
  CheSSStyle->SetStatBorderSize(1);
  CheSSStyle->SetFillColor(0);
  CheSSStyle->SetLegendBorderSize(1);
  CheSSStyle->SetTextSize(0.06);

  // set the paper & margin sizes
  CheSSStyle->SetPaperSize(20,26);
  CheSSStyle->SetPadTopMargin(0.15);
  CheSSStyle->SetPadRightMargin(0.05);
  CheSSStyle->SetPadBottomMargin(0.16);
  CheSSStyle->SetPadLeftMargin(0.2);

  // use large Times-Roman fonts
  CheSSStyle->SetTextFont(132);
  CheSSStyle->SetTextSize(0.08);
  CheSSStyle->SetNdivisions(505,"x");
  CheSSStyle->SetLabelFont(132,"x");
  CheSSStyle->SetLabelFont(132,"y");
  CheSSStyle->SetLabelFont(132,"z");
  CheSSStyle->SetLabelSize(0.08,"x");
  CheSSStyle->SetTitleSize(0.08,"x");
  CheSSStyle->SetNdivisions(505,"y");
  CheSSStyle->SetLabelSize(0.08,"y");
  CheSSStyle->SetTitleSize(0.08,"y");
  CheSSStyle->SetLabelSize(0.08,"z");
  CheSSStyle->SetTitleSize(0.08,"z");
  CheSSStyle->SetLabelFont(132,"t");
  CheSSStyle->SetTitleFont(132,"x");
  CheSSStyle->SetTitleFont(132,"y");
  CheSSStyle->SetTitleFont(132,"z");
  CheSSStyle->SetTitleFont(132,"t");
  CheSSStyle->SetTitleFillColor(0);
  CheSSStyle->SetTitleX(0.1);
  CheSSStyle->SetTitleFontSize(0.08);
  CheSSStyle->SetTitleFont(132,"pad");

  // use bold lines and markers
  CheSSStyle->SetMarkerStyle(20);
  CheSSStyle->SetMarkerSize(1.0);
  CheSSStyle->SetHistLineWidth(1.85);
  CheSSStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  CheSSStyle->SetErrorX(0.00001);



  // do not display any of the standard histogram decorations
  // CheSSStyle->SetOptTitle(0);
  CheSSStyle->SetOptStat(0);
  CheSSStyle->SetOptFit(111);

  // put tick marks on top and RHS of plots
  CheSSStyle->SetPadTickX(1);
  //CheSSStyle->SetPadTickY(1);


  // Add a greyscale palette for 2D plots
  int ncol=50;
  double dcol = 1./float(ncol);
  double gray = 1;
  TColor **theCols = new TColor*[ncol];
  for (int i=0;i<ncol;i++) theCols[i] = new TColor(999-i,0.0,0.7,0.7);
  for (int j = 0; j < ncol; j++) {
    theCols[j]->SetRGB(gray,gray,gray);
    gray -= dcol;
  }
  int ColJul[100];
  for  (int i=0; i<100; i++) ColJul[i]=999-i;
  CheSSStyle->SetPalette(ncol,ColJul);

  // Define a nicer color palette (red->blue)
  // Uncomment these lines for a color palette (default is B&W)
  CheSSStyle->SetPalette(1,0);  // use the nice red->blue palette

  //Custom palette
  Double_t r[]    = {0.0, 0.0};
  Double_t g[]    = {0.0, 1.0};
  Double_t b[]    = {0.0, 1.0};
  // Double_t r[]    = {0.0, 0.0};
  // Double_t g[]    = {1.0, 0.0};
  // Double_t b[]    = {1.0, 0.0};
  Double_t stop[] = {0.0, .3};
  TColor::CreateGradientColorTable(2, stop, r, g, b, 100);

  // End of definition of CheSSStyle




}
