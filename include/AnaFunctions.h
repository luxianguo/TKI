#ifndef _ANAFUNCTIONS_
#define _ANAFUNCTIONS_

using namespace std;

namespace AnaFunctions
{

//avoid using Nature mass in GiBUU internal kinematic relation
Double_t MuonMass(){ return 105.65837/1e3; }//in GeV //google = wiki
Double_t NeutronMass(){ return 939.565/1e3;}//in GeV //wiki
Double_t ProtonMass(){ return 938.272/1e3;}//in GeV //google = wiki
Double_t PionMass(){ return 139.570/1e3;}//in GeV //wiki
Double_t KaonMass(){ return 493.677/1e3;}//in GeV //wiki
Double_t ElectronMass(){ return 0.510998/1e3;}//in GeV //wiki
Double_t PiZeroMass(){return 134.976/1e3;}//in GeV//wiki

int getTargetA(const int targetZ)
{
  if(targetZ == 1){
    return 1;
  }
  else if(targetZ == 6){
    return 12;
  }
  else if(targetZ == 8){
    return 16;
  }
  else if(targetZ == 18){
    return 40;
  }
  else if(targetZ == 26){
    return 56;
  }
  else if(targetZ == 82){
    return 207;
  }
  else{
    printf("\nNuGenTKI/include/AnaFunctions.h Unknown targetZ %d\n\n", targetZ); exit(1);
  }
}

double nuclearMass(const int targetA, const int targetZ)
{
  if(targetZ==18 && targetA==40){
    //eb = 343.8;// MeV, binding energy for Ar-40 https://www.wolframalpha.com/input/?i=Argon-40+binding+energy+*+40
    return 37.2247242;//GeV directly return mass https://www.wolframalpha.com/input/?i=Argon-40+mass+in+gev 
  }
  else{
    const double eb = 92.162;//eb in MeV, use C12 binding energy for the time being https://www.wolframalpha.com/input/?i=carbon-12+binding+energy+*12

    const double MA = (targetA-targetZ)*NeutronMass() + targetZ*ProtonMass() - eb/1E3;//GeV

    //checked output: ma 11.174860 mastar 10.262425
    //https://www.wolframalpha.com/input/?i=carbon-12+mass+in+gev 11.18
    return MA;
  }
}

double nuclearMassStar(const int targetA, const int targetZ)
{
  if(targetZ==18 && targetA==40){
    //https://www.wolframalpha.com/input/?i=Argon-39+mass+in+gev                       36.295028
    //https://www.wolframalpha.com/input/?i=chlorine-39+mass+in+gev                    36.2984698
    //Ar-39 or Cl-39 both have mass 36.3 GeV, use mean 36.2967489
    return 36.2967489;//GeV average of the two
  }
  else{  
    //see MAstar() below
    const double Bin=27.13/1E3;//GeV, use C11 excitation energy for the time being
    const double MAstar = nuclearMass(targetA, targetZ) - NeutronMass() + Bin; //GeV

    //printf("testpneb10down bin %f ", Bin);
    //tested by varing Bin 10% up and down, the pn in 0pi and piZERO has (10up-10down)/default varing less than 1%, see
    //https://docs.google.com/spreadsheets/d/1EiQ-W-YeHJF4STdY_Gs_H7yEr4gCL2TgjxUjxV6eWXM/edit?usp=sharing

    //checked output: ma 11.174860 mastar 10.262425
    //https://www.wolframalpha.com/input/?i=carbon-11+mass+in+gev 10.2570855
    return MAstar;
  }
}

/*
//to-do: need to distinguish nu and nubar (neutron or proton mass), also Bin !
double MAstar()
{
  //The values used in this paper are: $E_m^n=28.7$~MeV and $E_m^p=26.1$~MeV, see \cite{Bodek:2018lmc}, Table 8.
  //doesn't make sense to distinguish proton and neutron in GiBUU after all

  const double Bin=27.13/1E3;//GeV //average excitation energy instead of sampling it as in the paper
  const double MAstar = MA() - NeutronMass() + Bin; //GeV
  return MAstar;
}
*/

double Energy(const TLorentzVector * vec, const double mm)
{
  //don't trust internal energy or mass; only use experimental momentum

  const double pp = vec->P();

  return TMath::Sqrt(pp*pp+mm*mm);
}

double Ekin(const TLorentzVector * vec, const double mm)
{
  return Energy(vec, mm) - mm;
}

double EkinToP(const double mass, const double ek)
{
  const double energy = mass+ek;
  return TMath::Sqrt(energy*energy-mass*mass);
}

double EnuCCH(const TLorentzVector * mufull)
{
  const double mn = NeutronMass();
  const double mp = ProtonMass();
  const double mm = MuonMass();
  const double em = mufull->E();
  const double pz = mufull->Pz();

  const double num = mn*mn-mp*mp-mm*mm+2*mp*em;
  const double den = 2*(mp-em+pz)+1E-10;

  return num/den;
}

double GetTrueCCQEQ2(const double muonP, const double muonTheta)
{
  //theta in radian
  const double bindingE  = 34E-3;//all in GeV 
  const double muonE = sqrt( pow( muonP, 2 ) + pow( MuonMass(), 2 ) );
  const double nu_energy_num = pow( ProtonMass(), 2 ) - pow( NeutronMass() - bindingE, 2 ) - pow( MuonMass(), 2 ) + 2.0 * ( NeutronMass() - bindingE ) * muonE;
  const double nu_energy_den = 2.0 * ( NeutronMass() - bindingE - muonE + muonP * cos( muonTheta ) );
  const double nuE = nu_energy_num / (nu_energy_den+1E-10);

  const double q2 = 2.0 * nuE * ( muonE - muonP * cos( muonTheta ) ) - pow( MuonMass(), 2 );
  
  return q2;
}

double GetW2(const double phi, const double asy0)
{
  return 1-(TMath::Pi()/2.)*asy0*TMath::Sin(phi*TMath::DegToRad());
}

 double GetTwoBoostAdlerPhi(TLorentzVector nufull, TLorentzVector muonfull, TLorentzVector pifull, TLorentzVector nucleonfull, TLorentzVector iniNfull)
{
  const bool kprint = false;

  /*
first to boost everything to the initial nucleon rest frame, record the nu-mu, or equivalently the delta, direction (called v0);
then boost everything (except v0) to the delta rest frame and use the initial-nucleon-rest-frame v0 as z-axis.
Then y should still be nu-mu in the delta-rest-frame, it should be equal to v0 cross mu.
   */
  //lab frame
  TLorentzVector delta = pifull + nucleonfull;

  const TLorentzVector nuold = nufull;
  const TLorentzVector muonold = muonfull;
  const TLorentzVector piold = pifull;
  const TLorentzVector nucleonold = nucleonfull;
  const TLorentzVector iniNold = iniNfull;
  const TLorentzVector delold = delta;

  if(kprint){
    printf("\n\n\n==============================================================================================================================\n");
    printf("******************************* AnaFunctions::GetTwoBoostAdlerPhi in the Lab Frame\n");
    cout<<"iniNfull "<<endl; iniNfull.Print(); cout<<endl;
    cout<<"nufull "<<endl; nufull.Print(); cout<<endl;
    cout<<"muonfull "<<endl; muonfull.Print(); cout<<endl;
    cout<<"pifull "<<endl; pifull.Print(); cout<<endl;
    cout<<"nucleonfull "<<endl; nucleonfull.Print(); cout<<endl;
    cout<<"delta "<<endl; delta.Print(); cout<<endl;
  }

  //check 4-momentum conserved at vertex
  const TLorentzVector p4balance = nufull + iniNfull - muonfull - pifull - nucleonfull;
  if(p4balance.P()>1E-10){
    printf("\n\n\nAnaFunctions::GetTwoBoostAdlerPhi 4 momentum not conserved at vertes!\n\n\n"); exit(1);
    cout<<"iniNfull "<<endl; iniNfull.Print();
    cout<<"nufull "<<endl; nufull.Print();
    cout<<"muonfull "<<endl; muonfull.Print();
    cout<<"pifull "<<endl; pifull.Print();
    cout<<"nucleonfull "<<endl; nucleonfull.Print();
    exit(1);
  }

  //go to initial nucleon rest frame first
  const TVector3 boostToIniN = -iniNfull.BoostVector();

  nufull.Boost(boostToIniN);
  muonfull.Boost(boostToIniN);
  pifull.Boost(boostToIniN);
  nucleonfull.Boost(boostToIniN);
  delta.Boost(boostToIniN);
  iniNfull.Boost(boostToIniN);

  if(iniNfull.P()>1E-10){
    printf("\n\n\nAnaFunctions::GetTwoBoostAdlerAngle something wrong with iniN boosting!\n\n\n");
    iniNfull.Print();
    exit(1);
  }

  const TVector3 Zaxis = delta.Vect().Unit();//should be equal to nu-mu

  if(kprint){
    printf("\n\n\n******************************* AnaFunctions::GetTwoBoostAdlerPhi in the Initial Nucleon Rest Frame\n");
    cout<<"iniNfull "<<endl; iniNfull.Print(); cout<<endl;
    cout<<"nufull "<<endl; nufull.Print(); cout<<endl;
    cout<<"muonfull "<<endl; muonfull.Print(); cout<<endl;
    cout<<"pifull "<<endl; pifull.Print(); cout<<endl;
    cout<<"nucleonfull "<<endl; nucleonfull.Print(); cout<<endl;
    cout<<"delta "<<endl; delta.Print(); cout<<endl;
    cout<<"Zaxis "<<endl; Zaxis.Print(); cout<<endl;
  }

  //from iniN rest frame to delta rest frame
  const TVector3 boostToDelta = -delta.BoostVector();

  nufull.Boost(boostToDelta);
  muonfull.Boost(boostToDelta);
  pifull.Boost(boostToDelta);
  nucleonfull.Boost(boostToDelta);
  delta.Boost(boostToDelta);
  iniNfull.Boost(boostToDelta);

  //boost to delta rest frame, check delta at rest after boost
  if(delta.P()>1E-10){
    printf("\n\n\nAnaFunctions::GetTwoBoostAdlerAngle something wrong with boosting!\n\n\n"); 
    delta.Print(); 
    exit(1);
  }

  const TVector3 Yaxis = (nufull.Vect().Cross(muonfull.Vect())).Unit();
  const double yzdot = Yaxis.Dot(Zaxis);
  if(fabs(yzdot)>1E-10){
    printf("\n\n\nAnaFunctions::GetTwoBoostAdlerPhi Y and Z not perpendicular! %f\n\n\n", yzdot); 
    cout<<"iniNfull "<<endl; iniNfull.Print();
    cout<<"nuold "<<endl; nuold.Print();
    cout<<"muonold "<<endl; muonold.Print();
    cout<<"piold "<<endl; piold.Print();
    cout<<"nucleonold "<<endl; nucleonold.Print();
    cout<<"nufull "<<endl; nufull.Print();
    cout<<"muonfull "<<endl; muonfull.Print();
    cout<<"pifull "<<endl; pifull.Print();
    cout<<"nucleonfull "<<endl; nucleonfull.Print();
    cout<<"Yaxis "<<endl; Yaxis.Print();
    cout<<"Zaxis "<<endl; Zaxis.Print();
    exit(1);
  }

  const TVector3 Xaxis = Yaxis.Cross(Zaxis);

  const double nucleonX = nucleonfull.Vect().Dot(Xaxis);
  const double nucleonY = nucleonfull.Vect().Dot(Yaxis);
  const double nucleonR = TMath::Sqrt(nucleonX*nucleonX + nucleonY*nucleonY);

  //in 0-180
  double phi = TMath::ACos(nucleonX/nucleonR)*TMath::RadToDeg();
  if(phi<0 || phi>180){
    printf("\n\n\nAnaFunctions::GetTwoBoostAdlerPhi wrong domain in ACos %f %f %f\n\n\n", nucleonX, nucleonY, phi); exit(1);
  }
  if(nucleonY<0){
    phi = 360-phi;
  }

  if(kprint){
    printf("\n\n\n******************************* AnaFunctions::GetTwoBoostAdlerPhi in Delta Rest Frame Adler Phi %f\n", phi);
    cout<<"iniNfull "<<endl; iniNfull.Print(); cout<<endl;
    cout<<"nufull "<<endl; nufull.Print(); cout<<endl;
    cout<<"muonfull "<<endl; muonfull.Print(); cout<<endl;
    cout<<"pifull "<<endl; pifull.Print(); cout<<endl;
    cout<<"nucleonfull "<<endl; nucleonfull.Print(); cout<<endl;
    cout<<"delta "<<endl; delta.Print(); cout<<endl;
    cout<<"Xaxis "<<endl; Xaxis.Print(); cout<<endl;
    cout<<"Yaxis "<<endl; Yaxis.Print(); cout<<endl;
    cout<<"Zaxis "<<endl; Zaxis.Print(); cout<<endl;
    cout<<"nu-mu"<<endl; (nufull-muonfull).Vect().Print(); cout<<endl;
    printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Check this final nucleon 3 component in The Two-Boost Frame\n");
    const double nucleonZ = nucleonfull.Vect().Dot(Zaxis);
    const TVector3 tmp(nucleonX, nucleonY, nucleonZ);
    tmp.Print();
    printf("==============================================================================================================================\n");
  }

  return phi;
}

double GetOneBoostAdlerPhi(TLorentzVector nufull, TLorentzVector muonfull, TLorentzVector pifull, TLorentzVector nucleonfull, TLorentzVector iniNfull)
{
  const bool kprint = false;

  //same as two-boost, just differ a Wigner Rotation which doesn't change the relative position of particles and therefore Adler angles
  //lab frame
  const TLorentzVector delta = pifull + nucleonfull;

  if(kprint){
    printf("\n\n\n==============================================================================================================================\n");
    printf("******************************* AnaFunctions::GetOneBoostAdlerPhi in the Lab Frame\n");
    cout<<"iniNfull "<<endl; iniNfull.Print(); cout<<endl;
    cout<<"nufull "<<endl; nufull.Print(); cout<<endl;
    cout<<"muonfull "<<endl; muonfull.Print(); cout<<endl;
    cout<<"pifull "<<endl; pifull.Print(); cout<<endl;
    cout<<"nucleonfull "<<endl; nucleonfull.Print(); cout<<endl;
    cout<<"delta "<<endl; delta.Print(); cout<<endl;
  }

  const TVector3 boost = -delta.BoostVector();
  
  nufull.Boost(boost);
  muonfull.Boost(boost);
  nucleonfull.Boost(boost);

  const TVector3 Zaxis = (nufull-muonfull).Vect().Unit();//only after boost! has to be q direction, otherwise won't be perpendicular to nu cross mu when there if Fermi motion
  const TVector3 Yaxis = (nufull.Vect().Cross(muonfull.Vect())).Unit();
  const TVector3 Xaxis = Yaxis.Cross(Zaxis);

  const double nucleonX = nucleonfull.Vect().Dot(Xaxis);
  const double nucleonY = nucleonfull.Vect().Dot(Yaxis);


  //in 0-180
  double phi = TMath::ATan2(nucleonY, nucleonX)*TMath::RadToDeg();
  if(phi<0){
    phi+=360;
  }

  if(kprint){
    printf("\n\n\n******************************* AnaFunctions::GetOneBoostAdlerPhi in Delta Rest Frame Adler Phi %f\n", phi);
    cout<<"iniNfull "<<endl; iniNfull.Print(); cout<<endl;
    cout<<"nufull "<<endl; nufull.Print(); cout<<endl;
    cout<<"muonfull "<<endl; muonfull.Print(); cout<<endl;
    cout<<"pifull "<<endl; pifull.Print(); cout<<endl;
    cout<<"nucleonfull "<<endl; nucleonfull.Print(); cout<<endl;
    cout<<"delta "<<endl; delta.Print(); cout<<endl;
    cout<<"Xaxis "<<endl; Xaxis.Print(); cout<<endl;
    cout<<"Yaxis "<<endl; Yaxis.Print(); cout<<endl;
    cout<<"Zaxis "<<endl; Zaxis.Print(); cout<<endl;
    cout<<"nu-mu"<<endl; (nufull-muonfull).Vect().Print(); cout<<endl;
    printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Check this final nucleon 3 component in The One-Boost Frame\n");
    const double nucleonZ = nucleonfull.Vect().Dot(Zaxis);
    const TVector3 tmp(nucleonX, nucleonY, nucleonZ);
    tmp.Print();
    printf("==============================================================================================================================\n");
  }

  return phi;
}


double GetPseudoPhi(TLorentzVector nufull, TLorentzVector muonfull, TLorentzVector pifull, TLorentzVector nucleonfull)
{
  TLorentzVector delta = pifull + nucleonfull;
  const TVector3 Zaxis = delta.Vect().Unit();

  const TVector3 boost = -delta.BoostVector();

  nufull.Boost(boost);
  muonfull.Boost(boost);
  pifull.Boost(boost);
  nucleonfull.Boost(boost);
  delta.Boost(boost);

  const TVector3 Yaxis = Zaxis.Cross(muonfull.Vect().Unit());
  const TVector3 Xaxis = Yaxis.Cross(Zaxis);

  const double nucleonX = nucleonfull.Vect().Dot(Xaxis);
  const double nucleonY = nucleonfull.Vect().Dot(Yaxis);
  const double nucleonR = TMath::Sqrt(nucleonX*nucleonX + nucleonY*nucleonY);

  double phi = TMath::ACos(nucleonX/nucleonR)*TMath::RadToDeg();
  if(phi<0 || phi>180){
    printf("\n\n\nAnaFunctions::GetPseudoPhi wrong domain in ACos %f %f %f\n\n\n", nucleonX, nucleonY, phi); exit(1);
  }
  if(nucleonY<0){
    phi = 360-phi;
  }

  return phi;
}

double SmeardpTT(const double sigma)
{
  static TF1 * fran = 0x0;

  if(!fran){
    printf("SmeardpTT initializing fran %f\n", sigma);
    fran = new TF1("fran", "TMath::CauchyDist(x,0,[0])", -0.5, 0.5);
  }

  //sigma in MeV  
  fran->SetParameter(0, sigma/1e3);

  const double scale = sigma/20;
  //need to set range, otherwise the sampling is very coarse
  fran->SetRange(-0.5*scale, 0.5*scale);

  //Random(xmin, xmax) is not working if range is not set
  return fran->GetRandom();
}

double GetThetaRef(const TVector3 &vold, const TVector3 &vreftmp)//in deg
{
  const TVector3 vrefUnit = vreftmp.Unit();

  const double oldDotRef = vrefUnit.Dot(vold);
  const TVector3 pL = vrefUnit * oldDotRef;
  const TVector3 pT = vold - pL;

  double theta = TMath::ATan(pT.Mag()/oldDotRef);
  if(theta<0){
    theta += TMath::Pi();
  }

  return theta;
}

TVector3 getPtVect(const TLorentzVector * fullp, const TLorentzVector * basevect)
{
  const TVector3 unitvect = (basevect->Vect()).Unit();
  
  const TVector3 pLvect = unitvect * fullp->Vect().Dot(unitvect);

  return fullp->Vect() - pLvect;
}

double getRecoilP(const double beamP, const double dPT, const double pLFS)
{
  const double dpl = pLFS - beamP;
  return TMath::Sqrt(dpl*dpl+dPT*dPT);
}

double getRecoilM(const double beamMass, const double beamP, const double dPT, const double pLFS, const double eFS, const double m1)
{
  const double beamEnergy = TMath::Sqrt( beamMass*beamMass + beamP*beamP );
  const double iniNp = getRecoilP(beamP, dPT, pLFS);
  const double BB = eFS - m1;
  const double mxSq = TMath::Power(beamEnergy - BB, 2) - iniNp*iniNp;

  if(mxSq<0){
    //do not print. Too many for GiBUU because no nucleus is included
    //printf("AnaFunctions::getRecoilM mxSq<0 beamEnergy %f BB %f dPT %f mxSq %f\n", beamEnergy, BB, dPT, mxSq); 
    return -999;
    //exit(1);
  }
  else{
    return TMath::Sqrt(mxSq); 
  }
}

double getdPL(const double beamMass, const double dPT, const double pLFS, const double eFS, const double m1, const double m2)
{
  //original idea by J. Sobczyk, incorporated in 14 Nov 2017
  //https://journals.aps.org/prc/abstract/10.1103/PhysRevC.95.065501
  //Eq. 11
  //generalized to massive beam by Xianguo Lu
  //all in GeV  

   //printf("debug beamMass %f m1 %f m2 %f dpT %f pLFS %f eFS %f\n", beamMass, m1, m2, dPT, pLFS, eFS);
   
   const double AA = pLFS;
   const double BB = eFS - m1;
   const double CC = m2*m2 + dPT*dPT;

   const double aa = -AA/BB;
   const double bb = ( AA*AA - BB*BB - CC + beamMass*beamMass )/2/BB;

   const double delta = 4*aa*aa*bb*bb-4*(aa*aa-1)*(bb*bb-CC);

   if(delta<0){
     printf("AnaFunctions::getdPL delta < 0!! aa %f bb %f CC %f delta %f\n", aa, bb, CC, delta);
     exit(1);
   }
   
   const double sol1 = (-2*aa*bb+TMath::Sqrt(delta))/2/(aa*aa-1);
   const double sol2 = (-2*aa*bb-TMath::Sqrt(delta))/2/(aa*aa-1);

   const double lhs1 = aa*sol1 + bb;
   const double lhs2 = aa*sol2 + bb;

   const double beamP1 = AA - sol1;
   const double beamP2 = AA - sol2;

   double dpL = -999;
   int kpass = 0;
   if(lhs1>0 && beamP1>0){
     kpass++;
     dpL = sol1;
   }
   if(lhs2>0 && beamP2>0){
     kpass++;
     dpL = sol2;
   }

   if(kpass==1){
     return dpL;
   }
   else if(kpass==0){
     //do not print, too many for GiBUU because no nucleus is included
     //printf("AnaFunctions::getdPL *no* solution AA %f sol1 %f sol2 %f lhs1 %f lhs2 %f\n", AA, sol1, sol2, lhs1, lhs2);
     return -999;
   }
   else{
     printf("AnaFunctions::getdPL bad solution AA %f sol1 %f sol2 %f lhs1 %f lhs2 %f\n", AA, sol1, sol2, lhs1, lhs2);
     exit(1);
   }
}

void getCommonTKI(const int targetA, const int targetZ, const TLorentzVector *beamfullp, const TLorentzVector *scatterfullp, const TLorentzVector *recoilfullp, double & dalphat, double & dphit, double & dpt, double & dpTT, double & beamCalcP, double & IApN, double & recoilM, double & recoilP)
{
  //
  //note that this is for general calculation, all particle energy is sqrt(p^2+m^2)!
  //
  const TLorentzVector allFSfullp =  scatterfullp ? ((*recoilfullp)+(*scatterfullp)) : (*recoilfullp) ;
  const TVector3 vdPt = getPtVect(&allFSfullp, beamfullp);
  dpt = vdPt.Mag();

  if(scatterfullp){
    const TVector3 pTscatter = getPtVect(scatterfullp, beamfullp);
    const TVector3 pTrecoil  = getPtVect(recoilfullp, beamfullp);

    const TVector3 unitqt = -pTscatter.Unit();
    dphit   = TMath::ACos(pTrecoil.Dot(unitqt)/pTrecoil.Mag())*TMath::RadToDeg(); //in Deg 

    //dpt cutoff for hydrogen res dpt is 1E-5 
    if(dpt>1E-5){
      dalphat = TMath::ACos(    vdPt.Dot(unitqt)/vdPt.Mag()    )*TMath::RadToDeg(); //in Deg
    }
    else{//hydrogen
      dalphat=-999;
    }

    //if dpt<1E-5, then dpTT is independent of dalphat anyway
    dpTT = dpt * sin(dalphat*TMath::DegToRad());
    const Double_t dotcross = recoilfullp->Vect().Dot( (beamfullp->Vect()).Cross( scatterfullp->Vect() ));
    if(dotcross<0){
      dpTT *= -1;
    }
  }

  //previous calculation has a bug in primL, where by definition both are positive, which could in fact be negative.
  //effect 506/126897 = 0.004; 50 out of 506 causing neutronmomentum difference larger than 10%; mx difference is all smaller than 1%
  //only affects events with backward final states which are not in the previous calcuation for MINERvA mu and proton
  const double pLFS = allFSfullp.Vect().Dot(beamfullp->Vect().Unit());
  const double ma     = nuclearMass(targetA, targetZ);

  //printf("testbug  P %f E %f M %f\n", beamfullp->P(), beamfullp->E(), beamfullp->M());

  //use block to separate the variable definitions
  {//(1)---> without knowledge of the beam momentum/energy, only direction and mass
    const double mastar = nuclearMassStar(targetA, targetZ);  //only assume one nucleon removal
    //double getdPL(const double beamMass, const double dPT, const double pLFS, const double eFS, const double m1, const double m2)
    const double dpL = getdPL(beamfullp->M(), dpt, pLFS, allFSfullp.E(), ma, mastar);

    beamCalcP = -999;
    IApN = -999;
    if(dpL!=-999){
      beamCalcP = pLFS - dpL;
      IApN = TMath::Sqrt(dpL*dpL + dpt*dpt);//implus approximation, emulated nucleon momentum assuming single nucleon knock-out
    }
    //printf("testpn ma %f mastar %f pL %f IApN %f\n", ma, mastar, pL, IApN);
  }//(1)<---
    
  {//(2)---> knowing beam 4-momentum
    //double getRecoilM(const double beamMass, const double beamP, const double dPT, const double pLFS, const double eFS, const double m1)
    recoilM = getRecoilM(beamfullp->M(), beamfullp->P(), dpt, pLFS, allFSfullp.E(), ma);

    //double getRecoilP(const double beamP, const double dPT, const double pLFS)
    recoilP = getRecoilP(beamfullp->P(), dpt, pLFS);
  }//(2)<---
}

//end of namespace
}

#endif
