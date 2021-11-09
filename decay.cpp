//Simulazione di un decadimento di un mesone in un pione + un kaone. Il mesone è ne sistema di riferimento del laboratorio, mentre i restanti nel canetro di massa.
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

#include "TRandom.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"

int main(int argc, char **argv){

	//dati del problema
	double m_B = 5.279; //GeV
	double m_pi = 0.140; //GeV
	double m_k = 0.500; //GeV
	double p_B = 0.300; //GeV (asse x)
	double p_cm = 2.614; //impulso dei prodotti nel centro di massa in GeV
	double mi; //massa invariante
	int N=10000;

	//apertura file
	TString rootfname("./output.root");
	TFile rfile(rootfname, "RECREATE");
	if(!rfile.IsOpen()){
		cout <<"Errore nell'apertura del file." << endl;
		exit(1);
	}

	//dichiarazione del vettore
	
	double angle;

	TLorentzVector p4_B, p4_pi, p4_k, p_t, p4_pi_O, p4_k_O; //qi del mesone
	p4_B.SetPxPyPzE(p_B,0,0,sqrt(p_B*p_B+m_B*m_B));
	TVector3 p1, p2;
	//istogramma
	TH1F h1("h1", "True Mass", 100, 5.2, 5.4);
	TH1F h2("h2", "Lab angles", 100, 2.9, 3.3);
	TH1F pi("pi", "Module", 100, 0, 6);
	


	//chiamata di TRandom3

	TRandom* gen = new TRandom();
	gen->SetSeed(0);
	double x, y, z;

	//dati per la simulazione del detector
	double p_k_O, p_pi_O, p_pi_meas, p_k_meas, resol = 0.03;

	for(int i=0; i < N; ++i){

			//punto casuale sulla sfera			
			gen->Sphere(x,y,z,p_cm);

			p4_pi.SetPxPyPzE(x,y,z,sqrt(p_cm*p_cm+m_pi*m_pi)); //qi pione

			p4_k.SetPxPyPzE(-x,-y,-z,sqrt(p_cm*p_cm+m_k*m_k)); //qi kaone
			
			//massa invariante
			p_t=p4_pi+p4_k;
			mi=sqrt(p_t.Dot(p_t)); 
			h1.Fill(mi);
			
			//boost
			p4_pi.Boost(p4_B.BoostVector());
			p4_k.Boost(p4_B.BoostVector());
		
			//calcoo angoli
			p1 = p4_pi.Vect();
			p2 = p4_k.Vect();

			angle = acos(p1.Dot(p2)/sqrt(p1.Dot(p1)*p2.Dot(p2)));
			h2.Fill(angle);

			//calcolo moduli (la loro distribuzione è uniforme)
			p_pi_O = sqrt(p1.Dot(p1));
			p_k_O = sqrt(p2.Dot(p2));

			//misura
			p_pi_meas = gen->Gaus(p_pi_O,resol*p_pi_O);
			p_k_meas = gen->Gaus(p_k_O,resol*p_k_O);
			pi.Fill(p_pi_meas);

	}			
	h1.Write();
	h2.Write();
	pi.Write();

	//plotting mi
	TCanvas canv1("canv1", "Canvas 1", 1280, 1024);
	h1.GetXaxis()->SetTitle("Distribution of invariant mass [GeV]");
	h1.Draw();
	canv1.SaveAs("true-mass.png");

	//plotting angle
	canv1.Clear();
	h2.GetXaxis()->SetTitle("Distribution of angles [rad]");
	h2.Draw();
	canv1.SaveAs("opening-angles.png");

	//plotting p_pi_meas
	canv1.Clear();
	pi.GetXaxis()->SetTitle("Distribution of pi [GeV]");
	pi.Draw();
	canv1.SaveAs("Momentum of pi after the measurement.png");

	//chiusure varie
	rfile.Close();
	delete gen;
	return 0;
}
