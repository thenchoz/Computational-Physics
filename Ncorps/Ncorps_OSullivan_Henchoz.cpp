#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <valarray>
#include <vector>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template


using namespace std;

typedef valarray<double> Position; // 4 dimension, sous la forme {x, y, vx, vy}

struct Corps
{
	string name;
	double masse;
	Position pos;
};


class Exercice4
{
  private:
	double t,dt,tFin;
	bool type_dt;
	double epsilon;
	size_t corps;
	vector<Corps> planete;
	double G;
	int sampling;
	int last;
	ofstream *outputFile;

	double d_position(bool if_x, Position const& p1, Position const& p2) const
	{
		if(if_x) return (p2 - p1)[0];
		return (p2 - p1)[1];
	};

	double distance(Position const& p1, Position const& p2) const
	{
		return sqrt(pow(d_position(true, p1, p2), 2) + pow(d_position(false, p1, p2), 2));
	};

	void printOut(bool force)
	{
		if((!force && last>=sampling) || (force && last!=1))
		{
			double ekin_tot(0);
			double epot_tot(0);

			*outputFile << t << ' ' << dt << ' ';

			for(auto cor : planete)
			{
				double ekin(cor.masse/2.0 * (pow(cor.pos[2], 2) + pow(cor.pos[3], 2)));
				double epot(0);

				for(auto gravit : planete)
				{
					if(cor.name != gravit.name)
					{
						epot -= G * cor.masse * gravit.masse / distance(cor.pos, gravit.pos);
					}
				}

				ekin_tot += ekin;
				epot_tot += epot;

				// *outputFile << cor.name << ' '; // pas lisible pour MatLab

				for(auto p : cor.pos)
				{
					*outputFile << p << ' ';
				}

				*outputFile << ekin << ' ' << epot << ' ';
			}

			*outputFile << ekin_tot << ' ' << epot_tot << endl;

			last=1;
		}else{
			last++;
		}
	};

	// force exercee sur le corps 1 par le corps 2, en composante if_x
	double Gravitation(bool if_x, double m1, double m2, Position const& p1, Position const& p2) const
	{
		return G * m1 * m2* d_position(if_x, p1, p2) / pow(distance(p1, p2), 3);
	};

	// force exercee par tous les autres corps sur le corps a la position p, en composante if_x
	double Somme_Gravitation(bool if_x, vector<Corps> const& corps, size_t p) const
	{
		double g(0);

		for(size_t i(0); i < corps.size(); ++i)
		{
			if(i != p)
				g += Gravitation(if_x, corps[p].masse, corps[i].masse, corps[p].pos, corps[i].pos);
		}

		// divise par sa masse pour avoir l'acceleration
		return g/corps[p].masse;
	};

	void f(Position& k, double dt, vector<Corps> const& planet, size_t i) const
	{
		k[2] = dt * Somme_Gravitation(true , planet, i);
		k[3] = dt * Somme_Gravitation(false, planet, i);
		k[0] = dt * planet[i].pos[2];
		k[1] = dt * planet[i].pos[3];
	};

	vector<Corps> step(vector<Corps> const& planet, double dt)
	{
		vector<Corps> retour(planet);

		vector<Corps> rk(planet);
		vector<Corps> rk_memory(planet);

		for(size_t i(0); i < corps; ++i)
		{
			// RK4 - K1
			Position k1(0., 4);
			f(k1, dt, rk, i);
			rk_memory[i].pos += k1/2.;

			retour[i].pos += k1/6.;
		}

		rk = rk_memory;

		for(size_t i(0); i < corps; ++i)
		{
			// RK4 - K2
			Position k2(0., 4);
			f(k2, dt, rk, i);
			rk_memory[i].pos = planet[i].pos + k2/2.;

			retour[i].pos += 2. * k2/6.;
		}

		rk = rk_memory;

		for(size_t i(0); i < corps; ++i)
		{
			// RK4 - K3
			Position k3(0., 4);
			f(k3, dt, rk, i);
			rk_memory[i].pos = planet[i].pos + k3;

			retour[i].pos += 2. * k3/6.;
		}

		rk = rk_memory;

		for(size_t i(0); i < corps; ++i)
		{
			// RK4 - K4
			Position k4(0., 4);
			f(k4, dt, rk, i);

			retour[i].pos += k4/6.;
		}

		return retour;
	};


  public:
	Exercice4(int argc, char* argv[])
	{
		string inputPath("configuration_systeme_complet.in"); // Fichier d'input par defaut
		if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice4 config_perso.in")
		inputPath = argv[1];

		ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.
		for(int i(2); i<argc; ++i) // Input complementaires ("./Onde config_perso.in input_scan=[valeur]")
			configFile.process(argv[i]);

		tFin	= configFile.get<double>("tFin");
		dt		= configFile.get<double>("dt");
		type_dt = configFile.get<bool>("type_dt");
		epsilon = configFile.get<double>("precision");
		corps	= configFile.get<size_t>("nb_corps");

		for(size_t i(1); i <= corps; ++i)
		{
			Corps cor;
			cor.name 	= configFile.get<string>("nom_corps" + to_string(i)).c_str();
			cor.masse 	= configFile.get<double>("masse_corps" + to_string(i));
			cor.pos.resize(4);
			for(size_t j(1); j <= 4; ++j)
			{
				string val("position" + to_string(j));
				val += "_corps";
				val += to_string(i);
				cor.pos[j - 1] = configFile.get<double>(val);
			}
			planete.push_back(cor);
		}

		G = 6.674 * pow(10,-11);
    // G = 1;

		sampling = configFile.get<int>("sampling");

		// Ouverture du fichier de sortie
		outputFile = new ofstream(configFile.get<string>("output").c_str());
		outputFile->precision(15);
	};

	~Exercice4()
	{
		outputFile->close();
		delete outputFile;
	};

	vector<Corps> dt_adapt(vector<Corps> const& planet, double& dt, double& t)
	{
		double d(0);
		vector<Corps> retour(step(planet,dt));
		vector<Corps> retourdt(step(step(planet, dt/2.), dt/2.));

		// stockage de la plus grand distance parmis toutes les planetes pour dt et dt/2, ainsi que toute les diff. de vitesse*dt
		for(size_t i(0); i < corps; ++i)
		{
			double d1(distance(retour[i].pos, retourdt[i].pos));
			if(d1 > d) d = d1;
		}

		// si la plus grande distance est <= a la precision
		if(d <= epsilon)
		{
			// augmentation du temps du pas de temps utilise
			t += dt;
			// modification du pas de temps
			dt *= pow(epsilon / d, 0.2);

			// retour du systeme apres un pas de temps dt precedent
			return retour;
		}

		// sinon, modification de dt et on recommence
		dt *= 0.99 * pow(epsilon / d, 0.2);
		return dt_adapt(planet, dt, t);
	};

	void run()
	{
		last = 0;
		t = 0.0;
		printOut(true);
		while( t<(tFin-0.5*dt))
		{
			// Si on veut utiliser un pas de temps adaptatif
			if(type_dt)
			{
				planete = dt_adapt(planete, dt, t);

				if(t + dt > tFin) dt = tFin - t;
			}
			// pas de temps standard
			else{
				planete = step(planete, dt);
				t += dt;

				if(t + dt > tFin) dt = tFin - t;
			}

			printOut(false);
		}
		printOut(true);
	};
};


int main(int argc, char* argv[])
{
	Exercice4 engine(argc, argv);
	engine.run();
	return 0;
}
