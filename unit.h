// internal units: Angstrom, GPa, atomic mass, e
// temperature to be expressed as thermal energy

const double meter = 1E10;
const double pascal = 1.e-9;
const double kg = 1./1.66053906660e-27;
const double coulomb = 1./1.602176634e-19;

const double joule = 1e21;
const double evolt = joule/coulomb;
const double hartree = 4.3597447222071e-18*joule;

const double second = meter*sqrt(kg/joule);		// [time] = 1.28862e-13 s
const double hPlanck = 6.62607015e-34*joule*second;
const double hBar = hPlanck/(2*M_PI);			// = 1.22194 in our units
const double aBohr = 5.29177210903e-11*meter;

const double inv4PiEps0 = hartree*aBohr;
const double epsilon0 = 1./(4*M_PI*inv4PiEps0);

const double kBoltzK = 1.380649e-23*joule;		// [temp] = 72.4297 kB*K
const double roomTemperature = 300*kBoltzK;		// = 4.14195 in our units

const double newton = joule/meter;
