#include <err.h>
#include <getopt.h>
#include <unistd.h>

#include "Rtypes.h"
#include "TGeoManager.h"

#include "HepMC/HepMC.h"
#include "GunGenerator.h"
#include "HepMCGenerator.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "TTabPhysProcess.h"
#include "CMSApplication.h"

static int n_events   = 1;
static int n_buffered = 1;
static int n_threads  = 1;
static int n_track_max = 64;
static int n_learn_steps = 32768;
static int max_memory = 2048; /* MB */

static bool monitor = false, score = false, debug = false;

static struct option options[] = {
	{"events",      required_argument,  0,          'e'},
	{"hepmc-event-file", required_argument, 0,      'E'},
	{"fstate",      required_argument,  0,          'f'},
	{"geometry",    required_argument,  0,          'g'},
	{"learn-steps", required_argument,  0,          'l'},
	{"max-tracks-per-basket", required_argument, 0, 'B'},
	{"monitor",           no_argument,  0,          'm'},
	{"debug",             no_argument,  0,          'd'},
	{"max-memory",  required_argument,  0,          'M'},
	{"nbuffer",     required_argument,  0,          'b'},
	{"score",             no_argument,  0,          's'},
	{"threads",     required_argument,  0,          't'},
	{"xsec",        required_argument,  0,          'x'},
	{0, 0, 0, 0}
};

void help()
{
    printf("\nUsage: cmsapp [OPTIONS] INPUT_FILE\n\n");

    for (int i = 0; options[i].name != NULL; i++) {
        printf("\t-%c  --%s\t%s\n",
            options[i].val, options[i].name,
            options[i].has_arg ? options[i].name : "");
    }
    printf("\n\n");
}

int main(int argc, char *argv[])
{
	std::string cms_geometry_filename("cms2015.root");
	std::string xsec_filename("xsec_FTFP_BERT_G496p02_1mev.root");
	std::string fstate_filename("fstate_FTFP_BERT_G496p02_1mev.root");
	std::string hepmc_event_filename; // "pp14TeVminbias.root";

	if (argc == 1) { help(); exit(0); }

	while (true) {
		int c, optidx = 0;

		c = getopt_long(argc, argv, "E:e:f:g:l:B:mM:b:t:x:", options, &optidx);

		if (c == -1) break;

		switch (c) {
			case 0:
				c = options[optidx].val;
				/* fall through */
			case 'e':
				n_events = (int) strtol(optarg, NULL, 10);

				if (n_events <= 0)
					errx(1, "number of events must be positive");

				break;

			case 'E':
				hepmc_event_filename = optarg;
				break;

			case 'f':
				fstate_filename = optarg;
				break;

			case 'g':
				cms_geometry_filename = optarg;
				break;

			case 'l':
				n_learn_steps = (int) strtol(optarg, NULL, 10);

				if (n_learn_steps <= 0)
					errx(1, "number of learning steps must be positive");

				break;

			case 'B':
				n_track_max = (int) strtol(optarg, NULL, 10);

				if (n_track_max < 1)
					errx(1, "max number of tracks per basket must be positive");

				break;

			case 'm': monitor = true; break;

			case 'M':
				max_memory = (int) strtol(optarg, NULL, 10);

				if (max_memory < 128)
					errx(1, "max memory is too low");

				break;

			case 'b':
				n_buffered = (int) strtol(optarg, NULL, 10);

				if (n_buffered < 1)
					errx(1, "number of buffered events must be positive");

				break;

			case 't':
				n_threads = (int) strtol(optarg, NULL, 10);

				if (n_threads < 1)
					errx(1, "number of threads must be positive");

				break;

			case 's': score = true; break;

			case 'x':
				xsec_filename = optarg;
				break;

			default:
				errx(1, "unknown option %c", c);
		}
	}

	TGeoManager::Import(cms_geometry_filename.c_str());
	WorkloadManager *wmanager   = WorkloadManager::Instance(n_threads);
	GeantPropagator *propagator = GeantPropagator::Instance(n_events, n_buffered);

	wmanager->SetNminThreshold(5 * n_threads);

	propagator->fUseMonitoring = monitor;

	wmanager->SetMonitored(WorkloadManager::kMonQueue,          monitor);
	wmanager->SetMonitored(WorkloadManager::kMonMemory,         monitor);
	wmanager->SetMonitored(WorkloadManager::kMonBasketsPerVol,  monitor);
	wmanager->SetMonitored(WorkloadManager::kMonVectors,        monitor);
	wmanager->SetMonitored(WorkloadManager::kMonConcurrency,    monitor);
	wmanager->SetMonitored(WorkloadManager::kMonTracksPerEvent, monitor);

	// Initial vector size, this is no longer an important model parameter,
	// because is gets dynamically modified to accomodate the track flow
	propagator->fNperBasket = 16;   // Initial vector size

	// This is now the most important parameter for memory considerations
	propagator->fMaxPerBasket = n_track_max;

	// Maximum user memory limit [MB]
	propagator->fMaxRes = max_memory;

	propagator->fEmin = 0.001; // [10 MeV] energy cut
	propagator->fEmax = 0.01; // 10 MeV

	propagator->fProcess = new TTabPhysProcess("tab_phys", xsec_filename.c_str(), fstate_filename.c_str());

	if (hepmc_event_filename.empty()) {
		propagator->fPrimaryGenerator = new GunGenerator(propagator->fNaverage, 11, propagator->fEmax, -8, 0, 0, 1, 0, 0);
	} else {
		// propagator->fPrimaryGenerator->SetEtaRange(-2.,2.);
		// propagator->fPrimaryGenerator->SetMomRange(0.,0.5);
		// propagator->fPrimaryGenerator = new HepMCGenerator("pp14TeVminbias.hepmc3");
		propagator->fPrimaryGenerator = new HepMCGenerator(hepmc_event_filename);
	}

	propagator->fLearnSteps = n_learn_steps;

	CMSApplication *CMSApp = new CMSApplication();

	if (score) {
		CMSApp->SetScoreType(CMSApplication::kScore);
	} else {
		CMSApp->SetScoreType(CMSApplication::kNoScore);
	}

	propagator->fApplication = CMSApp;

	if (debug) {
		propagator->fUseDebug = kTRUE;
		propagator->fDebugTrk = 1;
	}

	propagator->fUseMonitoring = monitor;
	propagator->PropagatorGeom(cms_geometry_filename.c_str(), n_threads, monitor);

	return 0;
}