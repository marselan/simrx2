#ifndef SIMRX_GLOBAL_STRUCTS_H
#define SIMRX_GLOBAL_STRUCTS_H

#define ELECTRON 0
#define PHOTON   1
#define POSITRON 2

#define CONIC_BEAM 0
#define PARALLEL_BEAM 1

extern int projectionCount;

extern "C"
{
    extern struct statecommon_t
    {
	double savestateevery;
    } statecommon_;

    extern struct processcommon_t
    {
	int shutdownafter;
	int continueexecution;
    } processcommon_;

    extern struct scannedarccommon_t
    {
        double scannedangle;
    } scannedarccommon_;

    extern struct detectorcommon_t
    {
        double detectorwidth;
        double detectorheight;
        double detectordistance;
    } detectorcommon_;

    extern struct imagecommon_t
    {
        int imagewidth;
        int imageheight;
        int bitresolution;
    } imagecommon_;

    extern struct beamcommon_t
    {
        double beamtype; // beam geometry type: conic = 0,  parallel = 1
        double halfangle; // half angle of the conic beam (degrees)
        double beamdistance; // distance (cm) of the beam spot from the origin for both conic and parallel types
        double beamwidth; // for parallel geometry, width of the emisor (cm)
        double beamheight; // for parallel geometry, height of the emisot (cm)
        double npmax; // number of (primary) photons to be simulated
        double emax; // initial energy of primary photons
    } beamcommon_;

    extern struct geometryfilecommon_t
    {
        char geometryfile[20];
    } geometryfilecommon_;

    extern struct materialcommon_t
    {
        int nmat; // number of materials
        char pmfile[10][100]; // file names of the materials (.pmf files)
	} materialcommon_;

    extern struct csimpa_t
    {
        double eabs[10][3];
        double c1[10];
        double c2[10];
        double wcc[10];
        double wcr[10];
    } csimpa_;

    extern struct dsmaxcommon_t
    {
        double dsmax[11];
    } dsmaxcommon_;
}

#endif // SIMRX_GLOBAL_STRUCTS_H
