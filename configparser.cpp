#include "configparser.h"
#include "tinyxml.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>

int parseTime(const char * sTime);

bool parseConfigFile(const char * xmlFile)
{
    TiXmlDocument doc(xmlFile);
    if(!doc.LoadFile())
    {
        std::cout<<"Could not open ["<<xmlFile<<"]"<<std::endl;
        return false;
    }
    TiXmlHandle hDoc(&doc);
    TiXmlElement * pElem;
    TiXmlHandle hRoot(0);
    TiXmlHandle hNode(0);

    // root
    {
        pElem = hDoc.FirstChildElement("simrx").Element();
        hRoot = TiXmlHandle(pElem);
    }
    // scanned arc
    {
	pElem = hRoot.FirstChildElement("scanned-arc").Element();
	scannedarccommon_.scannedangle = atof( pElem->Attribute("angle") );
    }
    // projections
    {
	pElem = hRoot.FirstChildElement("projections").Element();
	projectionCount = atoi( pElem->Attribute("count") );
    }
    // state
    {
	pElem = hRoot.FirstChildElement("state").Element();
	statecommon_.savestateevery = atof( pElem->Attribute("save-every") );
    }
    // process
    {
	pElem = hRoot.FirstChildElement("process").Element();
	int shutdownAfter = parseTime(pElem->Attribute("shutdown-after"));
	processcommon_.shutdownafter = shutdownAfter;
    }
    // detector
    {
	pElem = hRoot.FirstChildElement("detector").Element();
	detectorcommon_.detectorwidth = atof( pElem->Attribute("width") );
	detectorcommon_.detectorheight = atof( pElem->Attribute("height") );
	detectorcommon_.detectordistance = atof( pElem->Attribute("distance") );
    }
    // image
    {
        pElem = hRoot.FirstChildElement("image").Element();

        // resolution
        hNode = pElem->FirstChildElement("spatial-resolution");
        imagecommon_.imageheight = atoi( hNode.Element()->Attribute("height") );
        imagecommon_.imagewidth = atoi( hNode.Element()->Attribute("width") );
        // gray levels
        hNode = pElem->FirstChildElement("intensity-resolution");
        imagecommon_.bitresolution = atoi( hNode.Element()->Attribute("bit-depth") );
    }
    // beam
    {
       pElem = hRoot.FirstChildElement("beam").Element();

       if( strcmp(pElem->Attribute("type"), "conic") == 0 )
       {
         // conic beam
         beamcommon_.beamtype = CONIC_BEAM;
         beamcommon_.halfangle = atof( pElem->FirstChildElement("geometry")->Attribute("half-angle") );
         beamcommon_.beamdistance = atof( pElem->FirstChildElement("geometry")->Attribute("distance") );
       } 
       else if( strcmp(pElem->Attribute("type"), "parallel") == 0 )
       {
         // parallel beam
         beamcommon_.beamtype = PARALLEL_BEAM;
         beamcommon_.beamwidth = atof( pElem->FirstChildElement("geometry")->Attribute("width") );
         beamcommon_.beamheight = atof( pElem->FirstChildElement("geometry")->Attribute("height") );
         beamcommon_.beamdistance = atof( pElem->FirstChildElement("geometry")->Attribute("distance") );
       }
       else
       {
          // unknown geometry type
          return false;
       }
       beamcommon_.npmax = atof( pElem->FirstChildElement("photons")->GetText() );
       beamcommon_.emax = atof( pElem->FirstChildElement("energy")->GetText() );

    }
    // sample geometry
    {
        pElem = hRoot.FirstChildElement("geometry").Element();
        strcpy(geometryfilecommon_.geometryfile, pElem->Attribute("file") );
    }
    // materials
    {
        pElem = hRoot.FirstChildElement("materials").Element();
        int i = 0;
        for(pElem = pElem->FirstChildElement("material"); pElem != NULL; pElem = pElem->NextSiblingElement())
        {
                 strcpy(materialcommon_.pmfile[i], pElem->FirstChildElement("file")->GetText());
                printf("PARSER: MATERIAL: %s\n",materialcommon_.pmfile[i]); 
		csimpa_.eabs[i][ELECTRON] = atof( pElem->FirstChildElement("Eabs-electrons")->GetText() );
                 csimpa_.eabs[i][PHOTON] = atof( pElem->FirstChildElement("Eabs-photons")->GetText() );
                 csimpa_.eabs[i][POSITRON] = atof( pElem->FirstChildElement("Eabs-positrons")->GetText() );
                 csimpa_.c1[i] = atof( pElem->FirstChildElement("C1")->GetText() );
                 csimpa_.c2[i] = atof( pElem->FirstChildElement("C2")->GetText() );
                 csimpa_.wcc[i] = atof( pElem->FirstChildElement("Wcc")->GetText() );
                 csimpa_.wcr[i] = atof( pElem->FirstChildElement("Wcr")->GetText() );
                 dsmaxcommon_.dsmax[i+1] = atof( pElem->FirstChildElement("dSmax")->GetText() );
                 i++;
        }
        materialcommon_.nmat = i;
    }

	return true;
}

int parseTime(const char * sTime)
{
	char * sTimeCopy = strdup(sTime);
	char pattern[3] = "-:";
	char * token;	
	int days = 0;
	int hours = 0;
	int minutes = 0;
	
	days = atoi( strtok(sTimeCopy, pattern) );
	hours = atoi( strtok(NULL, pattern) );
	minutes = atoi( strtok(NULL, pattern) );
	free(sTimeCopy);
	printf("Shutdown after: days=[%d] hours=[%d] minutes=[%d]\n", days, hours, minutes);
	return days * 86400 + hours * 3600 + minutes * 60; 
}
