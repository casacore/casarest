//# tVisibilityIterator.cc: Tests the Synthesis MeasurementSet Iterator
//# Copyright (C) 1995,1999,2000,2001
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$

#include <msvis/MSVis/VisibilityIteratorAsync.h>
#include <msvis/MSVis/VisBufferAsync.h>
#include <casacore/tables/Tables.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/iostream.h>
#include <casacore/casa/iomanip.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <memory>

#include <casacore/tables/DataMan/ForwardCol.h>

#include <casacore/casa/namespace.h>

#define VIA ROVisibilityIteratorAsync
int
main(int argc, char **argv)
{
    // register forward col engine
    ForwardColumnEngine::registerClass();

    if (argc<2) {
        cout <<"Usage: tVisibilityIterator ms-table-name"<<endl;
        exit(1);
    }
    // try to iterate over a pre-existing MS table

    try {
        // read the MS from file
        cout << "Reading MS from file and checking type"<<endl;
        MeasurementSet synms(argv[1],Table::Update);

        cout << "Number of rows in ms: "<< synms.nrow()<<endl;

        cout << "Constructing Iterator"<<endl;
        Block<Int> bi(0); // create empty block with sortColumns

        {
            cout << "\n*** Test #1 ***\n" << endl;

            ROVisibilityIteratorAsync::PrefetchColumns prefetchColumns =
                    ROVisibilityIteratorAsync::prefetchColumns(VIA::Ant1, VIA::Ant2, VIA::Freq, VIA::Time,
							       VIA::ObservedCube, VIA::Sigma, VIA::Flag, VIA::Uvw,
                                                             -1);
            std::unique_ptr<ROVisibilityIterator> syniter (ROVisibilityIteratorAsync::create (synms,prefetchColumns, bi));
            VisBufferAutoPtr vb (syniter.get());

            // set iterator to start of data
            syniter->origin();

            cout << "Now testing data access functions"<<endl;
            cout << "Showing first 10 iterations:"<<endl;
            cout.flush();

            Double time0=vb->time()[0];
            for (Int i=0; i<10; i++) {
                cout << "Iteration "<<i<<" contains:"<<endl;
                cout << "Antenna1="<<vb->antenna1()[0];
                cout <<", Antenna2="<<vb->antenna2()[0];
                cout << setprecision(9);
                cout <<", freq="<<vb->frequency()[0];
                cout << setprecision(9);
                cout <<", time="<<vb->time()-time0;
                //				<<", uvw="<<vb->uvw();
                cout <<", vis="<<vb->visCube();
                cout <<", sigma="<<vb->sigma();
                cout <<", flag="<<vb->flag()<<endl;
                cout.flush();
                (* syniter) ++;
            }

            Double time=vb->time()[0];
            while( time == vb->time()[0]) (* syniter) ++;

            cout << "Next timeslot contains:"<<endl;
            cout << "Antenna1="<<vb->antenna1()[0]
                 <<",Antenna2="<<vb->antenna2()[0]
                 << setprecision(9)
                 <<", freq="<< vb->frequency()[0]
                 << setprecision(9)
                 <<", time="<< vb->time()[0]
                 //	 <<", uvw="<<vb->uvw()
                 <<", vis="<< vb->visCube()
                 <<", sigma="<<vb->sigma()
                 <<", flag="<<vb->flag()
                 <<", feed_pa="<< vb->feed_pa(vb->time()[0])<<endl;

            cout.flush();

            cout << " Reading through visibility data till end of file "<< endl;
            cout.flush();
            Complex tot=0;
            for (; syniter->more(); (* syniter) ++) {
                cout <<"."; tot+=sum(vb->visCube()); vb->uvw();
            }
            cout << " done : sum="<<tot<<endl;


        }

        {
            cout << "\n*** Test #2 ***\n" << endl;
            ROVisibilityIteratorAsync::PrefetchColumns prefetchColumns =
                    ROVisibilityIteratorAsync::prefetchColumns(VIA::Time, -1);


            cout << " Now try to iterate in time-intervals of 10s"<<endl;

            std::unique_ptr <ROVisibilityIterator>
            syniter2 (ROVisibilityIteratorAsync::create (synms,prefetchColumns, bi,10.));

            VisBufferAutoPtr vb2 (* syniter2);

            for (syniter2->originChunks();syniter2->moreChunks(); syniter2->nextChunk()) {
                syniter2->origin();
                cout <<" Interval start: "<< setprecision(10)<< vb2->time();
                Int n=0;
                Double time = 0;
                for (;syniter2->more(); (* syniter2) ++) {
                    n++;
                    time=vb2->time()[0];
                }
                cout << ". end: "<< setprecision(10)<< time;
                cout << ", # visibilities in interval: "<<n<<endl;
            }
        }

        {

            cout << "\n*** Test #3 ***\n" << endl;
            cout << " Try iterator with 1000s interval."<<endl;

            ROVisibilityIteratorAsync::PrefetchColumns prefetchColumns =
                    ROVisibilityIteratorAsync::prefetchColumns(VIA::Time, -1);

            std::unique_ptr <ROVisibilityIterator>
            syniter3 (ROVisibilityIteratorAsync::create(synms,prefetchColumns, bi,10.));

            syniter3->setInterval(1000.);
            VisBufferAutoPtr vb3 (* syniter3);

            for (int i = 0; i < 3; i++){ // see if rewinding works

                cout << "... Starting VI sweep #" << i << endl;

                for (syniter3->originChunks();syniter3->moreChunks(); syniter3->nextChunk()) {
                    syniter3->origin();
                    cout <<" Interval start: "<< setprecision(10)<< vb3->time()[0];
                    Int n=0;
                    //	vb.feed_pa(vb.time()[0]);
                    for (;syniter3->more(); (* syniter3) ++) {
                        n++;
                    }
                    cout << ", end: "<< setprecision(10)<< vb3->time();
                    cout << ", # visibilities in interval: "<<n<<endl;
                }
            }
        }


        //		{
        //
        //			ROVisibilityIteratorAsync syniter2(synms,bi,10.);
        //			VisBufferAsync vb2(syniter2);
        //
        //			ROVisibilityIteratorAsync syniter3 (synms,bi,10.);
        //			syniter3.setInterval(1000.);
        //			VisBufferAsync vb3(syniter3);
        //
        //			cout << " Two iterators running simultaneously"<<endl;
        //			for (syniter2.originChunks(), syniter3.originChunks();
        //					syniter2.moreChunks() || syniter3.moreChunks();
        //					syniter2.nextChunk(), syniter3.nextChunk()) {
        //				syniter2.origin();
        //				cout <<" Interval 1 start: "<< setprecision(10)<< vb2.time()[0];
        //				Int n=0;
        //				for (;syniter2.more(); syniter2++) {
        //					n++;
        //				}
        //				cout << ", end: "<< setprecision(10)<< vb2.time()[0];
        //				cout << ", # visibilities in interval: "<<n<<endl;
        //				syniter3.origin();
        //				cout <<" Interval 2 start: "<< setprecision(10)<< vb3.time()[0];
        //				n=0;
        //				for (;syniter3.more(); syniter3++) {
        //					n++;
        //				}
        //				cout << ", end: "<< setprecision(10)<< vb3.time()[0];
        //				cout << ", # visibilities in interval: "<<n<<endl;
        //			}
        //		}

        // try the polarization stuff
        // syniter.originChunks(); syniter.origin();
        //cout << " Polarization frame "<< syniter.polFrame()<<endl;
        //cout << " Circular="<< ROVisibilityIterator::Circular<<
        //	", Linear="<< ROVisibilityIterator::Linear<<endl;
        //cout << " CJones(antenna==0)="<< syniter.CJones(0).matrix()<<endl;

        // try assigment
        //VisibilityIterator myiter;
        //myiter=syniter;

        cout << "Exiting scope of tVisibilityIterator" << endl;
    } catch (AipsError x) {
        cout << "Caught exception " << endl;
        cout << x.getMesg() << endl;
        return 1;
    }
    cout << "Done." << endl;
    return 0;
}
