//# ImagerTask.cc: 
//# Copyright (C) 2005
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id: ImagerTask.cc,v 1.4 2005/06/16 17:14:23 ddebonis Exp $

#include <synthesis/MeasurementEquations/ImagerTask.h>
#include <casa/Arrays/Vector.h>
#include <ms/MeasurementSets/MeasurementSet.h>

/* CONVERT TO RDO WHEN FINISHED
#include <synthesis/MeasurementEquations/ImagerTool.h>
*/

namespace casa { //# NAMESPACE CASA - BEGIN

ImagerTask::ImagerTask()
 : CASATask(paramsDesc())
{
}

ImagerTask::ImagerTask(const Record &params)
 : CASATask(params)
{
}

ImagerTask::~ImagerTask()
{
}

RecordDesc ImagerTask::paramsDesc()
{
   RecordDesc psetDesc;
   psetDesc.addField("featherimage", TpString);
   psetDesc.addField("vis", TpString);
   psetDesc.addField("alg", TpString);
   psetDesc.addField("niter", TpInt);
   psetDesc.addField("gain", TpDouble);
   psetDesc.addField("threshold", TpDouble);
   psetDesc.addField("residual", TpArrayString);
   psetDesc.addField("image", TpArrayString);
   psetDesc.addField("model", TpArrayString);
   psetDesc.addField("mask", TpArrayString);
   psetDesc.addField("map", TpString);
   psetDesc.addField("beam", TpString);
   psetDesc.addField("mode", TpString);
   psetDesc.addField("gridfn", TpString);
   psetDesc.addField("grid", TpString);
   psetDesc.addField("nchan", TpArrayInt);
   psetDesc.addField("start", TpArrayInt);
   psetDesc.addField("width", TpArrayInt);
   psetDesc.addField("step", TpArrayInt);
   psetDesc.addField("imsize", TpArrayInt);
   psetDesc.addField("cell", TpArrayInt);
   psetDesc.addField("stokes", TpString);
   psetDesc.addField("fieldid", TpArrayInt);
   psetDesc.addField("reffieldid", TpInt);
   psetDesc.addField("spwid", TpArrayInt);
   psetDesc.addField("weighting", TpString);
   psetDesc.addField("rmode", TpString);
   psetDesc.addField("robust", TpDouble);
   psetDesc.addField("highres", TpString);
   psetDesc.addField("lowres", TpString);
   psetDesc.addField("lowpsf", TpString);
   psetDesc.addField("minpb", TpDouble);
   psetDesc.addField("cyclefactor", TpDouble);
   psetDesc.addField("cyclespeedup", TpDouble);
   psetDesc.addField("stoplargenegatives", TpInt);

   return psetDesc;
}

void ImagerTask::clean()
{
   const Record pset(getParams());

   // collect parameters
   String vis("");
   if(pset.isDefined("vis"))
      vis = pset.asString("vis");
   String alg("");
   if(pset.isDefined("alg"))
      alg = pset.asString("alg");
   Int niter(0);
   if(pset.isDefined("niter"))
      niter = pset.asInt("niter");
   Double gain(0);
   if(pset.isDefined("gain"))
      gain = pset.asDouble("gain");
   Double threshold(0);
   if(pset.isDefined("threshold"))
      threshold = pset.asDouble("threshold");
   Vector<String> residual;
   if(pset.isDefined("residual"))
      residual = pset.asArrayString("residual");
   else
      residual = "";
   Vector<String> image;
   if(pset.isDefined("image"))
      image = pset.asArrayString("image");
   else
      image = "";
   Vector<String> model;
   if(pset.isDefined("model"))
      model = pset.asArrayString("model");
   else
      model = "";
   Vector<String> mask;
   if(pset.isDefined("mask"))
      mask = pset.asArrayString("mask");
   else
      mask = "";
   String mode("");
   if(pset.isDefined("mode"))
      mode = pset.asString("mode");
   Vector<Int> nchan;
   if(pset.isDefined("nchan"))
      nchan = pset.asArrayInt("nchan");
   else
      nchan = 0;
   Vector<Int> start;
   if(pset.isDefined("start"))
      start = pset.asArrayInt("start");
   else
      start = 0;
   Vector<Int> width;
   if(pset.isDefined("width"))
      width = pset.asArrayInt("width");
   else
      width = 0;
   Vector<Int> step;
   if(pset.isDefined("step"))
      step = pset.asArrayInt("step");
   else
      step = 0;
   Vector<Int> imsize;
   if(pset.isDefined("imsize"))
      imsize = pset.asArrayInt("imsize");
   else
      imsize = 0;
   Vector<Int> cell;
   if(pset.isDefined("cell"))
      cell = pset.asArrayInt("cell");
   else
      cell = 0;
   String stokes("");
   if(pset.isDefined("stokes"))
      stokes = pset.asString("stokes");
   Vector<Int> fieldid;
   if(pset.isDefined("fieldid"))
      fieldid = pset.asArrayInt("fieldid");
   else
      fieldid = 0;
   Vector<Int> spwid;
   if(pset.isDefined("spwid"))
      spwid = pset.asArrayInt("spwid");
   else
      spwid = 0;
   String weighting("");
   if(pset.isDefined("weighting"))
      weighting = pset.asString("weighting");
   String rmode("");
   if(pset.isDefined("rmode"))
      rmode = pset.asString("rmode");
   Double robust(0);
   if(pset.isDefined("robust"))
      robust = pset.asDouble("robust");

   // create imager and run
   MeasurementSet ms(vis, Table::Update);

/* CONVERT TO RDO WHEN FINISHED
   ImagerTool imgr(ms);

   imgr.setdata(mode, nchan, start, step, Quantity(0, "km/s"),
                Quantity(0, "km/s"), spwid, fieldid);
   imgr.setimage(imsize[0], imsize[1], Quantity(cell[0], "arcsec"),
                 Quantity(cell[1], "arcsec"), stokes, False, MDirection(),
                 Quantity(0, "arcsec"), Quantity(0, "arcsec"), mode, nchan[0],
                 start[0], step[0], Quantity(0, "km/s"), Quantity(0,"km/s"),
                 spwid, fieldid[0], 1, Quantity(0,"m"));
   imgr.weight("uniform", rmode, Quantity(0.0, "Jy"), robust,
               Quantity(0,"rad"), 0);
   imgr.clean(alg, niter, gain, Quantity(threshold, "Jy"), False, model,
              Vector<Bool>(), "", mask, image, residual);
*/
}

void ImagerTask::feather()
{
   const Record pset(getParams());

   // collect parameters
   String vis("");
   if(pset.isDefined("vis"))
      vis = pset.asString("vis");
   String featherimage("");
   if(pset.isDefined("featherimage"))
      featherimage = pset.asString("featherimage");
   String highres("");
   if(pset.isDefined("highres"))
      highres = pset.asString("highres");
   String lowres("");
   if(pset.isDefined("lowres"))
      lowres = pset.asString("lowres");
   String lowpsf("");
   if(pset.isDefined("lowpsf"))
      lowpsf = pset.asString("lowpsf");

   // create imager and run
   MeasurementSet ms(vis, Table::Update);
/* CONVERT TO RDO WHEN FINISHED
   ImagerTool imgr(ms);

   imgr.setvp(True, True, "", False, Quantity(360, "deg"));
   imgr.feather(featherimage, highres, lowres, lowpsf);
*/
}

void ImagerTask::invert()
{
   Record pset(getParams());

   // collect parameters
   String vis("");
   if(pset.isDefined("vis"))
      vis = pset.asString("vis");
   String map("");
   if(pset.isDefined("map"))
      map = pset.asString("map");
   String beam("");
   if(pset.isDefined("beam"))
      beam = pset.asString("beam");
   String mode("");
   if(pset.isDefined("mode"))
      mode = pset.asString("mode");
   Vector<Int> nchan;
   if(pset.isDefined("nchan"))
      nchan = pset.asArrayInt("nchan");
   else
      nchan = 0;
   Vector<Int> start;
   if(pset.isDefined("start"))
      start = pset.asArrayInt("start");
   else
      start = 0;
   Vector<Int> width;
   if(pset.isDefined("width"))
      width = pset.asArrayInt("width");
   else
      width = 0;
   Vector<Int> step;
   if(pset.isDefined("step"))
      step = pset.asArrayInt("step");
   else
      step = 0;
   Vector<Int> imsize;
   if(pset.isDefined("imsize"))
      imsize = pset.asArrayInt("imsize");
   else
      imsize = 0;
   Vector<Int> cell;
   if(pset.isDefined("cell"))
      cell = pset.asArrayInt("cell");
   else
      cell = 0;
   String stokes("");
   if(pset.isDefined("stokes"))
      stokes = pset.asString("stokes");
   Vector<Int> fieldid;
   if(pset.isDefined("fieldid"))
      fieldid = pset.asArrayInt("fieldid");
   else
      fieldid = 0;
   Vector<Int> spwid;
   if(pset.isDefined("spwid"))
      spwid = pset.asArrayInt("spwid");
   else
      spwid = 0;
   String weighting("");
   if(pset.isDefined("weighting"))
      weighting = pset.asString("weighting");
   String rmode("");
   if(pset.isDefined("rmode"))
      rmode = pset.asString("rmode");
   Double robust(0);
   if(pset.isDefined("robust"))
      robust = pset.asDouble("robust");

   // create imager and run
   MeasurementSet ms(vis, Table::Update);
/* CONVERT TO RDO WHEN FINISHED
   ImagerTool imgr(ms);

   imgr.setdata(mode, nchan, start, step, Quantity(0, "km/s"),
                Quantity(0, "km/s"), spwid, fieldid);
   imgr.setimage(imsize[0], imsize[1], Quantity(cell[0], "arcsec"),
                 Quantity(cell[1], "arcsec"), stokes, False, MDirection(),
                 Quantity(0, "arcsec"), Quantity(0, "arcsec"), mode, nchan[0],
                 start[0], step[0], Quantity(0, "km/s"), Quantity(0,"km/s"),
                 spwid, fieldid[0], 1, Quantity(0,"m"));
   imgr.weight("uniform", rmode, Quantity(0.0, "Jy"), robust,
               Quantity(0,"rad"), 0);
   imgr.makeimage("corrected", map);
   imgr.makeimage("psf", beam);
*/
}

void ImagerTask::mosaic()
{
   const Record pset(getParams());

   // collect parameters
   String vis("");
   if(pset.isDefined("vis"))
      vis = pset.asString("vis");
   String alg("");
   if(pset.isDefined("alg"))
      alg = pset.asString("alg");
   Int niter(0);
   if(pset.isDefined("niter"))
      niter = pset.asInt("niter");
   Double gain(0);
   if(pset.isDefined("gain"))
      gain = pset.asDouble("gain");
   Double threshold(0);
   if(pset.isDefined("threshold"))
      threshold = pset.asDouble("threshold");
   Vector<String> residual;
   if(pset.isDefined("residual"))
      residual = pset.asArrayString("residual");
   else
      residual = "";
   Vector<String> image;
   if(pset.isDefined("image"))
      image = pset.asArrayString("image");
   else
      image = "";
   Vector<String> model;
   if(pset.isDefined("model"))
      model = pset.asArrayString("model");
   else
      model = "";
   Vector<String> mask;
   if(pset.isDefined("mask"))
      mask = pset.asArrayString("mask");
   else
      mask = "";
   String mode("");
   if(pset.isDefined("mode"))
      mode = pset.asString("mode");
   String gridfn("");
   if(pset.isDefined("gridfn"))
      gridfn = pset.asString("gridfn");
   String grid("");
   if(pset.isDefined("grid"))
      grid = pset.asString("grid");
   Vector<Int> nchan;
   if(pset.isDefined("nchan"))
      nchan = pset.asArrayInt("nchan");
   else
      nchan = 0;
   Vector<Int> start;
   if(pset.isDefined("start"))
      start = pset.asArrayInt("start");
   else
      start = 0;
   Vector<Int> width;
   if(pset.isDefined("width"))
      width = pset.asArrayInt("width");
   else
      width = 0;
   Vector<Int> step;
   if(pset.isDefined("step"))
      step = pset.asArrayInt("step");
   else
      step = 0;
   Vector<Int> imsize;
   if(pset.isDefined("imsize"))
      imsize = pset.asArrayInt("imsize");
   else
      imsize = 0;
   Vector<Int> cell;
   if(pset.isDefined("cell"))
      cell = pset.asArrayInt("cell");
   else
      cell = 0;
   String stokes("");
   if(pset.isDefined("stokes"))
      stokes = pset.asString("stokes");
   Vector<Int> fieldid;
   if(pset.isDefined("fieldid"))
      fieldid = pset.asArrayInt("fieldid");
   else
      fieldid = 0;
   Int reffieldid(0);
   if(pset.isDefined("reffieldid"))
      reffieldid = pset.asInt("reffieldid");
   Vector<Int> spwid;
   if(pset.isDefined("spwid"))
      spwid = pset.asArrayInt("spwid");
   else
      spwid = 0;
   String weighting("");
   if(pset.isDefined("weighting"))
      weighting = pset.asString("weighting");
   String rmode("");
   if(pset.isDefined("rmode"))
      rmode = pset.asString("rmode");
   Double robust(0);
   if(pset.isDefined("robust"))
      robust = pset.asDouble("robust");
   Double minpb(0);
   if(pset.isDefined("minpb"))
      minpb = pset.asDouble("minpb");
   String scaletype("");
   if(pset.isDefined("scaletype"))
      scaletype = pset.asString("scaletype");
   Double cyclefactor(0);
   if(pset.isDefined("cyclefactor"))
      cyclefactor = pset.asDouble("cyclefactor");
   Double cyclespeedup(0);
   if(pset.isDefined("cyclespeedup"))
      cyclespeedup = pset.asDouble("cyclespeedup");
   Int stoplargenegatives(0);
   if(pset.isDefined("stoplargenegatives"))
      stoplargenegatives = pset.asInt("stoplargenegatives");

   // create imager and run
   MeasurementSet ms(vis, Table::Update);
/* CONVERT TO RDO WHEN FINISHED
   ImagerTool imgr(ms);

   imgr.setvp(True, True, "", False, Quantity(360, "deg"));
   imgr.setdata(mode, nchan, start, step, Quantity(0, "km/s"),
                Quantity(0, "km/s"), spwid, fieldid);
   imgr.setimage(imsize[0], imsize[1], Quantity(cell[0], "arcsec"),
                 Quantity(cell[1], "arcsec"), stokes, False, MDirection(),
                 Quantity(0, "arcsec"), Quantity(0, "arcsec"), mode, nchan[0],
                 start[0], step[0], MRadialVelocity(Quantity(0, "km/s")),
                 MRadialVelocity(Quantity(0,"km/s")),
                 spwid, reffieldid, 1, Quantity(0,"m"));
   imgr.weight(weighting, rmode, Quantity(0.0, "Jy"), robust,
               Quantity(0,"rad"), 0);
   imgr.setoptions(grid, 0, 16, gridfn, MPosition(), 1.2);
   imgr.setmfcontrol(cyclefactor, cyclespeedup, stoplargenegatives, 
		     0, scaletype, minpb, 0.0, Vector<String>());
   imgr.clean(alg, niter, gain, Quantity(threshold, "Jy"), False, model,
              Vector<Bool>(), "", mask, image, residual);
*/
}

} //# NAMESPACE CASA - END
