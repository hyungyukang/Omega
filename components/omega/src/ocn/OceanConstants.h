namespace OMEGA {

   //const I4 NVertLevels = 30;

   const R8 pii = 3.141592653589793;
   const R8 Rearth = 6371220.0;
   //const R8 omega =  7.29212e-5;
   const R8 omega =  7.292e-5;
   const R8 gravity = 9.80616;
   const R8 rgas = 287.0;
   const R8 rv = 461.6;
   const R8 rvord = rv/rgas;
   const R8 cp = 7.0 * rgas / 2.0;
   const R8 cv = cp - rgas;
   const R8 cvpm = -cv/cp;
   const R8 prandtl = 1.0;

   const R8 bottomDragCoeff = 0.0;
   const R8 windStressCoeff = 0.0;
   const R8 mom_del2 = 0.0;
   const R8 mom_del4 = 0.0;

}
