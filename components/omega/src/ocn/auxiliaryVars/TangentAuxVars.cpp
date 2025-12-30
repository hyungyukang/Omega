#include "TangentAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

TangentAuxVars::TangentAuxVars(const std::string &AuxStateSuffix,
                               const HorzMesh *Mesh, const VertCoord *VCoord)
    : TangentialVelocity("TangentialVelocity" + AuxStateSuffix, Mesh->NEdgesSize,
                     VCoord->NVertLayers),
      NEdgesOnEdge(Mesh->NEdgesOnEdge), EdgesOnEdge(Mesh->EdgesOnEdge),
      WeightsOnEdge(Mesh->WeightsOnEdge),
      MinLayerEdgeTop(VCoord->MinLayerEdgeTop),
      MaxLayerEdgeBot(VCoord->MaxLayerEdgeBot) {}

void TangentAuxVars::registerFields(const std::string &AuxGroupName,
                                    const std::string &MeshName) const {

   // Create fields with metadata
   const Real FillValue = -9.99e30;
   int NDims            = 2;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   DimNames[0] = "NEdges" + DimSuffix; // for first three fields
   DimNames[1] = "NVertLayers";        // same for all fields below

   // Tangential velocity on edges
   auto TangentialVelocityField =
       Field::create(TangentialVelocity.label(), // field name
                     "azimuthal component of the horizontal velocity at "
                     "edges",                          // long name/describe
                     "m s^-1",                         // units
                     "",                               // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue, // scalar for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dimension names
       );

   // Add fields to Aux field group
   FieldGroup::addFieldToGroup(TangentialVelocity.label(), AuxGroupName);

   // Attach data to fields
   TangentialVelocityField->attachData<Array2DReal>(TangentialVelocity);
}

void TangentAuxVars::unregisterFields() const {
   Field::destroy(TangentialVelocity.label());
}

} // namespace OMEGA
